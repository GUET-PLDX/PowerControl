#pragma once

// clang-format off
/* === MODULE MANIFEST V2 ===
module_name: PowerControl
module_description: RM2024 chassis power control with internal energy budget
constructor_args:
  - superpower: '@&superpower'
  - chassis_static_power_loss: 3.5
  - motor_count_3508: 4
  - motor_count_6020: 0
template_args: []
required_hardware: []
depends:
  - pldx/SuperPower
=== END MANIFEST === */
// clang-format on

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

#include "RLS.hpp"
#include "SuperPower.hpp"
#include "app_framework.hpp"
#include "libxr_def.hpp"
#include "mutex.hpp"
#include "pid.hpp"
#include "timebase.hpp"

inline constexpr int POWER_CONTROL_MAX_MOTOR_COUNT = 6;
inline constexpr std::size_t POWER_CONTROL_MAX_TOTAL_MOTOR_COUNT = 12U;
inline constexpr float POWER_CONTROL_PI = static_cast<float>(LibXR::PI);
inline constexpr float POWER_CONTROL_RPM_TO_RAD_PER_SECOND =
    POWER_CONTROL_PI / 30.0f;
inline constexpr float M3508_COMMAND_TO_TORQUE_NM_PER_LSB =
    0.0156224f * 20.0f / 16384.0f;
inline constexpr float GM6020_COMMAND_TO_TORQUE_NM_PER_LSB =
    0.741f * 3.0f / 16384.0f;

/** @brief 功率数据源降级原因。 */
enum class DegradationReason {
  NONE,
  SUPER_CAP_OFFLINE,
  REFEREE_OFFLINE,
  BOTH_OFFLINE,
  INVALID_REFEREE_LIMIT,
};

enum class PowerRequest {
  NORMAL,
  BOOST,
};

/**
 * @brief 按 RM2024 模型预测单个电机功率，不包含底盘静态损耗。
 *
 * P = tau * omega + k1 * abs(omega) + k2 * tau^2
 */
inline float calculate_motor_model_power(
    float command_lsb, float rotor_rpm, float command_to_torque_nm_per_lsb,
    float speed_loss_w_per_rad_per_second,
    float torque_square_loss_w_per_nm_squared) {
  if (!std::isfinite(command_lsb) || !std::isfinite(rotor_rpm) ||
      !std::isfinite(command_to_torque_nm_per_lsb) ||
      !std::isfinite(speed_loss_w_per_rad_per_second) ||
      !std::isfinite(torque_square_loss_w_per_nm_squared) ||
      command_to_torque_nm_per_lsb < 0.0f ||
      speed_loss_w_per_rad_per_second < 0.0f ||
      torque_square_loss_w_per_nm_squared < 0.0f) {
    return std::numeric_limits<float>::max();
  }
  const float TAU = command_lsb * command_to_torque_nm_per_lsb;
  const float OMEGA = rotor_rpm * POWER_CONTROL_RPM_TO_RAD_PER_SECOND;
  const float POWER = TAU * OMEGA +
                      speed_loss_w_per_rad_per_second * std::fabs(OMEGA) +
                      torque_square_loss_w_per_nm_squared * TAU * TAU;
  return std::isfinite(POWER) ? POWER : std::numeric_limits<float>::max();
}

/**
 * @brief 在功率配额内选择最接近原请求的可行电流指令。
 */
inline float solve_current_for_power(float power_quota_w, float rotor_rpm,
                                     float command_to_torque_nm_per_lsb,
                                     float speed_loss_w_per_rad_per_second,
                                     float torque_square_loss_w_per_nm_squared,
                                     float original_command_lsb,
                                     float max_command_magnitude = 16384.0f) {
  constexpr float COEFFICIENT_FLOOR = 1.0e-12f;
  constexpr float DISCRIMINANT_RELATIVE_TOLERANCE = 1.0e-6f;
  constexpr std::size_t MAX_CANDIDATE_COUNT = 7U;

  const float MAX_COMMAND =
      std::isfinite(max_command_magnitude) && max_command_magnitude > 0.0f
          ? max_command_magnitude
          : 16384.0f;
  const float REQUESTED = std::clamp(
      std::isfinite(original_command_lsb) ? original_command_lsb : 0.0f,
      -MAX_COMMAND, MAX_COMMAND);
  const float COMMAND_LOWER_BOUND = std::min(0.0f, REQUESTED);
  const float COMMAND_UPPER_BOUND = std::max(0.0f, REQUESTED);
  const float QUOTA = std::isfinite(power_quota_w) && power_quota_w > 0.0f
                          ? power_quota_w
                          : 0.0f;
  const float REQUESTED_POWER = calculate_motor_model_power(
      REQUESTED, rotor_rpm, command_to_torque_nm_per_lsb,
      speed_loss_w_per_rad_per_second, torque_square_loss_w_per_nm_squared);
  if (REQUESTED_POWER <= QUOTA) {
    return REQUESTED;
  }

  const float OMEGA = rotor_rpm * POWER_CONTROL_RPM_TO_RAD_PER_SECOND;
  const float A = torque_square_loss_w_per_nm_squared *
                  command_to_torque_nm_per_lsb * command_to_torque_nm_per_lsb;
  const float B = command_to_torque_nm_per_lsb * OMEGA;
  const float C = speed_loss_w_per_rad_per_second * std::fabs(OMEGA) - QUOTA;
  std::array<float, MAX_CANDIDATE_COUNT> candidates{};
  std::size_t candidate_count = 0U;
  const auto ADD_CANDIDATE = [&](float command) {
    if (std::isfinite(command) && command >= COMMAND_LOWER_BOUND &&
        command <= COMMAND_UPPER_BOUND && candidate_count < candidates.size()) {
      candidates[candidate_count++] = command;
    }
  };

  ADD_CANDIDATE(REQUESTED);
  ADD_CANDIDATE(0.0f);
  ADD_CANDIDATE(-MAX_COMMAND);
  ADD_CANDIDATE(MAX_COMMAND);
  if (A > COEFFICIENT_FLOOR) {
    ADD_CANDIDATE(-B / (2.0f * A));
    float discriminant = B * B - 4.0f * A * C;
    const float SCALE = std::max(1.0f, B * B + std::fabs(4.0f * A * C));
    const float TOLERANCE = DISCRIMINANT_RELATIVE_TOLERANCE * SCALE;
    if (discriminant >= -TOLERANCE) {
      discriminant = std::max(0.0f, discriminant);
      const float ROOT = std::sqrt(discriminant);
      const float Q = -0.5f * (B + std::copysign(ROOT, B));
      if (std::fabs(Q) > COEFFICIENT_FLOOR) {
        ADD_CANDIDATE(Q / A);
        ADD_CANDIDATE(C / Q);
      } else {
        ADD_CANDIDATE(-B / (2.0f * A));
      }
    }
  } else if (std::fabs(B) > COEFFICIENT_FLOOR) {
    ADD_CANDIDATE(-C / B);
  }

  const float QUOTA_TOLERANCE = std::max(1.0e-4f, QUOTA * 1.0e-5f);
  float best_feasible_command = REQUESTED;
  float best_feasible_distance = std::numeric_limits<float>::max();
  float minimum_power_command = 0.0f;
  float minimum_power = calculate_motor_model_power(
      0.0f, rotor_rpm, command_to_torque_nm_per_lsb,
      speed_loss_w_per_rad_per_second, torque_square_loss_w_per_nm_squared);
  float minimum_command_magnitude = 0.0f;
  bool feasible_found = false;
  for (std::size_t index = 0U; index < candidate_count; ++index) {
    const float CANDIDATE = candidates[index];
    const float POWER = calculate_motor_model_power(
        CANDIDATE, rotor_rpm, command_to_torque_nm_per_lsb,
        speed_loss_w_per_rad_per_second, torque_square_loss_w_per_nm_squared);
    if (!std::isfinite(POWER)) {
      continue;
    }
    const float DISTANCE = std::fabs(CANDIDATE - REQUESTED);
    if (POWER < minimum_power ||
        (std::fabs(POWER - minimum_power) <= QUOTA_TOLERANCE &&
         std::fabs(CANDIDATE) < minimum_command_magnitude)) {
      minimum_power = POWER;
      minimum_power_command = CANDIDATE;
      minimum_command_magnitude = std::fabs(CANDIDATE);
    }
    if (POWER <= QUOTA + QUOTA_TOLERANCE &&
        (!feasible_found || DISTANCE < best_feasible_distance)) {
      feasible_found = true;
      best_feasible_command = CANDIDATE;
      best_feasible_distance = DISTANCE;
    }
  }
  return feasible_found ? best_feasible_command : minimum_power_command;
}

struct PowerControlData {
  float new_output_current_3508[POWER_CONTROL_MAX_MOTOR_COUNT] = {};
  float new_output_current_6020[POWER_CONTROL_MAX_MOTOR_COUNT] = {};
  bool is_power_limited = false;
  bool motor_input_valid = false;
  bool budget_feasible = true;
  float effective_budget_w = 0.0f;
  float requested_predicted_power_w = 0.0f;
  float limited_predicted_power_w = 0.0f;
  float measured_power_w = 0.0f;
  float cap_energy_normalized = 0.0f;
  uint8_t cap_energy_raw = 0U;
  float m3508_speed_loss = 0.22f;
  float m3508_torque_square_loss = 1.2f;
  float gm6020_speed_loss = 0.22f;
  float gm6020_torque_square_loss = 1.2f;
  bool supercap_online = false;
  bool supercap_healthy = false;
  bool referee_power_limit_online = false;
  bool referee_energy_buffer_online = false;
  bool referee_online = false;
  bool rls_updated = false;
  DegradationReason degradation_reason =
      DegradationReason::INVALID_REFEREE_LIMIT;
  PowerRequest requested_mode = PowerRequest::NORMAL;
};

/**
 * @brief RM2024 底盘功率控制模块。
 *
 * 模块内部计算能量预算、在线辨识损耗参数，并为全部 3508/6020 电机执行一次
 * 共享功率分配。所有工作区均为固定成员，不使用堆内存。
 */
class PowerControl : public LibXR::Application {
 public:
  static constexpr int MAX_MOTOR_COUNT = POWER_CONTROL_MAX_MOTOR_COUNT;

  struct AllocationBias3508 {
    bool enabled = false;
    float reserve_fraction = 0.0f;
    float reserve_weight[MAX_MOTOR_COUNT] = {};
    float allocation_weight_scale[MAX_MOTOR_COUNT] = {};
  };

  PowerControl(LibXR::HardwareContainer& hw, LibXR::ApplicationManager& app,
               SuperPower* superpower, float chassis_static_power_loss = 3.5f,
               int motor_count_3508 = 4, int motor_count_6020 = 0)
      : superpower_(superpower),
        base_energy_pid_(MakeEnergyPidParam()),
        full_energy_pid_(MakeEnergyPidParam()),
        rls_(RLS_INITIAL_COVARIANCE, RLS_FORGETTING_FACTOR,
             {DEFAULT_SPEED_LOSS, DEFAULT_TORQUE_SQUARE_LOSS}),
        chassis_static_power_loss_(
            NonnegativeOr(chassis_static_power_loss, 3.5f)),
        motor_count_config_valid_(
            MotorCountConfigValid(motor_count_3508, motor_count_6020)),
        motor_count_3508_(motor_count_config_valid_
                              ? static_cast<std::size_t>(motor_count_3508)
                              : 0U),
        motor_count_6020_(motor_count_config_valid_
                              ? static_cast<std::size_t>(motor_count_6020)
                              : 0U),
        requested_valid_3508_(motor_count_3508_ == 0U),
        requested_valid_6020_(motor_count_6020_ == 0U) {
    UNUSED(hw);
    UNUSED(app);
    rls_.SetParamBounds({RLS_SPEED_LOSS_MIN, RLS_TORQUE_LOSS_MIN},
                        {RLS_SPEED_LOSS_MAX, RLS_TORQUE_LOSS_MAX});
  }

  bool SetMotorData3508(const float* requested_command_lsb,
                        const float* rotor_rpm,
                        const float* tracking_error = nullptr, int count = -1,
                        const bool* active_mask = nullptr) {
    const int CONFIGURED_COUNT = static_cast<int>(motor_count_3508_);
    const int EFFECTIVE_COUNT = count == -1 ? CONFIGURED_COUNT : count;
    MotorInputs samples{};
    const bool INPUT_VALID =
        motor_count_config_valid_ &&
        ValidateMotorData(requested_command_lsb, rotor_rpm, tracking_error,
                          active_mask, EFFECTIVE_COUNT, CONFIGURED_COUNT,
                          samples);
    LibXR::Mutex::LockGuard lock(data_mutex_);
    requested_samples_3508_ = samples;
    requested_valid_3508_ = INPUT_VALID;
    return INPUT_VALID;
  }

  bool SetMotorData6020(const float* requested_command_lsb,
                        const float* rotor_rpm,
                        const float* tracking_error = nullptr, int count = -1,
                        const bool* active_mask = nullptr) {
    const int CONFIGURED_COUNT = static_cast<int>(motor_count_6020_);
    const int EFFECTIVE_COUNT = count == -1 ? CONFIGURED_COUNT : count;
    MotorInputs samples{};
    const bool INPUT_VALID =
        motor_count_config_valid_ &&
        ValidateMotorData(requested_command_lsb, rotor_rpm, tracking_error,
                          active_mask, EFFECTIVE_COUNT, CONFIGURED_COUNT,
                          samples);
    LibXR::Mutex::LockGuard lock(data_mutex_);
    requested_samples_6020_ = samples;
    requested_valid_6020_ = INPUT_VALID;
    return INPUT_VALID;
  }

  bool SetMotorFeedback3508(const float* command_current_lsb,
                            const float* rotor_rpm, int count,
                            const bool* active_mask = nullptr) {
    const int CONFIGURED_COUNT = static_cast<int>(motor_count_3508_);
    MotorInputs samples{};
    const bool INPUT_VALID =
        motor_count_config_valid_ &&
        ValidateMotorData(command_current_lsb, rotor_rpm, nullptr, active_mask,
                          count, CONFIGURED_COUNT, samples);
    LibXR::Mutex::LockGuard lock(data_mutex_);
    feedback_samples_3508_ = samples;
    feedback_ready_3508_ = INPUT_VALID;
    return INPUT_VALID;
  }

  bool SetMotorFeedback6020(const float* command_current_lsb,
                            const float* rotor_rpm, int count,
                            const bool* active_mask = nullptr) {
    const int CONFIGURED_COUNT = static_cast<int>(motor_count_6020_);
    MotorInputs samples{};
    const bool INPUT_VALID =
        motor_count_config_valid_ &&
        ValidateMotorData(command_current_lsb, rotor_rpm, nullptr, active_mask,
                          count, CONFIGURED_COUNT, samples);
    LibXR::Mutex::LockGuard lock(data_mutex_);
    feedback_samples_6020_ = samples;
    feedback_ready_6020_ = INPUT_VALID;
    return INPUT_VALID;
  }

  void SetAllocationBias3508(const AllocationBias3508& bias) {
    AllocationBias3508 sanitized{};
    sanitized.enabled = bias.enabled;
    sanitized.reserve_fraction =
        std::isfinite(bias.reserve_fraction)
            ? std::clamp(bias.reserve_fraction, 0.0f, 1.0f)
            : 0.0f;
    for (int index = 0; index < MAX_MOTOR_COUNT; ++index) {
      const float RESERVE_WEIGHT = bias.reserve_weight[index];
      const float WEIGHT_SCALE = bias.allocation_weight_scale[index];
      sanitized.reserve_weight[index] =
          std::isfinite(RESERVE_WEIGHT) && RESERVE_WEIGHT >= 0.0f
              ? std::min(RESERVE_WEIGHT, MAX_ALLOCATION_WEIGHT)
              : 0.0f;
      sanitized.allocation_weight_scale[index] =
          std::isfinite(WEIGHT_SCALE) && WEIGHT_SCALE > 0.0f
              ? std::min(WEIGHT_SCALE, MAX_ALLOCATION_WEIGHT)
              : 1.0f;
    }
    LibXR::Mutex::LockGuard lock(data_mutex_);
    allocation_bias_3508_ = sanitized;
  }

  void SetPowerRequest(PowerRequest power_request) {
    LibXR::Mutex::LockGuard lock(data_mutex_);
    requested_mode_ = power_request == PowerRequest::BOOST
                          ? PowerRequest::BOOST
                          : PowerRequest::NORMAL;
  }

  void SetBoostRequested(bool boost_requested) {
    SetPowerRequest(boost_requested ? PowerRequest::BOOST
                                    : PowerRequest::NORMAL);
  }

  /** @brief 执行一次完整的同步功率控制周期。 */
  void OutputLimit() {
    if (cycle_active_.test_and_set(std::memory_order_acquire)) {
      return;
    }
    const uint32_t CURRENT_TIMESTAMP_MS =
        static_cast<uint32_t>(LibXR::Timebase::GetMilliseconds());
    const float CYCLE_TIME_S = CalculateCycleTimeSeconds(CURRENT_TIMESTAMP_MS);
    SuperPower::TelemetrySnapshot telemetry{};
    if (superpower_ != nullptr) {
      telemetry = superpower_->GetTelemetrySnapshot();
    }
    SnapshotInputs();
    ProcessCycle(telemetry, CYCLE_TIME_S, CURRENT_TIMESTAMP_MS);
    cycle_active_.clear(std::memory_order_release);
  }

  PowerControlData GetPowerControlData() const {
    LibXR::Mutex::LockGuard lock(data_mutex_);
    return power_control_data_;
  }

  float GetMeasuredPower() const {
    LibXR::Mutex::LockGuard lock(data_mutex_);
    return power_control_data_.measured_power_w;
  }

  float GetCapEnergy() const {
    LibXR::Mutex::LockGuard lock(data_mutex_);
    return power_control_data_.cap_energy_normalized;
  }

  bool IsOnline() const {
    LibXR::Mutex::LockGuard lock(data_mutex_);
    return power_control_data_.supercap_online;
  }

  void OnMonitor() override {}

 private:
  static constexpr float DEFAULT_SPEED_LOSS = 0.22f;
  static constexpr float DEFAULT_TORQUE_SQUARE_LOSS = 1.2f;
  static constexpr float DEFAULT_MAX_COMMAND = 16384.0f;
  static constexpr float MAX_ALLOCATION_WEIGHT = 10.0f;
  static constexpr float MAX_ROTOR_RPM = 100000.0f;
  static constexpr float MAX_TRACKING_ERROR = 1000000.0f;
  static constexpr float POWER_MARGIN_W = 4.0f;
  static constexpr float CAP_FULL_TARGET_RAW = 230.0f;
  static constexpr float CAP_BASE_TARGET_RAW = 30.0f;
  static constexpr float REFEREE_FULL_TARGET_J = 60.0f;
  static constexpr float REFEREE_BASE_TARGET_J = 50.0f;
  static constexpr float MAX_EXTRA_CAP_POWER_W = 300.0f;
  static constexpr float MINIMUM_POWER_FRACTION = 0.8f;
  static constexpr float BOTH_OFFLINE_POWER_FRACTION = 0.85f;
  static constexpr float REFEREE_LIMIT_MIN_W = 45.0f;
  static constexpr float REFEREE_LIMIT_MAX_W = 120.0f;
  static constexpr float CONSERVATIVE_FALLBACK_LIMIT_W = 60.0f;
  static constexpr float ENERGY_PROPORTIONAL_GAIN = 50.0f;
  static constexpr float ENERGY_DERIVATIVE_GAIN = 0.2f;
  static constexpr float DEFAULT_CONTROL_PERIOD_S = 0.001f;
  static constexpr float RLS_INITIAL_COVARIANCE = 1000.0f;
  static constexpr float RLS_FORGETTING_FACTOR = 0.9999f;
  static constexpr float RLS_MIN_SUM_ABS_OMEGA = 20.0f;
  static constexpr float RLS_MIN_SUM_TAU_SQUARED = 0.02f;
  static constexpr float RLS_MIN_INNOVATION_LIMIT_W = 30.0f;
  static constexpr uint16_t RLS_PERSISTENT_INNOVATION_SAMPLES = 32U;
  static constexpr uint32_t RLS_PERSISTENT_INNOVATION_WINDOW_MS = 1000U;
  static constexpr float RLS_MEASURED_POWER_MIN_W = 5.0f;
  static constexpr float RLS_MEASURED_POWER_MAX_W = 1000.0f;
  static constexpr float RLS_RESIDUAL_MIN_W = 0.0f;
  static constexpr float RLS_RESIDUAL_MAX_W = 1000.0f;
  static constexpr float RLS_SPEED_LOSS_MIN = 0.0f;
  static constexpr float RLS_SPEED_LOSS_MAX = 5.0f;
  static constexpr float RLS_TORQUE_LOSS_MIN = 0.0f;
  static constexpr float RLS_TORQUE_LOSS_MAX = 20.0f;
  static constexpr float DENOMINATOR_FLOOR = 1.0e-6f;
  static constexpr float AUDIT_ABSOLUTE_TOLERANCE_W = 0.05f;
  static constexpr float AUDIT_RELATIVE_TOLERANCE = 0.001f;
  static constexpr float AUDIT_CONTRACTION_FACTOR = 0.9f;
  static constexpr std::size_t MAX_AUDIT_ITERATIONS = 8U;

  struct MotorModel {
    float command_to_torque_nm_per_lsb;
    float speed_loss;
    float torque_square_loss;
    float max_command;
  };

  struct MotorInput {
    float command_lsb = 0.0f;
    float rotor_rpm = 0.0f;
    float tracking_error = 0.0f;
    bool active = false;
  };

  struct MotorSample : MotorInput {
    MotorModel model{M3508_COMMAND_TO_TORQUE_NM_PER_LSB, DEFAULT_SPEED_LOSS,
                     DEFAULT_TORQUE_SQUARE_LOSS, DEFAULT_MAX_COMMAND};
    float allocation_weight_scale = 1.0f;
    float reserve_weight = 0.0f;
  };

  using MotorInputs = std::array<MotorInput, MAX_MOTOR_COUNT>;

  struct BudgetStatus {
    float effective_budget_w = 0.0f;
    DegradationReason degradation_reason =
        DegradationReason::INVALID_REFEREE_LIMIT;
  };

  enum class EnergyFeedbackDomain {
    NONE,
    CAP,
    REFEREE,
  };

  struct Workspace {
    MotorInputs requested_3508{};
    MotorInputs requested_6020{};
    MotorInputs feedback_3508{};
    MotorInputs feedback_6020{};
    AllocationBias3508 allocation_bias{};
    std::size_t motor_count_3508 = 0U;
    std::size_t motor_count_6020 = 0U;
    bool motor_count_config_valid = false;
    bool requested_valid_3508 = false;
    bool requested_valid_6020 = false;
    bool feedback_pending = false;
    PowerRequest requested_mode = PowerRequest::NORMAL;
    std::array<MotorSample, POWER_CONTROL_MAX_TOTAL_MOTOR_COUNT> samples{};
    std::array<float, POWER_CONTROL_MAX_TOTAL_MOTOR_COUNT> commands{};
    std::array<float, POWER_CONTROL_MAX_TOTAL_MOTOR_COUNT> requested_power{};
    std::array<float, POWER_CONTROL_MAX_TOTAL_MOTOR_COUNT> quotas{};
    std::array<float, POWER_CONTROL_MAX_TOTAL_MOTOR_COUNT> residual_capacity{};
    std::array<float, POWER_CONTROL_MAX_TOTAL_MOTOR_COUNT> weights{};
    PowerControlData output{};
  };

  static LibXR::PID<float>::Param MakeEnergyPidParam() {
    LibXR::PID<float>::Param param{};
    param.k = 1.0f;
    param.p = ENERGY_PROPORTIONAL_GAIN;
    param.i = 0.0f;
    param.d = ENERGY_DERIVATIVE_GAIN;
    param.i_limit = 0.0f;
    param.out_limit = MAX_EXTRA_CAP_POWER_W;
    param.cycle = false;
    return param;
  }

  static float NonnegativeOr(float value, float fallback) {
    return std::isfinite(value) && value >= 0.0f ? value : fallback;
  }

  static bool MotorCountConfigValid(int motor_count_3508,
                                    int motor_count_6020) {
    return motor_count_3508 >= 0 && motor_count_3508 <= MAX_MOTOR_COUNT &&
           motor_count_6020 >= 0 && motor_count_6020 <= MAX_MOTOR_COUNT;
  }

  float CalculateCycleTimeSeconds(uint32_t current_timestamp_ms) {
    float cycle_time_s = DEFAULT_CONTROL_PERIOD_S;
    if (control_timestamp_initialized_) {
      const uint32_t ELAPSED_MS =
          current_timestamp_ms - last_control_timestamp_ms_;
      if (ELAPSED_MS > 0U) {
        cycle_time_s = static_cast<float>(ELAPSED_MS) / 1000.0f;
      }
    }
    last_control_timestamp_ms_ = current_timestamp_ms;
    control_timestamp_initialized_ = true;
    return std::isfinite(cycle_time_s) && cycle_time_s > 0.0f
               ? cycle_time_s
               : DEFAULT_CONTROL_PERIOD_S;
  }

  static bool ValidateMotorData(const float* command_lsb,
                                const float* rotor_rpm,
                                const float* tracking_error,
                                const bool* active_mask, int count,
                                int configured_count, MotorInputs& samples) {
    samples = {};
    if (configured_count < 0 || configured_count > MAX_MOTOR_COUNT ||
        count != configured_count) {
      return false;
    }
    if (configured_count == 0) {
      return true;
    }
    if (command_lsb == nullptr || rotor_rpm == nullptr) {
      return false;
    }

    for (int index = 0; index < configured_count; ++index) {
      MotorInput& sample = samples[static_cast<std::size_t>(index)];
      sample.active = active_mask == nullptr || active_mask[index];
      if (!sample.active) {
        continue;
      }
      const bool COMMAND_VALID = std::isfinite(command_lsb[index]);
      const bool RPM_VALID = std::isfinite(rotor_rpm[index]) &&
                             std::fabs(rotor_rpm[index]) <= MAX_ROTOR_RPM;
      const bool ERROR_VALID =
          tracking_error == nullptr ||
          (std::isfinite(tracking_error[index]) &&
           std::fabs(tracking_error[index]) <= MAX_TRACKING_ERROR);
      if (!COMMAND_VALID || !RPM_VALID || !ERROR_VALID) {
        samples = {};
        return false;
      }
      sample.command_lsb = std::clamp(command_lsb[index], -DEFAULT_MAX_COMMAND,
                                      DEFAULT_MAX_COMMAND);
      sample.rotor_rpm = rotor_rpm[index];
      sample.tracking_error =
          tracking_error == nullptr ? 0.0f : std::fabs(tracking_error[index]);
    }
    return true;
  }

  MotorModel Make3508Model() const {
    const RLS<2>::ParamVector& params = rls_.GetParamVector();
    return {M3508_COMMAND_TO_TORQUE_NM_PER_LSB, params[0], params[1],
            DEFAULT_MAX_COMMAND};
  }

  static MotorModel Make6020Model() {
    return {GM6020_COMMAND_TO_TORQUE_NM_PER_LSB, DEFAULT_SPEED_LOSS,
            DEFAULT_TORQUE_SQUARE_LOSS, DEFAULT_MAX_COMMAND};
  }

  static MotorSample MakeMotorSample(const MotorInput& input,
                                     const MotorModel& model) {
    MotorSample sample{};
    sample.command_lsb = input.command_lsb;
    sample.rotor_rpm = input.rotor_rpm;
    sample.tracking_error = input.tracking_error;
    sample.active = input.active;
    sample.model = model;
    return sample;
  }

  static bool AllActive(const MotorInputs& samples, std::size_t count) {
    for (std::size_t index = 0U; index < count; ++index) {
      if (!samples[index].active) {
        return false;
      }
    }
    return true;
  }

  static float PredictMotorPower(const MotorSample& sample, float command_lsb) {
    if (!sample.active) {
      return 0.0f;
    }
    return calculate_motor_model_power(
        std::clamp(command_lsb, -sample.model.max_command,
                   sample.model.max_command),
        sample.rotor_rpm, sample.model.command_to_torque_nm_per_lsb,
        sample.model.speed_loss, sample.model.torque_square_loss);
  }

  static float SolveCommand(const MotorSample& sample, float quota_w) {
    return solve_current_for_power(
        quota_w, sample.rotor_rpm, sample.model.command_to_torque_nm_per_lsb,
        sample.model.speed_loss, sample.model.torque_square_loss,
        sample.command_lsb, sample.model.max_command);
  }

  void SnapshotInputs() {
    LibXR::Mutex::LockGuard lock(data_mutex_);
    workspace_.requested_3508 = requested_samples_3508_;
    workspace_.requested_6020 = requested_samples_6020_;
    workspace_.feedback_3508 = feedback_samples_3508_;
    workspace_.feedback_6020 = feedback_samples_6020_;
    workspace_.allocation_bias = allocation_bias_3508_;
    workspace_.motor_count_3508 = motor_count_3508_;
    workspace_.motor_count_6020 = motor_count_6020_;
    workspace_.motor_count_config_valid = motor_count_config_valid_;
    workspace_.requested_valid_3508 = requested_valid_3508_;
    workspace_.requested_valid_6020 = requested_valid_6020_;
    workspace_.feedback_pending =
        feedback_ready_3508_ &&
        (motor_count_6020_ == 0U || feedback_ready_6020_);
    workspace_.requested_mode = requested_mode_;
    feedback_ready_3508_ = false;
    feedback_ready_6020_ = false;
  }

  bool UpdatePowerModel(const SuperPower::TelemetrySnapshot& telemetry,
                        uint32_t current_timestamp_ms) {
    ExpireInnovationRecovery(current_timestamp_ms);
    const bool CAP_USABLE =
        telemetry.supercap_online && telemetry.supercap_healthy;
    const bool NEW_POWER_SAMPLE =
        !rls_power_sample_consumed_ ||
        telemetry.chassis_power_sequence != last_rls_power_sample_sequence_;
    if (!CAP_USABLE || !std::isfinite(telemetry.chassis_power_w)) {
      return false;
    }
    if (!workspace_.feedback_pending || !NEW_POWER_SAMPLE) {
      return false;
    }

    last_rls_power_sample_sequence_ = telemetry.chassis_power_sequence;
    rls_power_sample_consumed_ = true;
    if (!AllActive(workspace_.feedback_3508, workspace_.motor_count_3508) ||
        !AllActive(workspace_.feedback_6020, workspace_.motor_count_6020)) {
      return false;
    }
    if (telemetry.chassis_power_w < RLS_MEASURED_POWER_MIN_W ||
        telemetry.chassis_power_w > RLS_MEASURED_POWER_MAX_W) {
      return false;
    }

    const MotorModel M3508_MODEL = Make3508Model();
    const MotorModel GM6020_MODEL = Make6020Model();
    float sum_abs_omega = 0.0f;
    float sum_tau_squared = 0.0f;
    float mechanical_power_w = 0.0f;
    float fixed_group_power_w = 0.0f;
    for (std::size_t index = 0U; index < workspace_.motor_count_3508; ++index) {
      const MotorInput& sample = workspace_.feedback_3508[index];
      if (!sample.active) {
        continue;
      }
      const float TAU =
          sample.command_lsb * M3508_MODEL.command_to_torque_nm_per_lsb;
      const float OMEGA =
          sample.rotor_rpm * POWER_CONTROL_RPM_TO_RAD_PER_SECOND;
      sum_abs_omega += std::fabs(OMEGA);
      sum_tau_squared += TAU * TAU;
      mechanical_power_w += TAU * OMEGA;
    }
    for (std::size_t index = 0U; index < workspace_.motor_count_6020; ++index) {
      MotorSample sample =
          MakeMotorSample(workspace_.feedback_6020[index], GM6020_MODEL);
      fixed_group_power_w += PredictMotorPower(sample, sample.command_lsb);
    }

    const float RESIDUAL = telemetry.chassis_power_w - mechanical_power_w -
                           chassis_static_power_loss_ - fixed_group_power_w;
    const bool SUFFICIENT_EXCITATION =
        sum_abs_omega >= RLS_MIN_SUM_ABS_OMEGA ||
        sum_tau_squared >= RLS_MIN_SUM_TAU_SQUARED;
    if (!SUFFICIENT_EXCITATION || !std::isfinite(RESIDUAL) ||
        RESIDUAL < RLS_RESIDUAL_MIN_W || RESIDUAL > RLS_RESIDUAL_MAX_W) {
      return false;
    }
    const RLS<2>::ParamVector& PARAMS = rls_.GetParamVector();
    const float PREDICTED_RESIDUAL =
        PARAMS[0] * sum_abs_omega + PARAMS[1] * sum_tau_squared;
    const float INNOVATION = RESIDUAL - PREDICTED_RESIDUAL;
    const float INNOVATION_LIMIT =
        std::max(RLS_MIN_INNOVATION_LIMIT_W, std::fabs(PREDICTED_RESIDUAL));
    if (!std::isfinite(INNOVATION)) {
      ResetInnovationRecovery();
      return false;
    }

    const bool LARGE_INNOVATION = std::fabs(INNOVATION) > INNOVATION_LIMIT;
    if (LARGE_INNOVATION) {
      if (INNOVATION < 0.0f) {
        ResetInnovationRecovery();
        return false;
      }
      last_persistent_innovation_timestamp_ms_ = current_timestamp_ms;
      persistent_innovation_timestamp_initialized_ = true;
      if (!innovation_recovery_active_) {
        if (persistent_innovation_count_ < RLS_PERSISTENT_INNOVATION_SAMPLES) {
          ++persistent_innovation_count_;
        }
        if (persistent_innovation_count_ < RLS_PERSISTENT_INNOVATION_SAMPLES) {
          return false;
        }
        rls_.Reset(PARAMS);
        innovation_recovery_active_ = true;
      }
    } else {
      ResetInnovationRecovery();
    }
    return rls_.Update({sum_abs_omega, sum_tau_squared}, RESIDUAL);
  }

  void ResetInnovationRecovery() {
    persistent_innovation_count_ = 0U;
    innovation_recovery_active_ = false;
    persistent_innovation_timestamp_initialized_ = false;
  }

  void ExpireInnovationRecovery(uint32_t current_timestamp_ms) {
    if (persistent_innovation_timestamp_initialized_ &&
        current_timestamp_ms - last_persistent_innovation_timestamp_ms_ >
            RLS_PERSISTENT_INNOVATION_WINDOW_MS) {
      ResetInnovationRecovery();
    }
  }

  float ComputeEnergyBound(LibXR::PID<float>& energy_pid, float referee_limit_w,
                           float target, float feedback, float minimum_power_w,
                           float cycle_time_s, bool initialize) {
    const float TARGET = std::max(0.0f, target);
    const float FEEDBACK = std::max(0.0f, feedback);
    const float TARGET_ROOT = std::sqrt(TARGET);
    const float FEEDBACK_ROOT = std::sqrt(FEEDBACK);
    const float CORRECTION =
        initialize
            ? energy_pid.Calculate(TARGET_ROOT, FEEDBACK_ROOT, 0.0f,
                                   cycle_time_s)
            : energy_pid.Calculate(TARGET_ROOT, FEEDBACK_ROOT, cycle_time_s);
    const float BOUND = referee_limit_w - CORRECTION;
    return std::max(minimum_power_w,
                    std::isfinite(BOUND) ? BOUND : minimum_power_w);
  }

  static DegradationReason DetermineDegradation(
      bool supercap_online, bool referee_power_limit_online,
      bool referee_energy_buffer_online, bool used_fallback_limit) {
    if (used_fallback_limit) {
      return DegradationReason::INVALID_REFEREE_LIMIT;
    }
    if (!supercap_online && !referee_energy_buffer_online) {
      return DegradationReason::BOTH_OFFLINE;
    }
    if (!supercap_online) {
      return DegradationReason::SUPER_CAP_OFFLINE;
    }
    if (!referee_power_limit_online || !referee_energy_buffer_online) {
      return DegradationReason::REFEREE_OFFLINE;
    }
    return DegradationReason::NONE;
  }

  BudgetStatus CalculatePowerBudget(
      const SuperPower::TelemetrySnapshot& telemetry, float cycle_time_s) {
    const bool CAP_USABLE =
        telemetry.supercap_online && telemetry.supercap_healthy;
    const bool REFEREE_POWER_LIMIT_ONLINE =
        telemetry.referee_power_limit_online;
    const bool REFEREE_ENERGY_BUFFER_ONLINE =
        telemetry.referee_energy_buffer_online;
    const EnergyFeedbackDomain ENERGY_DOMAIN =
        CAP_USABLE
            ? EnergyFeedbackDomain::CAP
            : (REFEREE_ENERGY_BUFFER_ONLINE ? EnergyFeedbackDomain::REFEREE
                                            : EnergyFeedbackDomain::NONE);
    const float REFEREE_LIMIT =
        static_cast<float>(telemetry.referee_power_limit_w);
    const bool CURRENT_LIMIT_VALID = REFEREE_POWER_LIMIT_ONLINE &&
                                     std::isfinite(REFEREE_LIMIT) &&
                                     REFEREE_LIMIT >= REFEREE_LIMIT_MIN_W &&
                                     REFEREE_LIMIT <= REFEREE_LIMIT_MAX_W;
    const bool INVALID_REFEREE_LIMIT =
        REFEREE_POWER_LIMIT_ONLINE && !CURRENT_LIMIT_VALID;
    if (CURRENT_LIMIT_VALID) {
      latest_referee_limit_w_ = REFEREE_LIMIT;
      has_valid_referee_limit_ = true;
    }

    BudgetStatus status{};
    const float REFEREE_LIMIT_TO_USE =
        !has_valid_referee_limit_
            ? REFEREE_LIMIT_MIN_W
            : (INVALID_REFEREE_LIMIT ? std::min(latest_referee_limit_w_,
                                                CONSERVATIVE_FALLBACK_LIMIT_W)
                                     : latest_referee_limit_w_);
    const bool EXTRA_POWER_ALLOWED =
        has_valid_referee_limit_ && !INVALID_REFEREE_LIMIT;
    const float MINIMUM_CONFIGURED_POWER =
        REFEREE_LIMIT_TO_USE * MINIMUM_POWER_FRACTION;
    status.degradation_reason = DetermineDegradation(
        CAP_USABLE, REFEREE_POWER_LIMIT_ONLINE, REFEREE_ENERGY_BUFFER_ONLINE,
        INVALID_REFEREE_LIMIT);

    float effective_budget = 0.0f;

    if (!CAP_USABLE && !REFEREE_ENERGY_BUFFER_ONLINE) {
      base_energy_pid_.Reset();
      full_energy_pid_.Reset();
      energy_feedback_domain_ = EnergyFeedbackDomain::NONE;
      const float CONSERVATIVE_BUDGET =
          REFEREE_LIMIT_TO_USE * BOTH_OFFLINE_POWER_FRACTION;
      effective_budget = CONSERVATIVE_BUDGET;
    } else {
      float base_target = 0.0f;
      float full_target = 0.0f;
      float energy_feedback = 0.0f;
      float upper_request = 0.0f;
      if (CAP_USABLE) {
        energy_feedback = std::clamp(
            static_cast<float>(telemetry.cap_energy_raw), 0.0f, 255.0f);
        base_target = CAP_BASE_TARGET_RAW;
        full_target = CAP_FULL_TARGET_RAW;
        upper_request = REFEREE_LIMIT_TO_USE +
                        (EXTRA_POWER_ALLOWED ? MAX_EXTRA_CAP_POWER_W : 0.0f);
      } else {
        energy_feedback =
            std::clamp(static_cast<float>(telemetry.referee_energy_buffer_j),
                       0.0f, REFEREE_FULL_TARGET_J);
        base_target = REFEREE_BASE_TARGET_J;
        full_target = REFEREE_FULL_TARGET_J;
        const float EXTRA_REFEREE_BUFFER_POWER =
            ENERGY_PROPORTIONAL_GAIN * (std::sqrt(REFEREE_FULL_TARGET_J) -
                                        std::sqrt(REFEREE_BASE_TARGET_J));
        upper_request =
            REFEREE_LIMIT_TO_USE +
            (EXTRA_POWER_ALLOWED ? EXTRA_REFEREE_BUFFER_POWER : 0.0f);
      }
      const bool INITIALIZE_ENERGY_CONTROLLERS =
          ENERGY_DOMAIN != energy_feedback_domain_;
      float base_bound = ComputeEnergyBound(
          base_energy_pid_, REFEREE_LIMIT_TO_USE, base_target, energy_feedback,
          MINIMUM_CONFIGURED_POWER, cycle_time_s,
          INITIALIZE_ENERGY_CONTROLLERS);
      float full_bound = ComputeEnergyBound(
          full_energy_pid_, REFEREE_LIMIT_TO_USE, full_target, energy_feedback,
          MINIMUM_CONFIGURED_POWER, cycle_time_s,
          INITIALIZE_ENERGY_CONTROLLERS);
      energy_feedback_domain_ = ENERGY_DOMAIN;
      if (full_bound > base_bound) {
        std::swap(full_bound, base_bound);
      }
      const float REQUESTED_POWER =
          workspace_.requested_mode == PowerRequest::BOOST
              ? upper_request
              : REFEREE_LIMIT_TO_USE;
      const float ENERGY_CONSTRAINED_BUDGET =
          std::clamp(REQUESTED_POWER, full_bound, base_bound);
      const bool EXTRA_POWER_REQUESTED =
          workspace_.requested_mode == PowerRequest::BOOST;
      effective_budget =
          EXTRA_POWER_ALLOWED && EXTRA_POWER_REQUESTED
              ? ENERGY_CONSTRAINED_BUDGET
              : std::min(ENERGY_CONSTRAINED_BUDGET, REFEREE_LIMIT_TO_USE);
    }

    effective_budget =
        NonnegativeOr(effective_budget, CONSERVATIVE_FALLBACK_LIMIT_W);
    status.effective_budget_w = effective_budget;
    return status;
  }

  void PrepareMotorSamples() {
    workspace_.samples = {};
    const MotorModel M3508_MODEL = Make3508Model();
    const MotorModel GM6020_MODEL = Make6020Model();
    for (std::size_t index = 0U; index < workspace_.motor_count_3508; ++index) {
      workspace_.samples[index] =
          MakeMotorSample(workspace_.requested_3508[index], M3508_MODEL);
      if (workspace_.allocation_bias.enabled) {
        workspace_.samples[index].reserve_weight =
            workspace_.allocation_bias.reserve_weight[index];
        workspace_.samples[index].allocation_weight_scale =
            workspace_.allocation_bias.allocation_weight_scale[index];
      }
    }
    for (std::size_t index = 0U; index < workspace_.motor_count_6020; ++index) {
      const std::size_t OUTPUT_INDEX = workspace_.motor_count_3508 + index;
      workspace_.samples[OUTPUT_INDEX] =
          MakeMotorSample(workspace_.requested_6020[index], GM6020_MODEL);
    }
  }

  void UpdateLimitedPrediction(std::size_t motor_count, float static_loss_w) {
    float limited_total = static_loss_w;
    for (std::size_t index = 0U; index < motor_count; ++index) {
      if (!workspace_.samples[index].active) {
        continue;
      }
      workspace_.commands[index] =
          std::clamp(workspace_.commands[index],
                     -workspace_.samples[index].model.max_command,
                     workspace_.samples[index].model.max_command);
      limited_total += PredictMotorPower(workspace_.samples[index],
                                         workspace_.commands[index]);
    }
    workspace_.output.limited_predicted_power_w = limited_total;
  }

  void AllocatePower(float effective_budget_w) {
    const std::size_t MOTOR_COUNT =
        workspace_.motor_count_3508 + workspace_.motor_count_6020;
    const float BUDGET = NonnegativeOr(effective_budget_w, 0.0f);
    const float STATIC_LOSS = chassis_static_power_loss_;
    const float RESERVE_FRACTION =
        workspace_.allocation_bias.enabled
            ? std::clamp(workspace_.allocation_bias.reserve_fraction, 0.0f,
                         1.0f)
            : 0.0f;
    const float AUDIT_TOLERANCE =
        std::max(AUDIT_ABSOLUTE_TOLERANCE_W, BUDGET * AUDIT_RELATIVE_TOLERANCE);

    workspace_.commands = {};
    workspace_.requested_power = {};
    workspace_.quotas = {};
    workspace_.residual_capacity = {};
    workspace_.weights = {};
    float requested_total = STATIC_LOSS;
    float positive_required = 0.0f;
    float negative_power = 0.0f;
    float sum_error = 0.0f;
    for (std::size_t index = 0U; index < MOTOR_COUNT; ++index) {
      const MotorSample& sample = workspace_.samples[index];
      if (!sample.active) {
        continue;
      }
      workspace_.commands[index] = sample.command_lsb;
      const float POWER = PredictMotorPower(sample, sample.command_lsb);
      workspace_.requested_power[index] = POWER;
      requested_total += POWER;
      if (POWER > 0.0f) {
        positive_required += POWER;
        sum_error += sample.tracking_error;
      } else {
        negative_power += POWER;
      }
    }

    workspace_.output.effective_budget_w = BUDGET;
    workspace_.output.budget_feasible = true;
    workspace_.output.requested_predicted_power_w = requested_total;
    const float POSITIVE_POOL =
        std::max(0.0f, BUDGET - STATIC_LOSS - negative_power);
    const bool REQUIRES_LIMIT =
        requested_total > BUDGET + AUDIT_TOLERANCE &&
        positive_required > POSITIVE_POOL + DENOMINATOR_FLOOR;
    if (REQUIRES_LIMIT) {
      workspace_.output.is_power_limited = true;
      float reserve_weight_sum = 0.0f;
      for (std::size_t index = 0U; index < MOTOR_COUNT; ++index) {
        if (workspace_.requested_power[index] > 0.0f) {
          reserve_weight_sum += workspace_.samples[index].reserve_weight;
        }
      }
      const float RESERVE_POOL =
          std::min(POSITIVE_POOL, positive_required) * RESERVE_FRACTION;
      float reserved_sum = 0.0f;
      if (RESERVE_POOL > DENOMINATOR_FLOOR &&
          reserve_weight_sum > DENOMINATOR_FLOOR) {
        for (std::size_t index = 0U; index < MOTOR_COUNT; ++index) {
          if (workspace_.requested_power[index] <= 0.0f) {
            continue;
          }
          const float RESERVE_QUOTA = RESERVE_POOL *
                                      workspace_.samples[index].reserve_weight /
                                      reserve_weight_sum;
          workspace_.quotas[index] =
              std::min(workspace_.requested_power[index], RESERVE_QUOTA);
          reserved_sum += workspace_.quotas[index];
        }
      }

      float residual_required = 0.0f;
      for (std::size_t index = 0U; index < MOTOR_COUNT; ++index) {
        if (workspace_.requested_power[index] <= 0.0f) {
          continue;
        }
        workspace_.residual_capacity[index] = std::max(
            0.0f, workspace_.requested_power[index] - workspace_.quotas[index]);
        residual_required += workspace_.residual_capacity[index];
      }

      float error_confidence = 0.0f;
      if (sum_error > 20.0f) {
        error_confidence = 1.0f;
      } else if (sum_error > 15.0f) {
        error_confidence = (sum_error - 15.0f) / 5.0f;
      }
      for (std::size_t index = 0U; index < MOTOR_COUNT; ++index) {
        if (workspace_.requested_power[index] <= 0.0f ||
            workspace_.residual_capacity[index] <= DENOMINATOR_FLOOR) {
          continue;
        }
        const float ERROR_SHARE =
            sum_error > DENOMINATOR_FLOOR
                ? workspace_.samples[index].tracking_error / sum_error
                : 0.0f;
        const float PROPORTIONAL_SHARE =
            residual_required > DENOMINATOR_FLOOR
                ? workspace_.residual_capacity[index] / residual_required
                : 0.0f;
        workspace_.weights[index] =
            (error_confidence * ERROR_SHARE +
             (1.0f - error_confidence) * PROPORTIONAL_SHARE) *
            workspace_.samples[index].allocation_weight_scale;
      }

      float remaining_pool = std::max(0.0f, POSITIVE_POOL - reserved_sum);
      for (std::size_t pass = 0U;
           pass < MOTOR_COUNT && remaining_pool > DENOMINATOR_FLOOR; ++pass) {
        float weight_sum = 0.0f;
        float capacity_sum = 0.0f;
        for (std::size_t index = 0U; index < MOTOR_COUNT; ++index) {
          if (workspace_.residual_capacity[index] > DENOMINATOR_FLOOR) {
            weight_sum += workspace_.weights[index];
            capacity_sum += workspace_.residual_capacity[index];
          }
        }
        if (capacity_sum <= DENOMINATOR_FLOOR) {
          break;
        }

        float used_power = 0.0f;
        for (std::size_t index = 0U; index < MOTOR_COUNT; ++index) {
          if (workspace_.residual_capacity[index] <= DENOMINATOR_FLOOR) {
            continue;
          }
          const float SHARE =
              weight_sum > DENOMINATOR_FLOOR
                  ? remaining_pool * workspace_.weights[index] / weight_sum
                  : remaining_pool * workspace_.residual_capacity[index] /
                        capacity_sum;
          const float ADDITION =
              std::min(workspace_.residual_capacity[index], SHARE);
          workspace_.quotas[index] += ADDITION;
          workspace_.residual_capacity[index] -= ADDITION;
          used_power += ADDITION;
        }
        if (used_power <= DENOMINATOR_FLOOR) {
          break;
        }
        remaining_pool = std::max(0.0f, remaining_pool - used_power);
      }

      for (std::size_t index = 0U; index < MOTOR_COUNT; ++index) {
        if (workspace_.requested_power[index] > 0.0f) {
          const float QUOTA = std::min(workspace_.requested_power[index],
                                       workspace_.quotas[index]);
          workspace_.commands[index] =
              SolveCommand(workspace_.samples[index], QUOTA);
        }
      }
    }

    UpdateLimitedPrediction(MOTOR_COUNT, STATIC_LOSS);
    for (std::size_t iteration = 0U;
         iteration < MAX_AUDIT_ITERATIONS &&
         workspace_.output.limited_predicted_power_w > BUDGET + AUDIT_TOLERANCE;
         ++iteration) {
      const float PREVIOUS_POWER = workspace_.output.limited_predicted_power_w;
      for (std::size_t index = 0U; index < MOTOR_COUNT; ++index) {
        if (workspace_.requested_power[index] > 0.0f) {
          workspace_.quotas[index] = std::max(
              0.0f, workspace_.quotas[index] * AUDIT_CONTRACTION_FACTOR);
          workspace_.commands[index] =
              SolveCommand(workspace_.samples[index], workspace_.quotas[index]);
        }
      }
      UpdateLimitedPrediction(MOTOR_COUNT, STATIC_LOSS);
      workspace_.output.is_power_limited = true;
      if (workspace_.output.limited_predicted_power_w >=
          PREVIOUS_POWER - DENOMINATOR_FLOOR) {
        break;
      }
    }
    if (workspace_.output.limited_predicted_power_w >
        BUDGET + AUDIT_TOLERANCE) {
      workspace_.commands = {};
      UpdateLimitedPrediction(MOTOR_COUNT, STATIC_LOSS);
      workspace_.output.is_power_limited = true;
      workspace_.output.budget_feasible =
          workspace_.output.limited_predicted_power_w <=
          BUDGET + AUDIT_TOLERANCE;
    }
  }

  void ProcessCycle(const SuperPower::TelemetrySnapshot& telemetry,
                    float cycle_time_s, uint32_t current_timestamp_ms) {
    const bool RLS_UPDATED = UpdatePowerModel(telemetry, current_timestamp_ms);
    const BudgetStatus BUDGET_STATUS =
        CalculatePowerBudget(telemetry, cycle_time_s);
    const float EFFECTIVE_BUDGET =
        std::max(0.0f, BUDGET_STATUS.effective_budget_w - POWER_MARGIN_W);
    PrepareMotorSamples();
    workspace_.output = {};
    AllocatePower(EFFECTIVE_BUDGET);

    const bool REQUEST_INPUT_VALID = workspace_.motor_count_config_valid &&
                                     workspace_.requested_valid_3508 &&
                                     workspace_.requested_valid_6020;
    if (REQUEST_INPUT_VALID) {
      for (std::size_t index = 0U; index < workspace_.motor_count_3508;
           ++index) {
        workspace_.output.new_output_current_3508[index] =
            workspace_.samples[index].active ? workspace_.commands[index]
                                             : 0.0f;
      }
      for (std::size_t index = 0U; index < workspace_.motor_count_6020;
           ++index) {
        const std::size_t OUTPUT_INDEX = workspace_.motor_count_3508 + index;
        workspace_.output.new_output_current_6020[index] =
            workspace_.samples[OUTPUT_INDEX].active
                ? workspace_.commands[OUTPUT_INDEX]
                : 0.0f;
      }
    } else {
      workspace_.output.is_power_limited = true;
    }

    const RLS<2>::ParamVector& params = rls_.GetParamVector();
    workspace_.output.motor_input_valid = REQUEST_INPUT_VALID;
    workspace_.output.measured_power_w =
        std::isfinite(telemetry.chassis_power_w) ? telemetry.chassis_power_w
                                                 : 0.0f;
    workspace_.output.cap_energy_normalized =
        std::isfinite(telemetry.cap_energy_normalized)
            ? std::clamp(telemetry.cap_energy_normalized, 0.0f, 1.0f)
            : 0.0f;
    workspace_.output.cap_energy_raw = telemetry.cap_energy_raw;
    workspace_.output.m3508_speed_loss = params[0];
    workspace_.output.m3508_torque_square_loss = params[1];
    workspace_.output.gm6020_speed_loss = DEFAULT_SPEED_LOSS;
    workspace_.output.gm6020_torque_square_loss = DEFAULT_TORQUE_SQUARE_LOSS;
    workspace_.output.supercap_online = telemetry.supercap_online;
    workspace_.output.supercap_healthy = telemetry.supercap_healthy;
    workspace_.output.referee_power_limit_online =
        telemetry.referee_power_limit_online;
    workspace_.output.referee_energy_buffer_online =
        telemetry.referee_energy_buffer_online;
    workspace_.output.referee_online = telemetry.referee_power_limit_online &&
                                       telemetry.referee_energy_buffer_online;
    workspace_.output.rls_updated = RLS_UPDATED;
    workspace_.output.degradation_reason = BUDGET_STATUS.degradation_reason;
    workspace_.output.requested_mode = workspace_.requested_mode;

    LibXR::Mutex::LockGuard lock(data_mutex_);
    power_control_data_ = workspace_.output;
  }

  mutable LibXR::Mutex data_mutex_;
  std::atomic_flag cycle_active_ = ATOMIC_FLAG_INIT;
  SuperPower* superpower_;
  LibXR::PID<float> base_energy_pid_;
  LibXR::PID<float> full_energy_pid_;
  RLS<2> rls_;
  float chassis_static_power_loss_ = 0.0f;
  bool motor_count_config_valid_ = false;
  std::size_t motor_count_3508_ = 0U;
  std::size_t motor_count_6020_ = 0U;
  bool requested_valid_3508_ = false;
  bool requested_valid_6020_ = false;
  MotorInputs requested_samples_3508_{};
  MotorInputs requested_samples_6020_{};
  MotorInputs feedback_samples_3508_{};
  MotorInputs feedback_samples_6020_{};
  bool feedback_ready_3508_ = false;
  bool feedback_ready_6020_ = false;
  AllocationBias3508 allocation_bias_3508_{};
  PowerRequest requested_mode_ = PowerRequest::NORMAL;
  uint32_t last_rls_power_sample_sequence_ = 0U;
  bool rls_power_sample_consumed_ = false;
  uint16_t persistent_innovation_count_ = 0U;
  bool innovation_recovery_active_ = false;
  uint32_t last_persistent_innovation_timestamp_ms_ = 0U;
  bool persistent_innovation_timestamp_initialized_ = false;
  float latest_referee_limit_w_ = 0.0f;
  bool has_valid_referee_limit_ = false;
  uint32_t last_control_timestamp_ms_ = 0U;
  bool control_timestamp_initialized_ = false;
  EnergyFeedbackDomain energy_feedback_domain_ = EnergyFeedbackDomain::NONE;
  Workspace workspace_{};
  PowerControlData power_control_data_{};
};
