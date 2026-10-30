#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "PowerControl.hpp"

class TestCan final : public LibXR::CAN {
 public:
  LibXR::ErrorCode SetConfig(const Configuration&) override {
    return LibXR::ErrorCode::OK;
  }

  uint32_t GetClockFreq() const override { return 0U; }

  LibXR::ErrorCode AddMessage(const ClassicPack&) override {
    return LibXR::ErrorCode::OK;
  }

  void InjectPower(float power) {
    ClassicPack frame{};
    frame.id = SuperPowerProtocol::FEEDBACK_ID;
    frame.type = Type::STANDARD;
    frame.dlc = 8U;
    const SuperPowerProtocol::FeedbackData data{0U, power, 120U, 128U};
    std::memcpy(frame.data, &data, sizeof(data));
    OnMessage(frame, false);
  }
};

int main() {
  LibXR::PlatformInit();
  const auto referee_topic =
      LibXR::Topic::CreateTopic<Referee::ChassisPack>("chassis_ref");
  auto* can = new TestCan();
  auto* superpower = new SuperPower(*can);
  auto* power = new PowerControl(superpower);

  const auto initial = power->GetPowerControlData();
  assert(std::fabs(initial.rls_diagnostics.covariance_trace - 2.0e-5f) <
         1.0e-10f);

  const float commands[4] = {};
  const float rpm[4] = {1000.0f, 1000.0f, 1000.0f, 1000.0f};
  for (int i = 0; i < 300; ++i) {
    assert(power->SetMotorFeedback3508(commands, rpm, 4));
    if (i % 5 == 0) {
      can->InjectPower(12.0f);
    }
    LibXR::Thread::Sleep(1);
  }
  power->OutputLimit();
  const auto result = power->GetPowerControlData();
  assert(result.supercap_healthy);
  assert(result.rls_diagnostics.covariance_trace > 0.0f);
  assert(result.rls_diagnostics.covariance_trace <
         initial.rls_diagnostics.covariance_trace);
  assert(result.m3508_speed_loss > initial.m3508_speed_loss);
  assert(result.rls_diagnostics.numerical_rejection_count == 0U);

  // PowerControl owns a lifetime-long background thread.
  std::_Exit(EXIT_SUCCESS);
}
