#pragma once

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

/**
 * @brief 固定维度递归最小二乘估计器
 *
 * 更新采用临时参数和临时协方差，只有全部结果有限且协方差有效时才提交，
 * 因此异常样本不会破坏上一组可信参数。
 */
template <std::size_t DIMENSION>
class RLS {
 public:
  static_assert(DIMENSION > 0U, "RLS dimension must be positive");

  using ParamVector = Eigen::Matrix<float, DIMENSION, 1>;
  using Matrix = Eigen::Matrix<float, DIMENSION, DIMENSION>;

  RLS() = delete;

  explicit RLS(float initial_covariance, float forgetting_factor,
               const ParamVector& initial_params = ParamVector::Zero())
      : initial_covariance_(PositiveOr(initial_covariance, 1000.0f)),
        forgetting_factor_(ValidForgettingFactor(forgetting_factor)) {
    Reset(initial_params);
  }

  explicit RLS(float initial_covariance, float forgetting_factor,
               const float (&initial_params)[DIMENSION])
      : RLS(initial_covariance, forgetting_factor,
            ParamVectorFromArray(initial_params)) {}

  /** @brief 重置参数和协方差。 */
  void Reset(const ParamVector& initial_params = ParamVector::Zero()) {
    params_ = ProjectParams(SanitizeVector(initial_params));
    covariance_ = Matrix::Identity() * initial_covariance_;
  }

  void Reset(const float (&initial_params)[DIMENSION]) {
    Reset(ParamVectorFromArray(initial_params));
  }

  /** @brief 设置参数物理边界，并立即投影当前参数。 */
  void SetParamBounds(const ParamVector& lower, const ParamVector& upper) {
    for (std::size_t index = 0; index < DIMENSION; ++index) {
      lower_bounds_[index] = std::isfinite(lower[index])
                                 ? lower[index]
                                 : -std::numeric_limits<float>::max();
      upper_bounds_[index] = std::isfinite(upper[index])
                                 ? upper[index]
                                 : std::numeric_limits<float>::max();
      if (upper_bounds_[index] < lower_bounds_[index]) {
        std::swap(upper_bounds_[index], lower_bounds_[index]);
      }
    }
    params_ = ProjectParams(params_);
  }

  void SetParamBounds(const float (&lower)[DIMENSION],
                      const float (&upper)[DIMENSION]) {
    SetParamBounds(ParamVectorFromArray(lower), ParamVectorFromArray(upper));
  }

  /**
   * @brief 使用一个样本更新参数
   * @return 更新被接受时返回 true；无效或退化样本返回 false。
   */
  bool Update(const ParamVector& sample, float actual_output) {
    if (!sample.allFinite() || !std::isfinite(actual_output)) {
      return false;
    }

    const float EXCITATION = sample.squaredNorm();
    if (!std::isfinite(EXCITATION) || EXCITATION <= 1.0e-12f) {
      return false;
    }

    const ParamVector FULL_COVARIANCE_TIMES_SAMPLE = covariance_ * sample;
    const float FULL_DENOMINATOR =
        forgetting_factor_ + sample.dot(FULL_COVARIANCE_TIMES_SAMPLE);
    if (!std::isfinite(FULL_DENOMINATOR) || FULL_DENOMINATOR <= 1.0e-12f) {
      return false;
    }

    const ParamVector UNCONSTRAINED_GAIN =
        FULL_COVARIANCE_TIMES_SAMPLE / FULL_DENOMINATOR;
    const float PREDICTION = sample.dot(params_);
    const float ERROR = actual_output - PREDICTION;
    const ParamVector UNCONSTRAINED_DELTA = UNCONSTRAINED_GAIN * ERROR;
    if (!UNCONSTRAINED_GAIN.allFinite() || !std::isfinite(ERROR) ||
        !UNCONSTRAINED_DELTA.allFinite()) {
      return false;
    }

    ParamVector free_mask = ParamVector::Ones();
    for (std::size_t index = 0; index < DIMENSION; ++index) {
      const bool BLOCKED_AT_UPPER = params_[index] >= upper_bounds_[index] &&
                                    UNCONSTRAINED_DELTA[index] > 0.0f;
      const bool BLOCKED_AT_LOWER = params_[index] <= lower_bounds_[index] &&
                                    UNCONSTRAINED_DELTA[index] < 0.0f;
      if (BLOCKED_AT_UPPER || BLOCKED_AT_LOWER) {
        free_mask[index] = 0.0f;
      }
    }

    const ParamVector BLOCKED_MASK = ParamVector::Ones() - free_mask;
    const Matrix FREE_SELECTOR = free_mask.asDiagonal();
    const Matrix BLOCKED_SELECTOR = BLOCKED_MASK.asDiagonal();
    const Matrix FREE_COVARIANCE =
        (FREE_SELECTOR * covariance_ * FREE_SELECTOR).eval();
    const ParamVector FREE_SAMPLE = sample.cwiseProduct(free_mask);
    const ParamVector FREE_COVARIANCE_TIMES_SAMPLE =
        FREE_COVARIANCE * FREE_SAMPLE;
    const float DENOMINATOR =
        forgetting_factor_ + FREE_SAMPLE.dot(FREE_COVARIANCE_TIMES_SAMPLE);
    if (!std::isfinite(DENOMINATOR) || DENOMINATOR <= 1.0e-12f) {
      return false;
    }

    const ParamVector GAIN = FREE_COVARIANCE_TIMES_SAMPLE / DENOMINATOR;
    const ParamVector DELTA = GAIN * ERROR;
    if (!GAIN.allFinite() || !DELTA.allFinite()) {
      return false;
    }

    ParamVector component_scale = ParamVector::Ones();
    for (std::size_t index = 0; index < DIMENSION; ++index) {
      const float SPAN = upper_bounds_[index] - lower_bounds_[index];
      if (std::isfinite(SPAN) && SPAN >= 0.0f) {
        const float MAX_STEP = SPAN * MAX_PARAMETER_STEP_FRACTION;
        if (std::fabs(DELTA[index]) > MAX_STEP) {
          component_scale[index] = std::min(component_scale[index],
                                            MAX_STEP / std::fabs(DELTA[index]));
        }
        if (DELTA[index] > 0.0f) {
          const float REMAINING = upper_bounds_[index] - params_[index];
          const float BOUNDARY_SCALE = std::max(0.0f, REMAINING / DELTA[index]);
          component_scale[index] =
              std::min(component_scale[index], BOUNDARY_SCALE);
        } else if (DELTA[index] < 0.0f) {
          const float REMAINING = params_[index] - lower_bounds_[index];
          const float BOUNDARY_SCALE =
              std::max(0.0f, REMAINING / -DELTA[index]);
          component_scale[index] =
              std::min(component_scale[index], BOUNDARY_SCALE);
        }
      }
      if (!std::isfinite(component_scale[index]) ||
          component_scale[index] < 0.0f) {
        return false;
      }
    }

    const ParamVector EFFECTIVE_GAIN = GAIN.cwiseProduct(component_scale);
    const ParamVector CANDIDATE_PARAMS =
        ProjectParams(params_ + EFFECTIVE_GAIN * ERROR);
    if (!(EFFECTIVE_GAIN.array().abs() > 0.0f).any() ||
        !EFFECTIVE_GAIN.allFinite() || !CANDIDATE_PARAMS.allFinite()) {
      return false;
    }

    const Matrix ATTENUATION =
        Matrix::Identity() - EFFECTIVE_GAIN * FREE_SAMPLE.transpose();
    const Matrix RETAINED_BLOCKED_COVARIANCE =
        (BLOCKED_SELECTOR * covariance_ * BLOCKED_SELECTOR).eval();
    Matrix candidate_covariance =
        (ATTENUATION * FREE_COVARIANCE * ATTENUATION.transpose()) /
            forgetting_factor_ +
        EFFECTIVE_GAIN * EFFECTIVE_GAIN.transpose() +
        RETAINED_BLOCKED_COVARIANCE;
    candidate_covariance =
        (0.5f * (candidate_covariance + candidate_covariance.transpose()))
            .eval();
    if (!candidate_covariance.allFinite()) {
      return false;
    }
    Eigen::LLT<Matrix> covariance_decomposition(candidate_covariance);
    if (covariance_decomposition.info() != Eigen::Success) {
      return false;
    }
    params_ = CANDIDATE_PARAMS;
    covariance_ = candidate_covariance;
    return true;
  }

  bool Update(const float (&sample)[DIMENSION], float actual_output) {
    return Update(ParamVectorFromArray(sample), actual_output);
  }

  void SetParamVector(const ParamVector& params) {
    params_ = ProjectParams(SanitizeVector(params));
  }

  void SetParamVector(const float (&params)[DIMENSION]) {
    SetParamVector(ParamVectorFromArray(params));
  }

  const ParamVector& GetParamVector() const { return params_; }

 private:
  static constexpr float MAX_PARAMETER_STEP_FRACTION = 0.05f;

  static float PositiveOr(float value, float fallback) {
    return std::isfinite(value) && value > 0.0f ? value : fallback;
  }

  static float ValidForgettingFactor(float value) {
    return std::isfinite(value) && value > 0.0f && value <= 1.0f ? value
                                                                 : 0.9999f;
  }

  static ParamVector ParamVectorFromArray(const float (&values)[DIMENSION]) {
    return Eigen::Map<const ParamVector>(values);
  }

  static ParamVector SanitizeVector(const ParamVector& vector) {
    return vector.unaryExpr(
        [](const float VALUE) { return std::isfinite(VALUE) ? VALUE : 0.0f; });
  }

  ParamVector ProjectParams(const ParamVector& params) const {
    return params.cwiseMax(lower_bounds_).cwiseMin(upper_bounds_);
  }

  float initial_covariance_;
  float forgetting_factor_;
  ParamVector lower_bounds_ =
      ParamVector::Constant(-std::numeric_limits<float>::max());
  ParamVector upper_bounds_ =
      ParamVector::Constant(std::numeric_limits<float>::max());
  ParamVector params_ = ParamVector::Zero();
  Matrix covariance_ = Matrix::Zero();
};
