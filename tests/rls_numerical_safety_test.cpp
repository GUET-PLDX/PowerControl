#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>

#include "RLS.hpp"

static void RankOneExcitationKeepsCovarianceBounded() {
  constexpr float MAX_TRACE = 1.0e-3f;
  constexpr float MIN_EIGENVALUE = 1.0e-9f;
  RLS<2> rls(1.0e-5f, 0.99999f, {0.22f, 1.2f}, {MAX_TRACE, MIN_EIGENVALUE});
  rls.SetParamBounds({0.0f, 0.0f}, {5.0f, 20.0f});

  for (uint32_t i = 0; i < 2000000U; ++i) {
    assert(rls.Update({20.0f, 0.0f}, 4.4f));
    const auto& covariance = rls.GetCovariance();
    assert(covariance.allFinite());
    assert(covariance.trace() <= MAX_TRACE);
    Eigen::Matrix2f shifted = covariance;
    shifted.diagonal().array() -= MIN_EIGENVALUE;
    assert(Eigen::LLT<Eigen::Matrix2f>(shifted).info() == Eigen::Success);
  }

  const auto diagnostics = rls.GetDiagnostics();
  assert(diagnostics.trace_limit_count > 0U);
  assert(diagnostics.numerical_rejection_count == 0U);
  assert(std::fabs(diagnostics.covariance_trace - rls.GetCovariance().trace()) <
         1.0e-10f);
  assert(rls.GetCovariance()(0, 0) > 1.0e-8f);

  assert(rls.Update({100000.0f, 0.0f}, 22000.0f));
  assert(rls.GetDiagnostics().regularization_count > 0U);
  assert(rls.GetCovariance().trace() <= MAX_TRACE);

  for (uint32_t i = 0; i < 100U; ++i) {
    assert(rls.Update({0.0f, 10.0f}, 15.0f));
  }
  assert(rls.GetParamVector()[1] > 1.4f);
  assert(rls.GetParamVector()[1] <= 1.5f);
}

static void LimitedParameterStepsUseTheirActualGain() {
  RLS<2> rls(1.0f, 1.0f, {0.0f, 0.0f}, {100.0f, 1.0e-9f});
  rls.SetParamBounds({0.0f, 0.0f}, {1.0f, 100.0f});
  assert(rls.Update({1.0f, 1.0f}, 100.0f));

  assert(std::fabs(rls.GetParamVector()[0] - 0.05f) < 1.0e-6f);
  assert(std::fabs(rls.GetParamVector()[1] - 5.0f) < 1.0e-5f);
  const auto& covariance = rls.GetCovariance();
  assert(std::fabs(covariance(0, 0) - 0.99900075f) < 1.0e-6f);
  assert(std::fabs(covariance(0, 1) - -0.050425f) < 1.0e-6f);
  assert(std::fabs(covariance(1, 1) - 0.9075f) < 1.0e-6f);
}

static void RegularizationKeepsSmallEigenvaluesPositive() {
  RLS<2> rls(1.0e-5f, 1.0f, {0.22f, 1.2f}, {1.0e-3f, 1.0e-9f});
  assert(rls.Update({100000.0f, 0.0f}, 22000.0f));
  const auto diagnostics = rls.GetDiagnostics();
  assert(diagnostics.regularization_count > 0U);
  assert(diagnostics.numerical_rejection_count == 0U);
  Eigen::Matrix2f shifted = rls.GetCovariance();
  shifted.diagonal().array() -= 1.0e-9f;
  assert(Eigen::LLT<Eigen::Matrix2f>(shifted).info() == Eigen::Success);
}

static void InvalidSamplesPreserveTheLastTrustedState() {
  RLS<2> rls(1.0e-5f, 0.99999f, {0.22f, 1.2f}, {1.0e-3f, 1.0e-9f});
  assert(rls.Update({20.0f, 0.0f}, 4.4f));
  const auto params = rls.GetParamVector();
  const auto covariance = rls.GetCovariance();
  const float NAN_VALUE = std::numeric_limits<float>::quiet_NaN();
  const float MAX_VALUE = std::numeric_limits<float>::max();

  assert(!rls.Update({NAN_VALUE, 1.0f}, 1.0f));
  assert(!rls.Update({1.0f, 1.0f}, NAN_VALUE));
  assert(!rls.Update({MAX_VALUE, 1.0f}, 1.0f));
  assert(!rls.Update({0.0f, 0.0f}, 1.0f));
  assert(rls.GetParamVector() == params);
  assert(rls.GetCovariance() == covariance);
  assert(rls.GetDiagnostics().numerical_rejection_count == 3U);

  rls.Reset(params);
  assert(rls.GetCovariance().trace() <= 1.0e-3f);
  assert(rls.GetDiagnostics().numerical_rejection_count == 3U);
  assert(rls.Update({20.0f, 0.0f}, 4.4f));
}

static void BoundParametersRetainTheirCovariance() {
  RLS<2> rls(1.0f, 0.9f, {0.0f, 1.0f}, {100.0f, 1.0e-9f});
  rls.SetParamBounds({0.0f, 0.0f}, {1.0f, 100.0f});
  assert(rls.Update({1.0f, 1.0f}, 0.0f));
  assert(rls.GetParamVector()[0] == 0.0f);
  assert(rls.GetParamVector()[1] < 1.0f);
  assert(rls.GetCovariance()(0, 0) == 1.0f);
  assert(rls.GetCovariance()(0, 1) == 0.0f);
  assert(rls.GetCovariance()(1, 1) > 0.0f);
}

static void InvalidSafetySettingsRejectUpdates() {
  RLS<2> rls(1.0e-5f, 1.0f, {0.22f, 1.2f}, {1.0e-6f, 1.0e-9f});
  assert(!rls.IsCovarianceSafetyValid());
  const auto params = rls.GetParamVector();
  const auto covariance = rls.GetCovariance();
  assert(!rls.Update({20.0f, 0.0f}, 4.4f));
  assert(rls.GetParamVector() == params);
  assert(rls.GetCovariance() == covariance);
  assert(rls.GetDiagnostics().numerical_rejection_count == 1U);
}

static void OtherFixedDimensionsCanApplyTheSameLimit() {
  RLS<1> one(1.0e-5f, 0.99f, {0.0f}, {1.0e-3f, 1.0e-9f});
  RLS<3> three(1.0e-5f, 0.99f, {0.22f, 1.2f, 0.5f}, {1.0e-3f, 1.0e-9f});
  for (uint32_t i = 0; i < 3000U; ++i) {
    assert(one.Update({0.001f}, 0.0f));
    assert(three.Update({20.0f, 0.0f, 0.0f}, 4.4f));
    assert(one.GetCovariance().trace() <= 1.0e-3f);
    assert(three.GetCovariance().trace() <= 1.0e-3f);
  }
  assert(one.GetDiagnostics().trace_limit_count > 0U);
  assert(three.GetDiagnostics().trace_limit_count > 0U);
}

static void RotatedRankOneExcitationKeepsCovarianceBounded() {
  constexpr float MAX_TRACE = 1.0e-3f;
  constexpr float MIN_EIGENVALUE = 1.0e-9f;
  RLS<2> rls(1.0e-5f, 0.99999f, {0.22f, 1.2f}, {MAX_TRACE, MIN_EIGENVALUE});
  rls.SetParamBounds({0.0f, 0.0f}, {5.0f, 20.0f});
  for (uint32_t i = 0; i < 1000000U; ++i) {
    assert(rls.Update({14.0f, 14.0f}, 19.88f));
    const auto& covariance = rls.GetCovariance();
    assert(covariance.allFinite());
    assert(covariance.trace() <= MAX_TRACE);
    Eigen::Matrix2f shifted = covariance;
    shifted.diagonal().array() -= MIN_EIGENVALUE;
    assert(Eigen::LLT<Eigen::Matrix2f>(shifted).info() == Eigen::Success);
  }
  assert(rls.GetDiagnostics().trace_limit_count > 0U);
}

int main() {
  RankOneExcitationKeepsCovarianceBounded();
  LimitedParameterStepsUseTheirActualGain();
  RegularizationKeepsSmallEigenvaluesPositive();
  InvalidSamplesPreserveTheLastTrustedState();
  BoundParametersRetainTheirCovariance();
  InvalidSafetySettingsRejectUpdates();
  OtherFixedDimensionsCanApplyTheSameLimit();
  RotatedRankOneExcitationKeepsCovarianceBounded();
}
