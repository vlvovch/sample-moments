/*
 * sample-moments library
 *
 * Copyright (c) 2021 Volodymyr Vovchenko
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>
 */
#include "gtest/gtest.h"
#include "SampleMoments.h"

using namespace SampleMoments;

// Helper to compare doubles with a tight tolerance appropriate for exact rational expectations
static void ExpectNear(double value, double expected, double tol = 1e-12) {
  ASSERT_NEAR(value, expected, tol);
}

TEST(NumberStatistics, MomentsAndCumulantsSmallSample) {
  // Deterministic small sample: observations = {1, 2, 3}
  std::vector<int> obs = {1, 2, 3};
  NumberStatistics stats(obs, /*nmax=*/6);

  // Basic moments
  ExpectNear(stats.GetMean(), 2.0);
  ExpectNear(stats.GetVariance(), 2.0 / 3.0);

  // Central moments
  ExpectNear(stats.GetCentralMoment(3), 0.0);
  ExpectNear(stats.GetCentralMoment(4), 2.0 / 3.0); // for this sample, mu4 = 2/3

  // Cumulants (for univariate, first cumulant is mean, second is variance, third is mu3, fourth is k4)
  ExpectNear(stats.GetCumulant(1), 2.0);
  ExpectNear(stats.GetCumulant(2), 2.0 / 3.0);
  ExpectNear(stats.GetCumulant(3), 0.0);
  ExpectNear(stats.GetCumulant(4), -2.0 / 3.0); // for this sample, kappa4 = mu4 - 3*mu2^2 = 2/3 - 3*(4/9) = -2/3

  // Ratios
  ExpectNear(stats.GetCumulantRatio(2, 1), (2.0 / 3.0) / 2.0);
  ExpectNear(stats.GetCumulantRatio(4, 2), (-2.0 / 3.0) / (2.0 / 3.0));

  // Factorial moments / cumulants for non-negative data
  ExpectNear(stats.GetFactorialMoment(1), stats.GetMean());
  ExpectNear(stats.GetFactorialMoment(2), 8.0 / 3.0); // E[N(N-1)] over {1,2,3} is (0+2+6)/3 = 8/3
  ExpectNear(stats.GetFactorialMomentRatio(2, 1), (8.0 / 3.0) / 2.0);

  ExpectNear(stats.GetFactorialCumulant(1), stats.GetMean());
  ExpectNear(stats.GetFactorialCumulant(2), -4.0 / 3.0); // for one variable: FC2 = C2 - C1
  ExpectNear(stats.GetFactorialCumulantRatio(2, 1), (-4.0 / 3.0) / 2.0);
}

TEST(NumberStatistics, MeanShiftDoesNotChangeResults) {
  // Same observations but accumulate with a mean shift equal to the true mean
  std::vector<double> obs = {1.0, 2.0, 3.0};

  NumberStatistics stats_no_shift(obs, /*nmax=*/6);

  NumberStatistics stats_shifted(/*nmax=*/6);
  stats_shifted.SetMeanShift(2.0);
  stats_shifted.AddObservations(obs);

  // Verify key quantities agree
  ExpectNear(stats_no_shift.GetMean(), stats_shifted.GetMean());
  ExpectNear(stats_no_shift.GetVariance(), stats_shifted.GetVariance());
  ExpectNear(stats_no_shift.GetCumulant(4), stats_shifted.GetCumulant(4));
}
