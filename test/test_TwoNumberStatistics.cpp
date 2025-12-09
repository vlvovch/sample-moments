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

static void ExpectNear(double value, double expected, double tol = 1e-12) {
  ASSERT_NEAR(value, expected, tol);
}

TEST(TwoNumberStatistics, JointMomentsAndCumulantsSmallSample) {
  // Deterministic pairs: (1,0), (2,1), (3,1)
  std::vector<std::pair<int,int>> obs = { {1,0}, {2,1}, {3,1} };

  TwoNumberStatistics stats(/*nmax=*/4);
  stats.AddObservations(obs);

  // Means
  ExpectNear(stats.GetMean1(), 2.0);
  ExpectNear(stats.GetMean2(), 2.0 / 3.0);

  // Second joint moments
  ExpectNear(stats.GetJointMoment(2, 0), (1*1 + 2*2 + 3*3) / 3.0);
  ExpectNear(stats.GetJointMoment(0, 2), (0*0 + 1*1 + 1*1) / 3.0);
  ExpectNear(stats.GetJointMoment(1, 1), (1*0 + 2*1 + 3*1) / 3.0);

  // Central moments
  ExpectNear(stats.GetJointCentralMoment(2, 0), 2.0 / 3.0); // variance of first variable
  ExpectNear(stats.GetJointCentralMoment(0, 2), 2.0 / 9.0); // variance of second variable

  // Joint cumulants: kappa_{1,1} equals covariance
  ExpectNear(stats.GetJointCumulant(1, 1), stats.GetJointCentralMoment(1, 1));

  // A simple ratio to exercise ratio path (should be finite and match direct division)
  double k20 = stats.GetJointCumulant(2, 0);
  double k11 = stats.GetJointCumulant(1, 1);
  ExpectNear(stats.GetJointCumulantRatio(2, 0, 1, 1), k20 / k11);

  // Joint factorial moments and cumulants
  // F_{1,0} = <N1>, F_{0,1} = <N2>, F_{1,1} = <N1 N2>
  ExpectNear(stats.GetJointFactorialMoment(1, 0), stats.GetMean1());
  ExpectNear(stats.GetJointFactorialMoment(0, 1), stats.GetMean2());
  ExpectNear(stats.GetJointFactorialMoment(1, 1), stats.GetJointMoment(1, 1));

  // Joint factorial cumulant kappa_{1,1} matches covariance for these orders
  ExpectNear(stats.GetJointFactorialCumulant(1, 1), stats.GetJointCumulant(1, 1));
  // Ratio exercises factorial cumulant ratio path
  ExpectNear(stats.GetJointFactorialCumulantRatio(1, 1, 1, 0),
             stats.GetJointFactorialCumulant(1, 1) / stats.GetJointFactorialCumulant(1, 0));
}
