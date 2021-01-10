/*
 * sample-moments library
 *
 * Copyright (c) 2021 Volodymyr Vovchenko
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>
 */
#include <SampleMoments.h>  // The library
#include <random>           // Random number generation
#include <iostream>         // Output

// Sample usage with the normal distribution
int main() {
  // Parameters of the normal distribution
  double mean  = 10.;
  double sigma = 0.5;

  // Random number generator from the normal distribution
  std::mt19937 rng;
  std::normal_distribution<double> distribution(mean, sigma);

  // The number of observations
  int number_of_observations = 100000;

  // Populate the observations
  std::vector<double> observations;
  for(int i = 0; i < number_of_observations; ++i) {
    double N = distribution(rng);
    observations.push_back(N);
  }

  // Analyze the observations
  SampleMoments::NumberStatistics stats(observations);

  std::cout << "                      Sample mean: " << stats.GetMean()
    << " +- " << stats.GetMeanError() << std::endl;
  std::cout << "                  Sample variance: " << stats.GetVariance()
    << " +- " << stats.GetVarianceError() << std::endl;
  std::cout << "                  Sample skewness: " << stats.GetCentralMoment(3)
    << " +- " << stats.GetCentralMomentError(3) << std::endl;
  std::cout << "           Sample excess kurtosis: " << stats.GetCumulant(4)
    << " +- " << stats.GetCumulantError(4) << std::endl;
  std::cout << "Sample normalized excess kurtosis: " << stats.GetCumulantRatio(4, 2)
            << " +- " << stats.GetCumulantRatioError(4, 2) << std::endl;

  return 0;
}