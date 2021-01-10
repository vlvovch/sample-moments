/*
 * sample-moments library
 *
 * Copyright (c) 2021 Volodymyr Vovchenko
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>
 */
#include "../include/SampleMoments.h"
#include <random>
#include <iostream>
#include <iomanip>

double normal_distribution_cumulant(int k, double mean, double sigma)
{
  if (k == 1)
    return mean;
  if (k == 2)
    return sigma * sigma;
  return 0.;
}

/// Calculates cumulants and error estimates from observations sampled from the normal distribution
/// Counts how many times the deviation of the sampled cumulants from the cumulants is within one-sigma level
/// For accurately estimated standard errors the fraction of sample cumulants within one sigma of the true values should be about 68.2%
int main(int argc, char* argv[]) {
  // Iterations
  int iters = 10000;
  if (argc > 1)
    iters = atoi(argv[1]);

  // Size of the sample
  int sample_size = 10000;
  if (argc > 2)
    sample_size = atoi(argv[2]);

    
  // Parameters of the normal distribution
  double mean  = 10.;
  double sigma = 0.5;
  if (argc > 2)
    mean = atof(argv[2]);

  if (argc > 3)
    sigma = atof(argv[3]);

  std::cout << "Generating " << iters << " iterations of sample size " << sample_size << " numbers from the normal distribution with "
    << "\\mu = " << mean << " and " << "\\sigma = " << sigma << std::endl;
  std::cout << std::endl;

  // Random number generator
  std::mt19937 rng;
  rng.seed(1);
  std::normal_distribution<double> distribution(mean, sigma);

  // Gather statistics with a shifted mean to avoid round-off errors
  SampleMoments::NumberStatistics stats(12);

  int kmax = 6;
  std::vector<int> one_sigma_counts(kmax + 1, 0);

  for(int ii = 0; ii < iters; ++ii) {
    stats.Clear();
    stats.SetMeanShift(mean);

    for (int i = 0; i < sample_size; i++) {
      const double observation = distribution(rng);

      stats.AddObservation(observation);
    }

    for(int k = 1; k <= kmax; ++k) {
      double dev = (stats.GetCumulant(k) - normal_distribution_cumulant(k, mean, sigma)) / stats.GetCumulantError(k);
      if (abs(dev) < 1.)
        one_sigma_counts[k]++;
    }

  }

  std::cout << "Fraction of samples within one-sigma: " << std::endl;
  for(int k = 1; k <= 6; ++k) {
    std::cout << std::setw(14) << "\\kappa_" << k << ": ";
    std::cout << std::setw(14) << static_cast<double>(one_sigma_counts[k]) / iters << std::endl;
    std::cout << std::endl;
  }


  return 0;
}