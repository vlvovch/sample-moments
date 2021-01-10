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

/// Cumulants of the binomial distribution, up to 6th order
double binomial_cumulant(int k, int n, double p)
{
  if (k == 1)
    return n * p;
  if (k == 2)
    return n * p * (1. - p);
  if (k == 3)
    return n * p * (1. - p) * (1. - 2.*p);
  if (k == 4)
    return n * (1 - p)*p*((1 - 2*p)*(1 - p) - (1 - 2*p)*p - 2*(1 - p)*p);
  if (k == 5)
    return n * (1 - p)*p*((1 - p)*p*(-2 - 4*(1 - p) + 8*p) 
    + (1 - p)*((1 - 2*p)*(1 - p) - (1 - 2*p)*p - 2*(1 - p)*p) 
    - p*((1 - 2*p)*(1 - p) - (1 - 2*p)*p - 2*(1 - p)*p));
  if (k == 6)
    return n * (1 - p)*p*((1 - p)*p*(-2*(1 - 2*p)*(1 - p) 
    + 2*(1 - 2*p)*p + 16*(1 - p)*p 
    + 2*(1 - p)*(-2 - 4*(1 - p) + 8*p) 
    - 2*p*(-2 - 4*(1 - p) + 8*p)) 
    + (1 - p)*((1 - p)*p*(-2 - 4*(1 - p) + 8*p) 
    + (1 - p)*((1 - 2*p)*(1 - p) - (1 - 2*p)*p - 2*(1 - p)*p) 
    - p*((1 - 2*p)*(1 - p) - (1 - 2*p)*p - 2*(1 - p)*p)) 
    - p*((1 - p)*p*(-2 - 4*(1 - p) + 8*p) + (1 - p)*((1 - 2*p)*(1 - p) 
    - (1 - 2*p)*p - 2*(1 - p)*p) - p*((1 - 2*p)*(1 - p) - (1 - 2*p)*p 
    - 2*(1 - p)*p)));
  return 0.;
}

/// Calculates cumulants and error estimates from observations sampled from the binomial distribution
/// Counts how many times the deviation of the sampled cumulants from the cumulants is within one-sigma level
/// For accurately estimated standard errors the fraction of sample cumulants within one sigma of the true values should be about 68.2%
int main(int argc, char* argv[]) {
  // Number of observation batches
  int iters = 1000;
  if (argc > 1)
    iters = atoi(argv[1]);

  // Size of the sample
  int sample_size = 10000;
  if (argc > 2)
    sample_size = atoi(argv[2]);

    
  // Parameters of the binomial distribution
  int n = 10;
  double p = 0.7;
  if (argc > 3)
    n = atoi(argv[3]);

  if (argc > 4)
    p = atof(argv[4]);

  std::cout << "Generating " << iters << " iterations of sample size " << sample_size << " numbers from the binomial distribution with "
    << "n = " << n << " and " << "p = " << p << std::endl;
  

  // Random number generator
  std::mt19937 rng;
  rng.seed(1);
  std::binomial_distribution<int> distribution(n, p);


  // Create and populate the statistics
  SampleMoments::NumberStatistics stats(12);

  //stats.SetMeanShift(n * p);

  int kmax = 6;
  std::vector<int> one_sigma_counts(kmax + 1, 0);

  for(int ii = 0; ii < iters; ++ii) {
    stats.Clear();
    for (int i = 0; i < sample_size; i++) {
      const int observation = distribution(rng);
      stats.AddObservation(observation);
    }

    for(int k = 1; k <= kmax; ++k) {
      double dev = (stats.GetCumulant(k) - binomial_cumulant(k, n, p)) / stats.GetCumulantError(k);
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