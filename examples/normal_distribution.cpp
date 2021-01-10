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
/// and compares the results with the true cumulant values
/// Due to the fact that higher-order (>3) cumulants of the normal distribution are zero
/// this creates numerical issues with floating-point operations where one gets large cancellations
/// This issues is addressed using a technique of shifted means as illustrated in this example
int main(int argc, char* argv[]) {

  // Size of the sample
  int sample_size = 1000000;
  if (argc > 1)
    sample_size = atoi(argv[1]);

    
  // Parameters of the normal distribution
  double mean  = 10.;
  double sigma = 0.5;
  if (argc > 2)
    mean = atof(argv[2]);

  if (argc > 3)
    sigma = atof(argv[3]);

  std::cout << "Generating " << sample_size << " numbers from the normal distribution with "
    << "\\mu = " << mean << " and " << "\\sigma = " << sigma << std::endl;
  std::cout << std::endl;

  // Random number generator
  std::mt19937 rng;
  rng.seed(1);
  std::normal_distribution<double> distribution(mean, sigma);


  // Create and populate the statistics
  SampleMoments::NumberStatistics stats(12);

  // Gather statistics with a shifted mean to avoid round-off errors
  SampleMoments::NumberStatistics stats_with_mean_shift(12);


  stats_with_mean_shift.SetMeanShift(mean);

  for (int i = 0; i < sample_size; i++) {
    const double observation = distribution(rng);

    stats.AddObservation(observation);

    stats_with_mean_shift.AddObservation(observation);
  }

  // Check the first six cumulants
  // First without the mean shift
  std::cout << "Sampled cumulants and their errors without mean shift: " << std::endl;
  for(int k = 1; k <= 6; ++k) {
    std::cout << std::setw(14) << "\\kappa_" << k << ": ";
    std::cout << std::setw(14) << "Expected:" << " ";
    std::cout << std::setw(14) << normal_distribution_cumulant(k, mean, sigma) << std::endl;
    std::cout << std::setw(14) << " " << "   ";
    std::cout << std::setw(14) << "Observed:" << " ";
    std::cout << std::setw(14) << stats.GetCumulant(k) << " +- ";
    std::cout << std::setw(14) << stats.GetCumulantError(k) << " ";
    std::cout << std::setw(14) << "Deviation " << "(sigmas):" << " ";
    std::cout << (normal_distribution_cumulant(k, mean, sigma) - stats.GetCumulant(k)) /stats.GetCumulantError(k) << " ";
    std::cout << std::endl;
    std::cout << std::endl;
  }

  // Now the statistics gathered with the mean shift
  std::cout << "Sampled cumulants and their errors with mean shift: " << std::endl;
  for(int k = 1; k <= 6; ++k) {
    std::cout << std::setw(14) << "\\kappa_" << k << ": ";
    std::cout << std::setw(14) << "Expected:" << " ";
    std::cout << std::setw(14) << normal_distribution_cumulant(k, mean, sigma) << std::endl;
    std::cout << std::setw(14) << " " << "   ";
    std::cout << std::setw(14) << "Observed:" << " ";
    std::cout << std::setw(14) << stats_with_mean_shift.GetCumulant(k) << " +- ";
    std::cout << std::setw(14) << stats_with_mean_shift.GetCumulantError(k) << " ";
    std::cout << std::setw(14) << "Deviation " << "(sigmas):" << " ";
    std::cout << (normal_distribution_cumulant(k, mean, sigma) - stats_with_mean_shift.GetCumulant(k)) /stats_with_mean_shift.GetCumulantError(k) << " ";
    std::cout << std::endl;
    std::cout << std::endl;
  }

  return 0;
}