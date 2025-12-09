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
#include <vector>
#include <cmath>
#include <limits>

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

double binomial_factorial_cumulant(int k, int n, double p) {
  if (k < 1)
    return 0.0;
  double fact = 1.0;
  for (int i = 1; i < k; ++i)
    fact *= static_cast<double>(i);
  double sign = ((k + 1) % 2 == 0) ? 1.0 : -1.0; // (-1)^{k+1}
  return n * fact * sign * std::pow(p, k);
}

/// Calculates cumulants and error estimates from observations sampled from the normal distribution
/// and compares the results with the true cumulant values
int main(int argc, char* argv[]) {

  // Size of the sample
  int sample_size = 1000000;
  if (argc > 1)
    sample_size = atoi(argv[1]);


  // Parameters of the binomial distribution
  int n = 10;
  double p = 0.7;
  if (argc > 2)
    n = atoi(argv[2]);

  if (argc > 3)
    p = atof(argv[3]);

  std::cout << "Generating " << sample_size << " numbers from the binomial distribution with "
    << "n = " << n << " and " << "p = " << p << std::endl;


  // Random number generator
  std::mt19937 rng;
  rng.seed(1);
  std::binomial_distribution<int> distribution(n, p);


  // Create and populate the statistics
  SampleMoments::NumberStatistics stats(12);

  //stats.SetMeanShift(n * p);

  std::map<int,int> hist;

  for (int i = 0; i < sample_size; i++) {
    const int observation = distribution(rng);
    hist[observation]++;
    stats.AddObservation(observation);
  }

  for(auto el : hist) {
    std::cout << el.first << " " << el.second / static_cast<double>(sample_size)
    << " +- " << sqrt(el.second) / static_cast<double>(sample_size)
    << std::endl;
  }

  // Check the first four cumulants
  for(int k = 1; k <= 6; ++k) {
    std::cout << std::setw(14) << "\\kappa_" << k << ": ";
    std::cout << std::setw(14) << "Expected:" << " ";
    std::cout << std::setw(14) << binomial_cumulant(k, n, p) << std::endl;
    std::cout << std::setw(14) << " " << "   ";
    std::cout << std::setw(14) << "Observed:" << " ";
    std::cout << std::setw(14) << stats.GetCumulant(k) << " +- ";
    std::cout << std::setw(14) << stats.GetCumulantError(k) << " ";
    std::cout << std::setw(14) << "Deviation " << "(sigmas):" << " ";
    std::cout << (binomial_cumulant(k, n, p) - stats.GetCumulant(k)) /stats.GetCumulantError(k) << " ";
    std::cout << std::endl;
    std::cout << std::endl;
  }

  // Cumulant ratios of interest
  struct RatioSpec { int num; int den; };
  std::vector<RatioSpec> ratio_specs = { {2, 1}, {3, 2}, {4, 2}, {6, 2} };

  std::cout << "Cumulant ratios:" << std::endl;
  for (const auto &ratio : ratio_specs) {
    const double expected_num = binomial_cumulant(ratio.num, n, p);
    const double expected_den = binomial_cumulant(ratio.den, n, p);
    double expected_ratio = (expected_den == 0.0) ? std::numeric_limits<double>::quiet_NaN() : expected_num / expected_den;
    const double observed = stats.GetCumulantRatio(ratio.num, ratio.den);
    const double observed_err = stats.GetCumulantRatioError(ratio.num, ratio.den);
    double deviation = std::numeric_limits<double>::quiet_NaN();
    if (observed_err != 0.0 && !std::isnan(observed_err))
      deviation = (expected_ratio - observed) / observed_err;

    std::cout << std::setw(10) << "\\kappa_" << ratio.num << "/\\kappa_" << ratio.den << ": ";
    std::cout << std::setw(14) << "Expected:" << " ";
    std::cout << std::setw(14) << expected_ratio << std::endl;
    std::cout << std::setw(14) << " " << "   ";
    std::cout << std::setw(14) << "Observed:" << " ";
    std::cout << std::setw(14) << observed << " +- ";
    std::cout << std::setw(14) << observed_err << " ";
    std::cout << std::setw(14) << "Deviation " << "(sigmas):" << " ";
    std::cout << deviation << " ";
    std::cout << std::endl;
    std::cout << std::endl;
  }

  // Factorial cumulants (more natural for count data) and a few ratios
  std::cout << "Factorial cumulants:" << std::endl;
  for (int k = 1; k <= 4; ++k) {
    const double expected_fc = binomial_factorial_cumulant(k, n, p);
    std::cout << std::setw(14) << "C~" << k << ": ";
    std::cout << std::setw(14) << "Expected:" << " ";
    std::cout << std::setw(14) << expected_fc << std::endl;
    std::cout << std::setw(14) << " " << "   ";
    std::cout << std::setw(14) << "Observed:" << " ";
    std::cout << std::setw(14) << stats.GetFactorialCumulant(k) << " +- ";
    std::cout << std::setw(14) << stats.GetFactorialCumulantError(k) << " ";
    std::cout << std::setw(14) << "Deviation " << "(sigmas):" << " ";
    std::cout << (expected_fc - stats.GetFactorialCumulant(k)) / stats.GetFactorialCumulantError(k) << " ";
    std::cout << std::endl;
    std::cout << std::endl;
  }

  std::vector<RatioSpec> fc_ratio_specs = { {2, 1}, {3, 1}, {4, 1} };
  std::cout << "Factorial cumulant ratios:" << std::endl;
  for (const auto &ratio : fc_ratio_specs) {
    const double expected_num = binomial_factorial_cumulant(ratio.num, n, p);
    const double expected_den = binomial_factorial_cumulant(ratio.den, n, p);
    double expected_ratio = (expected_den == 0.0) ? std::numeric_limits<double>::quiet_NaN() : expected_num / expected_den;
    const double observed = stats.GetFactorialCumulantRatio(ratio.num, ratio.den);
    const double observed_err = stats.GetFactorialCumulantRatioError(ratio.num, ratio.den);
    double deviation = std::numeric_limits<double>::quiet_NaN();
    if (observed_err != 0.0 && !std::isnan(observed_err))
      deviation = (expected_ratio - observed) / observed_err;

    std::cout << std::setw(10) << "C~" << ratio.num << "/C~" << ratio.den << ": ";
    std::cout << std::setw(14) << "Expected:" << " ";
    std::cout << std::setw(14) << expected_ratio << std::endl;
    std::cout << std::setw(14) << " " << "   ";
    std::cout << std::setw(14) << "Observed:" << " ";
    std::cout << std::setw(14) << observed << " +- ";
    std::cout << std::setw(14) << observed_err << " ";
    std::cout << std::setw(14) << "Deviation " << "(sigmas):" << " ";
    std::cout << deviation << " ";
    std::cout << std::endl;
    std::cout << std::endl;
  }

  return 0;
}
