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

  return 0;
}