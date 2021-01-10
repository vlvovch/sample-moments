# sample-moments

**sample-moments** is a header-only C++ library 
to calculate sample averages and standard error estimates of 
(joint) (higher-order) moments, central moments and 
cumulants for one or two random variables.

Such a task arises when one tries to estimate various properties 
of a distribution from a finite set of observations sampled from
that distribution.
Common examples of such properties are given by the leading 
moments of the distribution, for instance the mean, variance,
skewness, kurtosis etc.
And while the calculation of the sample mean and its standard error 
is well known and straightforward, estimation of sample statistics 
(and their errors) that are functions of higher-order moments is more 
involved.

This library calculates sample average of any (joint) (central) moment
or cumulants for one or two variables from a set of observations.
The standard error estimate is calculated as well, 
to leading order in 1/√n<sub>obs</sub>
via the Delta method (see e.g. M. Kendall and A. Stuart, 
[''The advanced theory of statistics''](https://www.amazon.com/Kendalls-Advanced-Theory-Statistics-Distribution/dp/0470665300)).

## Prerequisites

A compiler with C++11 support. 

CMake v3.8+ to build samples and tests.

## Usage

Include [`SampleMoments.h`](include/SampleMoments.h) header file
to use the library in C++ code.
Then use either `NumberStatistics` or `TwoNumberStatistics` objects
to calculate sample moments and/or cumulants 
from single or two variable observations.
For example (see [`example_basic.cpp`](examples/example_basic.cpp):
```cpp
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
  int number_of_observations = 100;

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
```

A sample output is the following:
```
                      Sample mean: 9.99956 +- 0.00157491
                  Sample variance: 0.248034 +- 0.00110527
                  Sample skewness: -0.000286803 +- 0.00093758
           Sample excess kurtosis: -0.000880853 +- 0.000896272
Sample normalized excess kurtosis: -0.00355133 +- 0.00361298
```

See [examples](examples/) for more use cases and annotated header files for all the features.


## Notes

- By default, calculation of sample moments is supported up to 16th order,
and standard error estimates up to 8th order. 
  To increase this limit, e.g. increase 16/8 to 32/16 call the constructor as follows: 
  `SampleMoments::NumberStatistics(observations, 32)`.
- The library uses double-precision floating-point arithmetic to store the values of the moments.
Beware of the round-off errors, especially when the first moment is large but the (higher-order) central moments are small.
  See [normal_distribution.cpp](examples/normal_distribution.cpp) for an example on how to deal with this.

*Copyright (C) 2021 Volodymyr Vovchenko*
