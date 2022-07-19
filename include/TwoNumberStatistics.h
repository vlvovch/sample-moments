/*
 * sample-moments library
 *
 * Copyright (c) 2021 Volodymyr Vovchenko
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>
 */
#ifndef SAMPLEMOMENTS_TWONUMBERSTATISTICS_H
#define SAMPLEMOMENTS_TWONUMBERSTATISTICS_H

#include "MomentsTransformations.h"
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>

namespace SampleMoments {
  /**
   * \brief A class implementing calculation of various sample statistics and standard errors for two variables.
   *        For example this includes joint (higher-order) moments, central moments, cumulants, and ratios of such quantities.
   *
   *
   */
  class TwoNumberStatistics {
    /// The maximum order of (joint) moments stored
    int m_Nmax;

    /// The cumulative sums of <N_1^k N_2^m> values collected from all observations
    std::vector<std::vector<double> > m_MomentSums;

    /// The sample means of the moments, <N_1^k N_2^m>. Populated by CalculateMoments()
    std::vector<std::vector<double> > m_Moments;

    /// The sample means of the central moments, <(N_1-<N_1>)^k (N_2-<N_2>)^m>. Populated by CalculateMoments()
    std::vector<std::vector<double> > m_CentralMoments;

    /// Shifts all moments by a consant values
    /// Can be useful to avoid large round-off errors if the expected mean is known
    /// Does not affect the central moments
    double m_MeanShift1, m_MeanShift2;

    /// The total number of observations
    int64_t m_NumberOfObservations;

    /// Same as m_NumberOfObservations
    double m_nE;

    /// Whether the moments have been computed from the sums
    bool m_MomentsComputed;

    /// The pre-computed binomial coefficients
    std::vector<std::vector<int64_t> > m_BinomialCoefficients;
  public:
    /**
     * \brief Construct a new NumberStatistics object.
     * 
     * \param nmax   The maximum order of moments to be stored
     * 
     * For the error estimation nmax must be at least twice larger than the order of the moment to be estimated
     */
    TwoNumberStatistics(int nmax = 16) {
      Reset(nmax);
    }

    /**
     * \brief Construct a new NumberStatistics object.
     *
     * \param observations   A vector of observations
     * \param nmax           The maximum order of moments to be stored
     *
     * For the error estimation nmax must be at least twice larger than the order of the moment to be estimated
     */
    template<class T1, class T2>
    TwoNumberStatistics(const std::vector<std::pair<T1, T2>> &observations, int nmax = 16) {
      Reset(nmax);
      AddObservations(observations);
    }

    /// Destructor. Does nothing
    virtual ~TwoNumberStatistics() {}

    /**
     * \brief Resets the observations to zero.
     * 
     * \param nmax   The maximum order of moments to be stored
     * 
     * For the error estimation nmax must be at least twice larger than the order of the moment to be estimated
     */
    virtual void Reset(int nmax) {
      assert(nmax >= 1);
      m_Nmax = nmax;
      m_MomentSums = std::vector<std::vector<double> >(nmax + 1, std::vector<double>(nmax + 1, 0.));
      m_Moments = std::vector<std::vector<double> >(nmax + 1, std::vector<double>(nmax + 1, 0.));
      m_CentralMoments = std::vector<std::vector<double> >(nmax + 1, std::vector<double>(nmax + 1, 0.));
      m_NumberOfObservations = 0;
      m_nE = 0.;
      m_MeanShift1 = 0.;
      m_MeanShift2 = 0.;

      // Pre-compute the binomial coefficients
      CalculateBinomialCoefficients();

      m_MomentsComputed = false;
    }

    bool IsMeanShifted() {
      return (m_MeanShift1 != 0.0 || m_MeanShift2 != 0.0);
    }

    /// Adds an observation with a pair integer-valued numbers N1, N2. 
    /// These values are cast to double-precision floating-point numbers
    template<typename T1, typename T2>
    void AddObservation(const T1& N1, const T2& N2) {
      double dN1 = static_cast<double>(N1);
      double dN2 = static_cast<double>(N2);
      AddObservation(dN1, dN2);
    }

    /// Adds an observation  with a pair real-valued numbers N1, N2. 
    void AddObservation(double N1, double N2) {
      UpdateMomentSumsWithObservation(N1, N2, m_MomentSums);

      m_NumberOfObservations++;

      m_MomentsComputed = false;
    }

    /// Adds several integer-valued observations
    /// \param observations A vector of observations
    template<typename T1, typename T2>
    void AddObservations(const std::vector<std::pair<T1, T2>> &observations) {
      std::vector<std::pair<double, double>> observations_to_double;
      observations_to_double.reserve(observations.size());
      for (const auto &N : observations)
        observations_to_double.push_back({static_cast<double>(N.first), static_cast<double>(N.second)});
      AddObservations(observations_to_double);
    }

    /// Adds several real-valued observations
    /// \param observations A vector of observations
    void AddObservations(const std::vector<std::pair<double, double>> &observations) {
      // Compute moments_sum from the whole observations vector to minimize round-off errors
      auto moment_sums = std::vector<std::vector<double> >(m_Nmax, std::vector<double>(m_Nmax, 0.));

      for (const auto &N : observations) {
        UpdateMomentSumsWithObservation(N.first, N.second, moment_sums);
      }

      // Now add to the whole sums
      for (int i = 0; i < m_Nmax + 1; ++i) {
        for (int j = 0; j < m_Nmax + 1; ++j) {
          m_MomentSums[i][j] += moment_sums[i][j];
        }
      }

      m_NumberOfObservations += static_cast<int64_t>(observations.size());

      m_MomentsComputed = false;
    }


    /// Returns the joint moment <N1^i N2^j>
    double GetJointMoment(int i, int j) {
      assert(i <= m_Nmax && j <= m_Nmax);
      CalculateMoments();
      return m_Moments[i][j];
    }

    double GetJointMoment(const std::pair<int,int> &index) { return GetJointMoment(index.first, index.second); }

    /// Returns the sample covariance cov(<N1^r N2^s>,<N1^u N2^v>)
    /// In accordance with M. Kendall and A. Stuart, ``The advanced theory of statistics''
    double GetJointMomentsSampleCovariance(int r, int s, int u, int v) {
      assert(r + u <= m_Nmax && s + v <= m_Nmax);
      CalculateMoments();
      return (m_Moments[r + u][s + v] - m_Moments[r][s] * m_Moments[u][v]) / m_nE;
    }

    double GetJointMomentsSampleCovariance(const std::pair<int,int> &index1, const std::pair<int,int> &index2)
    {
      return GetJointMomentsSampleCovariance(index1.first, index1.second, index2.first, index2.second);
    }


    /// Returns the error estimate for the joint moment <N1^i N2^j>
    double GetJointMomentError(int i, int j) {
      return std::sqrt(GetJointMomentsSampleCovariance(i, j, i, j));
    }

    double GetJointMomentError(const std::pair<int,int> &index) { return GetJointMomentError(index.first, index.second); }


    /// Returns the ratio of joint moments <N1^i1 N2^j1> / <N1^i2 N2^j2>
    double GetJointMomentRatio(int i1, int j1, int i2, int j2) {
      return GetJointCentralMoment(i1, j1) / GetJointCentralMoment(i2, j2);
    }

    double GetJointMomentRatio(const std::pair<int,int> &index1, const std::pair<int,int> &index2)
    {
      return GetJointMomentRatio(index1.first, index1.second, index2.first, index2.second);
    }

    /// Returns the ratio of joint moments <N1^i1 N2^j1> / <N1^i2 N2^j2>
    double GetJointMomentRatioError(int i1, int j1, int i2, int j2) {
      assert(i1 <= m_Nmax && j1 <= m_Nmax);
      assert(i2 <= m_Nmax && j2 <= m_Nmax);

      double c1 = GetJointMoment(i1, j1);
      double c2 = GetJointMoment(i2, j2);
      double c1Dev = GetJointMomentsSampleCovariance(i1, j1, i1, j1);
      double c2Dev = GetJointMomentsSampleCovariance(i2, j2, i2, j2);
      double c1c2cov = GetJointMomentsSampleCovariance(i1, j1, i2, j2);

      return std::abs(c1 / c2) * std::sqrt(c1Dev / c1 / c1 + c2Dev / c2 / c2 - 2. * c1c2cov / c1 / c2);
    }

    double GetJointMomentRatioError(const std::pair<int,int> &index1, const std::pair<int,int> &index2)
    {
      return GetJointMomentRatioError(index1.first, index1.second, index2.first, index2.second);
    }

    /// Returns the r-th central moment, \mu_i,j =  <(N1-<N1>)^i (N2-<N2>)^j>
    double GetJointCentralMoment(int i, int j) {
      assert(i <= m_Nmax && j <= m_Nmax);
      CalculateMoments();
      return m_CentralMoments[i][j];
    }

    double GetJointCentralMoment(const std::pair<int,int> &index) { return GetJointCentralMoment(index.first, index.second); }

    /// Returns the sample covariance for joint central moments cov(\mu_r,s, \mu_u,v)
    /// In accordance with M. Kendall and A. Stuart, ``The advanced theory of statistics''
    double GetJointCentralMomentsSampleCovariance(int r, int s, int u, int v) {
      CalculateMoments();

      assert(m_Nmax >= 2);
      assert(r + u <= m_Nmax);
      assert(s + v <= m_Nmax);

      double ret = 0.;

      ret += m_CentralMoments[r + u][s + v] - m_CentralMoments[r][s] * m_CentralMoments[u][v];
      if (r > 0 && u > 0)
        ret += r * u * m_CentralMoments[2][0] * m_CentralMoments[r - 1][s] * m_CentralMoments[u - 1][v];
      if (s > 0 && v > 0)
        ret += s * v * m_CentralMoments[0][2] * m_CentralMoments[r][s - 1] * m_CentralMoments[u][v - 1];
      if (r > 0 && v > 0)
        ret += r * v * m_CentralMoments[1][1] * m_CentralMoments[r - 1][s] * m_CentralMoments[u][v - 1];
      if (s > 0 && u > 0)
        ret += s * u * m_CentralMoments[1][1] * m_CentralMoments[r][s - 1] * m_CentralMoments[u - 1][v];
      if (u > 0)
        ret -= u * m_CentralMoments[r + 1][s] * m_CentralMoments[u - 1][v];
      if (v > 0)
        ret -= v * m_CentralMoments[r][s + 1] * m_CentralMoments[u][v - 1];
      if (r > 0)
        ret -= r * m_CentralMoments[r - 1][s] * m_CentralMoments[u + 1][v];
      if (s > 0)
        ret -= s * m_CentralMoments[r][s - 1] * m_CentralMoments[u][v + 1];

      return ret / m_nE;
    }

    double GetJointCentralMomentsSampleCovariance(const std::pair<int,int> &index1, const std::pair<int,int> &index2)
    {
      return GetJointCentralMomentsSampleCovariance(index1.first, index1.second, index2.first, index2.second);
    }

    /// Returns the error estimate for the r-th central moment \mu_i,j =  <(N1-<N1>)^i (N2-<N2>)^j>
    double GetJointCentralMomentError(int i, int j) {
      return std::sqrt(GetJointCentralMomentsSampleCovariance(i, j, i, j));
    }

    double GetJointCentralMomentError(const std::pair<int,int> &index) { return GetJointCentralMomentError(index.first, index.second); }


    /// Returns the ratio of joint central moments \mu_i1,j1 / \mu_i2,j2
    double GetJointCentralMomentRatio(int i1, int j1, int i2, int j2) {
      return GetJointCentralMoment(i1, j1) / GetJointCentralMoment(i2, j2);
    }

    double GetJointCentralMomentRatio(const std::pair<int,int> &index1, const std::pair<int,int> &index2)
    {
      return GetJointCentralMomentRatio(index1.first, index1.second, index2.first, index2.second);
    }

    /// Returns the error estimate for the ratio of joint central moments \mu_i1,j1 / \mu_i2,j2
    double GetJointCentralMomentRatioError(int i1, int j1, int i2, int j2) {
      double c1 = GetJointCentralMoment(i1, j1);
      double c2 = GetJointCentralMoment(i2, j2);
      double c1Dev = GetJointCentralMomentsSampleCovariance(i1, j1, i1, j1);
      double c2Dev = GetJointCentralMomentsSampleCovariance(i2, j2, i2, j2);
      double c1c2cov = GetJointCentralMomentsSampleCovariance(i1, j1, i2, j2);

      return std::abs(c1 / c2) * std::sqrt(c1Dev / c1 / c1 + c2Dev / c2 / c2 - 2. * c1c2cov / c1 / c2);
    }

    double GetJointCentralMomentRatioError(const std::pair<int,int> &index1, const std::pair<int,int> &index2)
    {
      return GetJointCentralMomentRatioError(index1.first, index1.second, index2.first, index2.second);
    }

    /// Returns the joint cumulant \kappa_i,j
    double GetJointCumulant(int i, int j) {
      return CalculateJointCumulantFromMoments(i, j);
    }

    double GetJointCumulant(const std::pair<int,int> &index) { return GetJointCumulant(index.first, index.second); }

    /// Returns the error estimate for the joint cumulant \kappa_i,j
    double GetJointCumulantError(int i, int j) {
      return CalculateJointCumulantErrorFromMoments(i, j);
    }

    double GetJointCumulantError(const std::pair<int,int> &index) { return GetJointCumulantError(index.first, index.second); }


    /// Returns the joint cumulant ratio \kappa_i1,j2 / \kappa_i2,j2
    double GetJointCumulantRatio(int i1, int j1, int i2, int j2) {
      return GetJointCumulant(i1, j1) / GetJointCumulant(i2, j2);
    }

    double GetJointCumulantRatio(const std::pair<int,int> &index1, const std::pair<int,int> &index2)
    {
      return GetJointCumulantRatio(index1.first, index1.second, index2.first, index2.second);
    }

    /// Returns the error estimate for the joint cumulant ratio \kappa_i1,j2 / \kappa_i2,j2
    double GetJointCumulantRatioError(int i1, int j1, int i2, int j2) {
      return CalculateJointCumulantRatioErrorFromMoments(i1, j1, i2, j2);
    }

    double GetJointCumulantRatioError(const std::pair<int,int> &index1, const std::pair<int,int> &index2)
    {
      return GetJointCumulantRatioError(index1.first, index1.second, index2.first, index2.second);
    }

    /// Returns the sample mean for the first variable
    double GetMean1() {
      return GetJointMoment(1, 0);
    }

    /// Returns the error estimate for the sample mean for the first variable
    double GetMean1Error() {
      return GetJointMomentError(1, 0);
    }

    /// Returns the sample mean for the second variable
    double GetMean2() {
      CalculateMoments();
      return GetJointMoment(0, 1);
    }

    /// Returns the error estimate for the sample mean for the second variable
    double GetMean2Error() {
      CalculateMoments();
      return GetJointMomentError(0, 1);
    }

    /// Set the shifts of means to use when accumulating statistics
    /// Can be useful to avoid round-off error when dealing with higher-order central moments and cumulants
    void SetMeanShifts(double meanshift1, double meanshift2) {
      if (m_NumberOfObservations == 0) {
        m_MeanShift1 = meanshift1;
        m_MeanShift2 = meanshift2;
        return;
      }

      auto shifted_moments = m_Moments;
      for (int i = 0; i < m_Nmax + 1; ++i)
        for (int j = 0; j < m_Nmax + 1; ++j)
          shifted_moments[i][j] = m_MomentSums[i][j] / m_NumberOfObservations;

      double delta1 = m_MeanShift1 - meanshift1;
      double delta2 = m_MeanShift2 - meanshift2;
      m_MeanShift1 = meanshift1;
      m_MeanShift2 = meanshift2;

      auto new_shifted_moments = shifted_moments;
      for (int n = 0; n <= m_Nmax; ++n) {
        for (int m = 0; m <= m_Nmax; ++m) {

          new_shifted_moments[n][m] = 0.;
          double Nshift1 = 1.;
          for (int k = 0; k <= n; ++k) {

            double Nshift2 = 1.;
            for (int l = 0; l <= m; ++l) {

              new_shifted_moments[n][m] += BinomialCoefficient(n, k) * BinomialCoefficient(m, l)
                                           * shifted_moments[n - k][m - l] * Nshift1 * Nshift2;

              Nshift2 *= delta2;

            }
            Nshift1 *= delta1;
          }
        }
      }

      for (int i = 0; i < m_Nmax + 1; ++i)
        for (int j = 0; j < m_Nmax + 1; ++j)
          m_MomentSums[i][j] = new_shifted_moments[i][j] * m_NumberOfObservations;

      m_MomentsComputed = false;
      CalculateMoments();
    }

  protected:
    /// Pre-compute the needed binomial coefficients
    void CalculateBinomialCoefficients() {
      m_BinomialCoefficients = std::vector<std::vector<int64_t> >(m_Nmax + 1);

      m_BinomialCoefficients[0] = std::vector<int64_t>(1, 1);
      m_BinomialCoefficients[1] = std::vector<int64_t>(2, 1);

      for (int n = 2; n < m_Nmax + 1; ++n) {
        m_BinomialCoefficients[n] = std::vector<int64_t>(n + 1, 1);
        m_BinomialCoefficients[n][0] = m_BinomialCoefficients[n][n] = 1;
        for (int mm = 1; mm < n; ++mm)
          m_BinomialCoefficients[n][mm] = m_BinomialCoefficients[n - 1][mm] + m_BinomialCoefficients[n - 1][mm - 1];
      }
    }

    /// Return the binomial coefficient C_n,k. Use pre-computed values if possible, or use recursion otherwise
    int64_t BinomialCoefficient(int n, int k) const {
      if (k < 0 || k > n) return 0;
      if (k == 0 || k == n) return 1;
      if (n < m_Nmax + 1) return m_BinomialCoefficients[n][k];
      return BinomialCoefficient(n - 1, k) + BinomialCoefficient(n - 1, k - 1);
    }

    /// Returns the (central) moment of order i,j
    /// \param  i               Order of the first variable
    /// \param  j               Order of the second variable
    /// \param  central_moment  Whether computed moment is central (true) or ordinary (false)
    double CalculateJointMoment(int i, int j, bool central_moment) {
      if (central_moment)
        return GetJointCentralMoment(i, j);
      else
        return GetJointMoment(i, j);
    }

    /// Returns the sample covariance for {r,s} and {u,v} (central) moments
    /// \param  r               Order of the first moment
    /// \param  s               Order of the second moment
    /// \param  u               Order of the first moment
    /// \param  v               Order of the second moment
    /// \param  central_moment  Whether the covariance is for central (true) or ordinary (false) moments
    double CalculateJointMomentsSampleCovariance(int r, int s, int u, int v, bool central_moment) {
      if (central_moment)
        return GetJointCentralMomentsSampleCovariance(r, s, u, v);
      else
        return GetJointMomentsSampleCovariance(r, s, u, v);
    }

    /// Returns the i,j-th order cumulant \kappa_i,j computed from (central) moments
    /// \param  i                     Order of the first variable
    /// \param  j                     Order of the second variable
    /// \param  use_central_moments   Whether to use central (true) or ordinary (false) moments for the calculation
    double CalculateJointCumulantFromMoments(int i, int j, bool use_central_moments = true) {
      /// If cumulant order is less than two, must use ordinary moments
      if (i + j < 2 && use_central_moments)
        return CalculateJointCumulantFromMoments(i, j, false);

      std::vector<int> indices;
      for (int ii = 0; ii < i; ++ii)
        indices.push_back(0);
      for (int ii = 0; ii < j; ++ii)
        indices.push_back(1);
      auto cumulant_vs_centralmoments = JointCumulantToCentralMoments(indices, 2, !use_central_moments);
      double ret = 0.;
      for (const auto &term : cumulant_vs_centralmoments) {
        double tmp = static_cast<double>(term.second);
        for (const auto &multiplier : term.first) {
          tmp *= CalculateJointMoment(multiplier[0], multiplier[1], use_central_moments);
        }
        ret += tmp;
      }
      return ret;
    }

    /// Returns the error estimate i,j-th order cumulant \kappa_i,j computed from (central) moments
    /// \param  i                     Order of the first variable
    /// \param  j                     Order of the second variable
    /// \param  use_central_moments   Whether to use central (true) or ordinary (false) moments for the calculation
    double CalculateJointCumulantErrorFromMoments(int i, int j, bool use_central_moments = true) {
      if (i + j < 2 && use_central_moments)
        return CalculateJointCumulantErrorFromMoments(i, j, false);

      std::vector<int> indices;
      for (int ii = 0; ii < i; ++ii)
        indices.push_back(0);
      for (int ii = 0; ii < j; ++ii)
        indices.push_back(1);
      auto cumulant_vs_centralmoments = JointCumulantToCentralMoments(indices, 2, !use_central_moments);
      double ret = 0.;
      for (const auto &term1 : cumulant_vs_centralmoments) {
        double tmp1 = static_cast<double>(term1.second);
        for (const auto &multiplier1 : term1.first) {
          tmp1 *= CalculateJointMoment(multiplier1[0], multiplier1[1], use_central_moments);
        }

        for (const auto &multiplier1 : term1.first) {
          double terr1 = tmp1 / CalculateJointMoment(multiplier1[0], multiplier1[1], use_central_moments);

          for (const auto &term2 : cumulant_vs_centralmoments) {
            double tmp2 = static_cast<double>(term2.second);
            for (const auto &multiplier2 : term2.first) {
              tmp2 *= CalculateJointMoment(multiplier2[0], multiplier2[1], use_central_moments);
            }

            for (const auto &multiplier2 : term2.first) {
              double terr2 = tmp2 / CalculateJointMoment(multiplier2[0], multiplier2[1], use_central_moments);

              ret += terr1 * terr2 *
                     CalculateJointMomentsSampleCovariance(multiplier1[0], multiplier1[1], multiplier2[0],
                                                           multiplier2[1], use_central_moments);
            }
          }

        }
      }
      return std::sqrt(ret);
    }

    /// Returns the joint cumulant ratio \kappa_i1,j1/\kappa_i2,j2 computed from (central) moments
    /// \param  i1                    Order of the first variable of the first cumulant
    /// \param  j1                    Order of the second variable of the first cumulant
    /// \param  i2                    Order of the first variable of the second cumulant
    /// \param  j2                    Order of the second variable of the second cumulant
    /// \param  use_central_moments   Whether to use central (true) or ordinary (false) moments for the calculation
    double CalculateJointCumulantRatioFromMoments(int i1, int j1, int i2, int j2, bool use_central_moments = true) {
      return CalculateJointCumulantFromMoments(i1, j1, use_central_moments) /
             CalculateJointCumulantFromMoments(i2, j2, use_central_moments);
    }

    /// Returns the error estimate for the joint cumulant ratio \kappa_i1,j1/\kappa_i2,j2 computed from (central) moments
    /// \param  i1                    Order of the first variable of the first cumulant
    /// \param  j1                    Order of the second variable of the first cumulant
    /// \param  i2                    Order of the first variable of the second cumulant
    /// \param  j2                    Order of the second variable of the second cumulant
    /// \param  use_central_moments   Whether to use central (true) or ordinary (false) moments for the calculation
    double
    CalculateJointCumulantRatioErrorFromMoments(int i1, int j1, int i2, int j2, bool use_central_moments = true) {
      /// If either of the cumulants' order is less than two, must use ordinary moments
      if ((i1 + j1 < 2 || i2 + j2 < 2) && use_central_moments)
        return CalculateJointCumulantRatioErrorFromMoments(i1, j1, i2, j2, false);

      double numerator = CalculateJointCumulantFromMoments(i1, j1, use_central_moments);
      double denominator = CalculateJointCumulantFromMoments(i2, j2, use_central_moments);

      std::vector<int> indices;
      for (int ii = 0; ii < i1; ++ii)
        indices.push_back(0);
      for (int ii = 0; ii < j1; ++ii)
        indices.push_back(1);
      auto cumulant_vs_centralmoments_numerator = JointCumulantToCentralMoments(indices, 2, !use_central_moments);

      indices.clear();
      for (int ii = 0; ii < i2; ++ii)
        indices.push_back(0);
      for (int ii = 0; ii < j2; ++ii)
        indices.push_back(1);
      auto cumulant_vs_centralmoments_denominator = JointCumulantToCentralMoments(indices, 2, !use_central_moments);

      double ret = 0.0;

      // First loop over the terms in the numerator
      for (const auto &term1num : cumulant_vs_centralmoments_numerator) {
        double tmp1 = static_cast<double>(term1num.second);
        for (const auto &multiplier1num : term1num.first) {
          tmp1 *= CalculateJointMoment(multiplier1num[0], multiplier1num[1], use_central_moments);
        }

        for (const auto &multiplier1num : term1num.first) {
          double terr1 = tmp1 / CalculateJointMoment(multiplier1num[0], multiplier1num[1], use_central_moments);

          // Second term from the numerator
          for (const auto &term2num : cumulant_vs_centralmoments_numerator) {
            double tmp2 = static_cast<double>(term2num.second);
            for (const auto &multiplier2num : term2num.first) {
              tmp2 *= CalculateJointMoment(multiplier2num[0], multiplier2num[1], use_central_moments);
            }

            for (const auto &multiplier2num : term2num.first) {
              double terr2 = tmp2 / CalculateJointMoment(multiplier2num[0], multiplier2num[1], use_central_moments);

              ret += (terr1 / denominator) * (terr2 / denominator)
                     * CalculateJointMomentsSampleCovariance(multiplier1num[0], multiplier1num[1], multiplier2num[0],
                                                             multiplier2num[1], use_central_moments);
            }
          }

          // Second term from the denominator, double the contribution due to symmetry
          for (const auto &term2den : cumulant_vs_centralmoments_denominator) {
            double tmp2 = static_cast<double>(term2den.second);
            for (const auto &multiplier2den : term2den.first) {
              tmp2 *= CalculateJointMoment(multiplier2den[0], multiplier2den[1], use_central_moments);
            }

            for (const auto &multiplier2den : term2den.first) {
              double terr2 = tmp2 / CalculateJointMoment(multiplier2den[0], multiplier2den[1], use_central_moments);

              ret += 2. * (terr1 / denominator) * (terr2 * numerator * (-1.) / denominator / denominator)
                     * CalculateJointMomentsSampleCovariance(multiplier1num[0], multiplier1num[1], multiplier2den[0],
                                                             multiplier2den[1], use_central_moments);
            }
          }

        }
      }

      // Denominator-denominator terms
      for (const auto &term1den : cumulant_vs_centralmoments_denominator) {
        double tmp1 = static_cast<double>(term1den.second);
        for (const auto &multiplier1den : term1den.first) {
          tmp1 *= CalculateJointMoment(multiplier1den[0], multiplier1den[1], use_central_moments);
        }

        for (const auto &multiplier1den : term1den.first) {
          double terr1 = tmp1 / CalculateJointMoment(multiplier1den[0], multiplier1den[1], use_central_moments);

          // Second term from the denominator
          for (const auto &term2den : cumulant_vs_centralmoments_denominator) {
            double tmp2 = static_cast<double>(term2den.second);
            for (const auto &multiplier2den : term2den.first) {
              tmp2 *= CalculateJointMoment(multiplier2den[0], multiplier2den[1], use_central_moments);
            }

            for (const auto &multiplier2den : term2den.first) {
              double terr2 = tmp2 / CalculateJointMoment(multiplier2den[0], multiplier2den[1], use_central_moments);

              ret += (terr1 * numerator * (-1.) / denominator / denominator) *
                     (terr2 * numerator * (-1.) / denominator / denominator)
                     * CalculateJointMomentsSampleCovariance(multiplier1den[0], multiplier1den[1], multiplier2den[0],
                                                             multiplier2den[1], use_central_moments);
            }
          }
        }
      }

      return std::sqrt(ret);
    }

    /// Calculate the values of the ordinary and central moments from the currently accumulated observations
    void CalculateMoments() {
      if (m_MomentsComputed)
        return;

      if (m_NumberOfObservations > 0) {
        for (int i = 0; i < m_Nmax + 1; ++i)
          for (int j = 0; j < m_Nmax + 1; ++j)
            m_Moments[i][j] = m_MomentSums[i][j] / m_NumberOfObservations;
      } else {
        m_MomentsComputed = true;
        m_nE = 0.0;
        return;
      }

      m_Moments[0][0] = 1.;

      for (int n = 0; n <= m_Nmax; ++n) {
        for (int m = 0; m <= m_Nmax; ++m) {
          m_CentralMoments[n][m] = 0.;
          double Nav1 = 1.;
          for (int k = 0; k <= n; ++k) {
            double Nav2 = 1.;
            for (int l = 0; l <= m; ++l) {

              m_CentralMoments[n][m] += BinomialCoefficient(n, k) * BinomialCoefficient(m, l)
                                 * m_Moments[n - k][m - l] * Nav1 * Nav2;

              Nav2 *= -m_Moments[0][1];
            }
            Nav1 *= -m_Moments[1][0];;
          }
        }
      }

      // if m_MomentSums are mean-shifted, need to recalculate moments from the shifted moments
      if (IsMeanShifted()) {
        auto shifted_moments = m_Moments;
        for (int n = 0; n <= m_Nmax; ++n) {
          for (int m = 0; m <= m_Nmax; ++m) {

            m_Moments[n][m] = 0.;
            double Nshift1 = 1.;
            for (int k = 0; k <= n; ++k) {

              double Nshift2 = 1.;
              for (int l = 0; l <= m; ++l) {

                m_Moments[n][m] += BinomialCoefficient(n, k) * BinomialCoefficient(m, l)
                                   * shifted_moments[n - k][m - l] * Nshift1 * Nshift2;

                Nshift2 *= m_MeanShift2;

              }
              Nshift1 *= m_MeanShift1;

            }
          }
        }
      }

      m_nE = static_cast<double>(m_NumberOfObservations);

      m_MomentsComputed = true;
    }

  private:
    void UpdateMomentSumsWithObservation(double N1, double N2, std::vector<std::vector<double> > &moment_sums) {
      if (IsMeanShifted()) {
        N1 -= m_MeanShift1;
        N2 -= m_MeanShift2;
      }

      double tmn1 = 1.;
      for (int i = 0; i <= m_Nmax; ++i) {
        double tmn2 = 1.;
        for (int j = 0; j <= m_Nmax; ++j) {
          moment_sums[i][j] += tmn1 * tmn2;
          tmn2 *= N2;
        }
        tmn1 *= N1;
      }
    }

    void UpdateMomentSumsWithObservation(int N1, int N2, std::vector<std::vector<double> > &moment_sums) {
      UpdateMomentSumsWithObservation(static_cast<double>(N1), static_cast<double>(N2), moment_sums);
    }
  };

} // namespace SampleMoments

#endif
