/*
 * sample-moments library
 *
 * Copyright (c) 2021 Volodymyr Vovchenko
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>
 */
#ifndef SAMPLEMOMENTS_MOMENTSTRANSFORMATIONS_H
#define SAMPLEMOMENTS_MOMENTSTRANSFORMATIONS_H

#include <vector>
#include <map>
#include <algorithm>
#include <cassert>

namespace SampleMoments {

  // Faa di Bruno formula (single variable)
  typedef std::vector<int> Block;
  typedef std::vector<Block> Partition;

  /// Returns all partitions of a set {1,2,...,n} via a recursive procedure
  /// Each element contains a vector of blocks, each block is a vector of indices
  /// See https://en.wikipedia.org/wiki/Partition_of_a_set
  /// \param n  The dimension of the set
  static std::vector<Partition> PartitionsOfSet(int n) {
    if (n == 1)
      return {{{0}}};

    auto subpartitions = PartitionsOfSet(n - 1);
    auto ret = std::vector<Partition>();

    for (auto &el : subpartitions) {
      auto temp_element = el;
      temp_element.push_back({n - 1});
      ret.push_back(temp_element);
      temp_element = el;
      for (int i = 0; i < el.size(); ++i) {
        temp_element[i].push_back(n - 1);
        ret.push_back(temp_element);
        temp_element[i] = el[i];
      }
    }

    return ret;
  }

  /// Express an arbitrary joint cumulant \kappa(X_i0,X_i1,X_i2,...) in terms of joint (central) moments.
  /// Uses the multivariate version of the Faa di Bruno's formula https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula
  /// \param indices            A vector of 0-based indices of the joint cumulant
  /// \param dimensions         The total number of distinct random variable. 
  ///                           If set to -1, its minimum possible value is inferred from the contents of the provided cumulant \p indices
  /// \param ordinary_moments   Whether to express the cumulants in terms ordinary moments (true) or central moments (false). The default is central moments
  ///
  /// \return                   Returns the joint cumulant as a sum of products of various (central) moments.
  ///                           The return value is a map. 
  ///                           Each key-value pair correspond to a single term.
  ///                           The key corresponds to a vector of (central) moments. Each element of this vector is a vector of the moments' indices.
  ///                           The value is a numeral factor in front of each term. 
  static std::map<std::vector<std::vector<int>>, int64_t> JointCumulantToCentralMoments(
          const std::vector<int> &indices,
          int dimensions = -1,
          bool ordinary_moments = false
  ) {
    assert(indices.size() > 0);

    if (dimensions == -1) {
      for (int ind : indices) {
        assert(ind >= 0);
        dimensions = std::max(dimensions, ind);
      }
      dimensions++;
    }

    // Use the multivariate version of the Faa di Bruno's formula https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula
    auto ret = std::map<std::vector<std::vector<int>>, int64_t>();

    std::vector<int64_t> factorials(indices.size(), 1);
    for (int k = 1; k < factorials.size(); ++k)
      factorials[k] = k * factorials[k - 1];

    auto partitions_FdB = PartitionsOfSet(indices.size());
    for (const auto &partition : partitions_FdB) {
      int64_t multiplier = 1;
      int num_blocks = partition.size();
      if (num_blocks % 2 == 0)
        multiplier = -1;
      multiplier *= factorials[num_blocks - 1];

      std::vector<std::vector<int>> mults;
      bool zero_term = false;
      for (const auto &block : partition) {
        if (block.size() == 1) {
          zero_term = true;
        }

        auto inds = std::vector<int>(dimensions, 0);
        for (int ind : block)
          inds[indices[ind]]++;

        mults.push_back(inds);
      }
      if (zero_term && !ordinary_moments)
        continue;
      sort(mults.begin(), mults.end());
      ret[mults] += multiplier;
    }

    return ret;
  }

  /// Express an arbitrary joint cumulant in terms of the joint moments
  /// Calls JointCumulantToCentralMoments()
  static std::map<std::vector<std::vector<int>>, int64_t> JointCumulantToMoments(
          const std::vector<int> &indices,
          int dimensions = -1
  ) {
    return JointCumulantToCentralMoments(indices, dimensions, true);
  }

  /// Express an arbitrary factorial cumulant C~n in terms of ordinary cumulants \\kappa_k.
  /// Uses the combinartorial version of the Faa di Bruno's formula https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula
  /// \param n            Order of the factorial cumulant
  /// \return                   Returns the factorial cumulants as linear combination of ordinary cumulants.
  ///                           The return value is a map..
  ///                           Each key-value pair correspond to a single term.
  ///                           The key corresponds to the ordinary cumualant order.
  ///                           The value is a numeral factor in front of each term.
  static std::map<int,int64_t> FactorialCumulantsToOrdinaryCumulants(
          int n
  ) {


    // Use the multivariate version of the Faa di Bruno's formula https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula
    auto ret = std::map<int, int64_t>();

    std::vector<int64_t> factorials(n, 1);
    for (int k = 1; k < n; ++k)
      factorials[k] = k * factorials[k - 1];

    auto partitions_FdB = PartitionsOfSet(n);
    for (const auto &partition : partitions_FdB) {
      int num_blocks = partition.size();
      int64_t multiplier = 1;

      std::vector<std::vector<int>> mults;
      for (const auto &block : partition) {
        if (block.size()%2 == 0)
          multiplier *= -1;
        multiplier *= factorials[block.size() - 1];
      }
      ret[num_blocks] += multiplier;
    }

    return ret;
  }

  // Faa di Bruno formula (multivariate), see Appendix in https://arxiv.org/pdf/2106.13775
  typedef std::pair<int, Block> ColoredBlock;
  typedef std::vector<ColoredBlock> ColoredPartition;

  /// Precompute partitions of a set into 2D blocks
  std::vector< std::vector<ColoredPartition> > partitionsOfSet2DPrecomputed;

  /// Returns all partitions of a set {1,2,...,n} into m distinct colors of blocks via a recursive procedure
  /// Each element contains a vector of blocks, each block is a pair where the first element
  /// is the block color (0-indexed) and the second element is the vector of indices
  /// This is the extended version of the standard partition of a set to accomodate multivariate chain rule
  /// e.g. d^n f(g_1({x}),...,g_m({x}) / (d x_i1 ... d_x_im)
  /// \param n  The dimension of the set
  /// \param m  The number of block colors
  std::vector<ColoredPartition> ColoredPartitionsOfSet(int n, int m) {
    if (m == 2 && n < partitionsOfSet2DPrecomputed.size())
      return partitionsOfSet2DPrecomputed[n];

    if (n == 0)
      return std::vector<ColoredPartition>(); // empty set

    // single blocks of m types, each contains the single element
    if (n == 1) {
      auto ret = std::vector<ColoredPartition>();
      for (int j = 0; j < m; ++j)
        //ret.push_back({{j,{0}}});
        ret.push_back({ ColoredBlock ({j,{0}}) });
      return ret;
    }

    auto subpartitions = ColoredPartitionsOfSet(n - 1, m);
    auto ret = std::vector<ColoredPartition>();

    for (auto& el : subpartitions) {
      for (int j = 0; j < m; ++j) {
        auto temp_element = el;
        temp_element.push_back({ j,{n - 1} });
        sort(temp_element.begin(), temp_element.end());
        ret.push_back(temp_element);
      }
      auto temp_element = el;
      for (int i = 0; i < el.size(); ++i) {
        temp_element[i].second.push_back(n - 1);
        ret.push_back(temp_element);
        sort(temp_element.begin(), temp_element.end());
        temp_element[i] = el[i];
      }
    }


    return ret;
  }

  /// Precompute partitions of a 2D set
  void PrecomputePartitions(int nmax)
  {
    for (int n = 0; n <= nmax; ++n) {
      if (n >= partitionsOfSet2DPrecomputed.size())
        partitionsOfSet2DPrecomputed.push_back(ColoredPartitionsOfSet(n, 2));
    }
  }

  /// Express an arbitrary joint factorial cumulant C~i1,...,ik in terms of ordinary joint cumulants
  /// Uses the combinatorial version of the Faa di Bruno's formula https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula
  /// \param indices            A vector of 0-based indices of the joint cumulant
  /// \param dimensions         The total number of distinct random variable.
  ///                           If set to -1, its minimum possible value is inferred from the contents of the provided cumulant \p indices
  /// \return                   Returns the factorial cumulants as linear combination of ordinary cumulants.
  ///                           The return value is a map.
  ///                           Each key-value pair correspond to a single term.
  ///                           The key corresponds to the ordinary cumulant order.
  ///                           The value is a numerical factor in front of each term.
  static std::map<std::vector<int>,int64_t> JointFactorialCumulantsToOrdinaryCumulants(
          const std::vector<int> &indices,
          int dimensions = -1
  ) {

    if (dimensions == -1) {
      for (int ind : indices) {
        assert(ind >= 0);
        dimensions = std::max(dimensions, ind);
      }
      dimensions++;
    }

    // Use the multivariate version of the Faa di Bruno's formula https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula
    auto ret = std::map<std::vector<int>, int64_t>();

    std::vector<int64_t> factorials(indices.size(), 1);
    for (int k = 1; k < factorials.size(); ++k)
      factorials[k] = k * factorials[k - 1];

    auto partitions_FdB = ColoredPartitionsOfSet(indices.size(), dimensions);
    for (const auto &partition : partitions_FdB) {
      int num_blocks = partition.size();

      // Indices of the cumulant
      std::vector<int> Cinds(dimensions, 0);
      for (const auto &block : partition) {
        Cinds[block.first]++;
      }

      int64_t multiplier = 1;

      for (const auto &block : partition) {
        int color = block.first;
        bool fl = true;
        for(auto &el : block.second) {
          if (indices[el] != color)
            fl = false;
        }

        if (!fl) {
          multiplier = 0;
          break;
        }

        if (block.second.size()%2 == 0)
          multiplier *= -1;
        multiplier *= factorials[block.second.size() - 1];
      }
      if (multiplier != 0)
        ret[Cinds] += multiplier;
    }

    return ret;
  }

} // namespace SampleMoments

#endif