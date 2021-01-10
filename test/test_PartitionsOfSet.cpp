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

namespace {

  // The number of partitions of a set must be equal to the Bell number https://en.wikipedia.org/wiki/Partition_of_a_set
  TEST(PartitionsOfSet, BellNumbers) {
    EXPECT_EQ(PartitionsOfSet(1).size(), 1);
    EXPECT_EQ(PartitionsOfSet(2).size(), 2);
    EXPECT_EQ(PartitionsOfSet(3).size(), 5);
    EXPECT_EQ(PartitionsOfSet(4).size(), 15);
    EXPECT_EQ(PartitionsOfSet(5).size(), 52);
    EXPECT_EQ(PartitionsOfSet(6).size(), 203);
    EXPECT_EQ(PartitionsOfSet(7).size(), 877);
    EXPECT_EQ(PartitionsOfSet(8).size(), 4140);
    EXPECT_EQ(PartitionsOfSet(9).size(), 21147);
    EXPECT_EQ(PartitionsOfSet(10).size(), 115975);
    // too slow for larger sets
    //EXPECT_EQ(PartitionsOfSet(11).size(), 678570);
    //EXPECT_EQ(PartitionsOfSet(12).size(), 4213597);
    //EXPECT_EQ(PartitionsOfSet(13).size(), 27644437);
    //EXPECT_EQ(PartitionsOfSet(14).size(), 190899322);
    //EXPECT_EQ(PartitionsOfSet(15).size(), 1382958545);
    //EXPECT_EQ(PartitionsOfSet(16).size(), 10480142147ll);
  }

  TEST(PartitionsOfSet, Partitions) {
    auto partition = PartitionsOfSet(1);
    auto reference_partition = partition;
    reference_partition = {{{0}}};
    std::sort(partition.begin(), partition.end());
    for (auto& blocks : reference_partition) {
      std::sort(blocks.begin(), blocks.end());
    }
    std::sort(reference_partition.begin(), reference_partition.end());
    EXPECT_EQ(partition, reference_partition);

    partition = PartitionsOfSet(2);
    reference_partition = {{{0,1}}, {{0},{1}}};
    std::sort(partition.begin(), partition.end());
    for (auto& blocks : reference_partition) {
      std::sort(blocks.begin(), blocks.end());
    }
    std::sort(reference_partition.begin(), reference_partition.end());
    EXPECT_EQ(partition, reference_partition);

    partition = PartitionsOfSet(3);
    reference_partition = {{{0,1,2}}, {{0,1},{2}}, {{0,2},{1}}, {{1,2},{0}}, {{0},{1},{2}}};
    std::sort(partition.begin(), partition.end());
    for (auto& blocks : reference_partition) {
      std::sort(blocks.begin(), blocks.end());
    }
    std::sort(reference_partition.begin(), reference_partition.end());
    EXPECT_EQ(partition, reference_partition);
  }

}