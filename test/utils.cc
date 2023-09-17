#include <gtest/gtest.h>

#include "spatial-tree/euclide.hh"
#include "spatial-tree/intervals.hh"

TEST(TestUtils, TestEuclideanDistanceDouble) {
    using namespace st;
    ASSERT_DOUBLE_EQ(euclidean_distance_squared(0.0, 0.0), 0.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared(0.0, 1.0), 1.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared(0.0, 2.0), 4.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared(0.0, 0.0, 0.0, 0.0), 0.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared(0.0, 0.0, 0.0, 1.0), 1.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared(0.0, 3.0, 0.0, 4.0), 25.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared(1.0, 4.0, 10.0, 14.0), 25.0);
}

TEST(TestUtils, TestEuclideanDistanceInt) {
    using namespace st;
    ASSERT_EQ(euclidean_distance_squared(0, 0), 0);
    ASSERT_EQ(euclidean_distance_squared(0, 1), 1);
    ASSERT_EQ(euclidean_distance_squared(0, 2), 4);
    ASSERT_EQ(euclidean_distance_squared(0, 0, 0, 0), 0);
    ASSERT_EQ(euclidean_distance_squared(0, 0, 0, 1), 1);
    ASSERT_EQ(euclidean_distance_squared(0, 3, 0, 4), 25);
    ASSERT_EQ(euclidean_distance_squared(1, 4, 10, 14), 25);
}