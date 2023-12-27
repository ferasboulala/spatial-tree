#include "spatial-tree/euclide.hh"

#include <gtest/gtest.h>

TEST(TestEuclide, TestEuclideanDistanceDouble) {
    using namespace st;
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(0.0, 0.0), 0.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(0.0, 1.0), 1.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(0.0, 2.0), 4.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(0.0, 0.0, 0.0, 0.0), 0.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(0.0, 0.0, 0.0, 1.0), 1.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(0.0, 3.0, 0.0, 4.0), 25.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(1.0, 4.0, 10.0, 14.0), 25.0);
}

TEST(TestEuclide, TestEuclideanDistanceInt) {
    using namespace st;
    ASSERT_EQ(euclidean_distance_squared<int>(0, 0), 0);
    ASSERT_EQ(euclidean_distance_squared<int>(0, 1), 1);
    ASSERT_EQ(euclidean_distance_squared<int>(0, 2), 4);
    ASSERT_EQ(euclidean_distance_squared<int>(0, 0, 0, 0), 0);
    ASSERT_EQ(euclidean_distance_squared<int>(0, 0, 0, 1), 1);
    ASSERT_EQ(euclidean_distance_squared<int>(0, 3, 0, 4), 25);
    ASSERT_EQ(euclidean_distance_squared<int>(1, 4, 10, 14), 25);
}

TEST(TestEuclide, TestEuclideanDistanceUnsigned) {
    using namespace st;
    ASSERT_EQ(euclidean_distance_squared<unsigned>(0, 0), 0);
    ASSERT_EQ(euclidean_distance_squared<unsigned>(0, 1), 1);
    ASSERT_EQ(euclidean_distance_squared<unsigned>(0, 2), 4);
    ASSERT_EQ(euclidean_distance_squared<unsigned>(0, 0, 0, 0), 0);
    ASSERT_EQ(euclidean_distance_squared<unsigned>(0, 0, 0, 1), 1);
    ASSERT_EQ(euclidean_distance_squared<unsigned>(0, 3, 0, 4), 25);
    ASSERT_EQ(euclidean_distance_squared<unsigned>(1, 4, 10, 14), 25);
}

TEST(TestEuclide, TestAbsDiff) {
    using namespace st;
    ASSERT_EQ(st::absdiff<unsigned>(0, 2), 2);
    ASSERT_EQ(st::absdiff<unsigned>(2, 0), 2);
    ASSERT_EQ(st::absdiff<int>(0, 2), 2);
    ASSERT_EQ(st::absdiff<int>(0, -2), 2);
    ASSERT_EQ(st::absdiff<int>(-2, 0), 2);
    ASSERT_EQ(st::absdiff<double>(0, 2), 2);
    ASSERT_EQ(st::absdiff<double>(0, -2), 2);
    ASSERT_EQ(st::absdiff<double>(-2, 0), 2);
}