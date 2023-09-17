#include "spatial-tree/kdtree.hh"

#include <gtest/gtest.h>

#include "test/fixture.hh"

using TestIntKDTree = TreeTest<st::KDTree<int, int>>;
using TestDoubleKDTree = TreeTest<st::KDTree<int, double>>;

#define TREE_TEST(TreeType)                                    \
    TEST_F(TreeType, Capacity) { test_capacity(); }            \
    TEST_F(TreeType, Insertions) {                             \
        test_unique_insertions();                              \
        test_random_insertions();                              \
    }                                                          \
    TEST_F(TreeType, Iterators) { test_iterators_coverage(); } \
    TEST_F(TreeType, Walk) { test_walk(); }                    \
    TEST_F(TreeType, Find) {                                   \
        test_single_find();                                    \
        test_bbox_find();                                      \
    }                                                          \
    TEST_F(TreeType, Nearest) { test_nearest(); }

TREE_TEST(TestIntKDTree);
TREE_TEST(TestDoubleKDTree);