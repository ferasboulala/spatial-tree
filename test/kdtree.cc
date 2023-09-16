#include "spatial-tree/kdtree.hh"

#include <gtest/gtest.h>

#include <random>
#include <unordered_set>
#include <vector>

#include "test/fixture.hh"

using TestIntKDTree = TreeTest<st::KDTree<int, int>>;
using TestDoubleKDTree = TreeTest<st::KDTree<int, double>>;

#define TREE_TEST(TreeType)        \
    TEST_F(TreeType, Test) {       \
        test_capacity();           \
        test_unique_insertions();  \
        test_random_insertions();  \
        test_iterators_coverage(); \
        test_walk();               \
        test_single_find();        \
        test_bbox_find();          \
        test_nearest();            \
    }

TREE_TEST(TestIntKDTree);
TREE_TEST(TestDoubleKDTree);