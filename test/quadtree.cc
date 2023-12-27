
#include "spatial-tree/quadtree.hh"

#include <gtest/gtest.h>

#include <random>
#include <string>
#include <unordered_set>

#include "test/fixture.hh"

/// TODO: Try different values.
/// TODO: Try rectangular domains.
/// TODO: Better error messages.
/// TODO: Better edges cases (small example with now powers of two)
TEST(TestQuadTree, UniqueInsertions) {
    const std::vector<int> test_sizes = {1000, 10000, 12345, 13, 256, 2048, 5, 999};
    for (int test_size : test_sizes) {
        st::QuadTree<int, int> tree({-test_size, test_size, test_size, -test_size});
        for (int i = 1; i <= test_size; ++i) {
            ASSERT_TRUE(tree.emplace(i, i, i).second);
            ASSERT_TRUE(tree.emplace(i, -i, i).second);
            ASSERT_TRUE(tree.emplace(-i, i, i).second);
            ASSERT_TRUE(tree.emplace(-i, -i, i).second);
        }
        ASSERT_EQ(tree.size(), 4 * test_size);
        tree.clear();
        ASSERT_TRUE(tree.empty());
        ASSERT_EQ(tree.size(), std::distance(tree.begin(), tree.end()));
    }
}

TEST(TestQuadTree, RandomInsertions) {
    const std::vector<int> test_sizes = {1000, 10000, 12345, 13, 256, 2048, 5, 999};
    const auto             test = [&]<typename CoordT>() {
        for (int test_size : test_sizes) {
            static constexpr int DISTRIBUTION_BEG = 1;
            static constexpr int DISTRIBUTION_END = 10000;
            static constexpr int RANGE_SIZE = DISTRIBUTION_END - DISTRIBUTION_BEG + 1;
            static_assert(RANGE_SIZE * RANGE_SIZE <= std::numeric_limits<int>::max());

            st::QuadTree<int, CoordT> tree(
                {DISTRIBUTION_BEG, DISTRIBUTION_END, DISTRIBUTION_END, DISTRIBUTION_BEG});
            std::random_device            device;
            std::uniform_int_distribution distribution(static_cast<int>(DISTRIBUTION_BEG),
                                                       static_cast<int>(DISTRIBUTION_END));
            std::unordered_set<int>       added_points;
            for (int i = 0; i < test_size; ++i) {
                const int     x = distribution(device);
                const int     y = distribution(device);
                const bool    was_inserted_in_tree = tree.emplace(x, y, i).second;
                const int64_t unique_hash = x + RANGE_SIZE * y;
                const bool    was_inserted_in_map = added_points.insert(unique_hash).second;
                ASSERT_EQ(was_inserted_in_map, was_inserted_in_tree);
            }

            ASSERT_EQ(added_points.size(), tree.size());
            ASSERT_EQ(tree.size(), std::distance(tree.begin(), tree.end()));
        }
    };

    /// TODO: Proper floating point tests.
    test.operator()<float>();
    test.operator()<double>();
    test.operator()<int>();
    test.operator()<unsigned>();
}

TEST(TestQuadTree, IteratorCoverage) {
    const std::vector<int> test_sizes = {1000, 10000, 12345, 13, 256, 2048, 5, 999};
    for (int test_size : test_sizes) {
        st::QuadTree<int, int> tree({-test_size, test_size, test_size, -test_size});
        ASSERT_FALSE(std::distance(tree.begin(), tree.end()));

        for (int i = 0; i < test_size; ++i) {
            ASSERT_TRUE(tree.emplace(i, i, i).second);
        }
        ASSERT_EQ(tree.size(), test_size);
        ASSERT_EQ(std::distance(tree.begin(), tree.end()), test_size);
    }
}

/// TODO: Random test
TEST(TestQuadTree, SingleFind) {
    const std::vector<int> test_sizes = {10};
    for (int test_size : test_sizes) {
        st::QuadTree<int, int> tree({-test_size, test_size, test_size, -test_size});
        ASSERT_EQ(tree.find(0, 0), tree.end());

        for (int i = 0; i < test_size; ++i) {
            ASSERT_TRUE(tree.emplace(i, i, i).second);
        }
        ASSERT_EQ(tree.size(), test_size);

        for (int i = 0; i < test_size; ++i) {
            ASSERT_NE(tree.find(i, i), tree.end());
        }
    }
}

TEST(TestQuadTree, BBoxFind) {
    const std::vector<int> test_sizes = {1000, 10000, 12345, 13, 256, 2048, 5, 999};
    const auto             test = [&]<typename CoordT>() {
        for (int test_size : test_sizes) {
            static constexpr int DISTRIBUTION_BEG = 1;
            static constexpr int DISTRIBUTION_END = 10000;
            static constexpr int RANGE_SIZE = DISTRIBUTION_END - DISTRIBUTION_BEG + 1;
            static_assert(RANGE_SIZE * RANGE_SIZE <= std::numeric_limits<int>::max());

            st::QuadTree<int, CoordT> tree(
                {DISTRIBUTION_BEG, DISTRIBUTION_END, DISTRIBUTION_END, DISTRIBUTION_BEG});
            std::random_device            device;
            std::uniform_int_distribution distribution(static_cast<int>(DISTRIBUTION_BEG),
                                                       static_cast<int>(DISTRIBUTION_END));
            std::unordered_set<int>       added_points;
            for (int i = 0; i < test_size; ++i) {
                const int x = distribution(device);
                const int y = distribution(device);
                tree.emplace(x, y, i);
            }

            static constexpr int N_BBOX = 100;
            for (int i = 0; i < N_BBOX; ++i) {
                const CoordT x_start = distribution(device);
                const CoordT y_start = distribution(device);

                std::uniform_int_distribution distribution_stop_x(int(x_start), DISTRIBUTION_END);
                std::uniform_int_distribution distribution_stop_y(int(y_start), DISTRIBUTION_END);

                const CoordT x_stop = distribution_stop_x(device);
                const CoordT y_stop = distribution_stop_y(device);

                const st::BoundingBox<CoordT> bbox = {x_start, y_stop, x_stop, y_start};

                int counter = std::count_if(tree.begin(), tree.end(), [&](auto it) {
                    const auto [x, y, val] = it;
                    return st::is_inside_bounding_box(x, y, bbox);
                });

                int rcounter = 0;
                tree.find(bbox, [&](auto) { ++rcounter; });

                ASSERT_EQ(counter, rcounter);
            }
        }
    };

    /// TODO: Proper floating point tests.
    test.operator()<float>();
    test.operator()<double>();
    test.operator()<int>();
    test.operator()<unsigned>();
}

TEST(TestQuadTree, Nearest) {
    const std::vector<int> test_sizes = {1000, 10000, 12345, 13, 256, 2048, 5, 999};
    const auto             test = [&]<typename CoordT>() {
        for (int test_size : test_sizes) {
            static constexpr int DISTRIBUTION_BEG = 1;
            static constexpr int DISTRIBUTION_END = 10000;
            static constexpr int RANGE_SIZE = DISTRIBUTION_END - DISTRIBUTION_BEG + 1;
            static_assert(RANGE_SIZE * RANGE_SIZE <= std::numeric_limits<int>::max());

            st::QuadTree<int, CoordT> tree(
                {DISTRIBUTION_BEG, DISTRIBUTION_END, DISTRIBUTION_END, DISTRIBUTION_BEG});
            std::random_device            device;
            std::uniform_int_distribution distribution(static_cast<int>(DISTRIBUTION_BEG),
                                                       static_cast<int>(DISTRIBUTION_END));
            std::unordered_set<int>       added_points;
            for (int i = 0; i < test_size; ++i) {
                const int x = distribution(device);
                const int y = distribution(device);
                tree.emplace(x, y, i);

                const CoordT x_ = distribution(device);
                const CoordT y_ = distribution(device);

                const auto [nearest_x1, nearest_y1, nearest_val1] =
                    *std::min_element(tree.begin(), tree.end(), [&](auto lhs, auto rhs) {
                        const auto [x1, y1, val1] = lhs;
                        const auto [x2, y2, val2] = rhs;

                        return st::euclidean_distance_squared(x1, x_, y1, y_) <
                               st::euclidean_distance_squared(x2, x_, y2, y_);
                    });
                const CoordT distance1 =
                    st::euclidean_distance_squared(x_, nearest_x1, y_, nearest_y1);

                // TODO: Check that it found them all.
                const auto nearest_points = tree.nearest(x_, y_);
                ASSERT_TRUE(std::all_of(nearest_points.begin(), nearest_points.end(), [&](auto it) {
                    const auto [nearest_x2, nearest_y2, nearest_val2] = *it;
                    const CoordT distance2 =
                        st::euclidean_distance_squared(x_, nearest_x2, y_, nearest_y2);
                    return distance2 == distance1;
                }));
            }
        }
    };

    /// TODO: Proper floating point tests.
    test.operator()<float>();
    test.operator()<double>();
    test.operator()<int>();
    test.operator()<unsigned>();
}