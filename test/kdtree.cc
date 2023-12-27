#include <gtest/gtest.h>

#ifdef NDEBUG
#undef NDEBUG
#include "spatial-tree/kdtree.hh"
#else
#include "spatial-tree/kdtree.hh"
#endif

#include <gtest/gtest.h>

#include <random>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

TEST(TestKDTree, ManyValueTypes) {
    st::KDTree<int>         integer;
    st::KDTree<float>       floating;
    st::KDTree<std::string> str;
    st::KDTree<std::vector<std::string>> vec;
    /// TODO
    // st::KDTree<void> nothing;
}

TEST(TestKDTree, UniqueInsertions) {
    const auto typed_test = [&]<typename CoordT>() {
        const std::vector<uint64_t> test_sizes = {1, 2, 3, 7, 10, 100, 9999, 10000};
        for (uint64_t test_size : test_sizes) {
            const auto inner_test = [&]<typename TreeType>(auto &tree) {
                for (uint64_t i = 1; i <= test_size; ++i) {
                    CoordT coord = i;
                    ASSERT_TRUE(tree.emplace(coord, coord, i).second);
                    ASSERT_TRUE(tree.emplace(coord, -coord, i).second);
                    ASSERT_TRUE(tree.emplace(-coord, coord, i).second);
                    ASSERT_TRUE(tree.emplace(-coord, -coord, i).second);
                }
                ASSERT_EQ(tree.size(), 4 * test_size);
                ASSERT_EQ(tree.size(), std::distance(tree.begin(), tree.end()));
                tree.clear();
            };
            st::KDTree<int, CoordT> tree;
            inner_test.template     operator()<st::KDTree<int, CoordT>>(tree);
        }
    };
    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestKDTree, RandomInsertions) {
    static constexpr int DISTRIBUTION_BEG = -10000;
    static constexpr int DISTRIBUTION_END = 10000;
    static constexpr int RANGE_SIZE = DISTRIBUTION_END - DISTRIBUTION_BEG + 1;
    static_assert(RANGE_SIZE * RANGE_SIZE <= std::numeric_limits<int>::max());
    std::random_device            device;
    std::uniform_int_distribution distribution(static_cast<int>(DISTRIBUTION_BEG),
                                               static_cast<int>(DISTRIBUTION_END));

    const std::vector<uint64_t> test_sizes = {1, 2, 3, 7, 10, 100, 9999, 10000};
    const auto                  typed_test = [&]<typename CoordT>() {
        const auto inner_test = [&]<typename TreeType>(auto &tree) {
            for (uint64_t test_size : test_sizes) {
                std::unordered_set<int> added_points;
                for (uint64_t i = 0; i < test_size; ++i) {
                    const int     x = distribution(device);
                    const int     y = distribution(device);
                    const bool    was_inserted_in_tree = tree.emplace(x, y, i).second;
                    const int64_t unique_hash = x + RANGE_SIZE * y;
                    const bool    was_inserted_in_map = added_points.insert(unique_hash).second;
                    ASSERT_EQ(was_inserted_in_map, was_inserted_in_tree);
                }

                ASSERT_EQ(added_points.size(), tree.size());
                ASSERT_EQ(tree.size(), std::distance(tree.begin(), tree.end()));
                tree.clear();
            }
        };
        st::KDTree<int, CoordT> tree;
        inner_test.template     operator()<st::KDTree<int, CoordT>>(tree);
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestKDTree, IteratorCoverage) {
    const std::vector<uint64_t> test_sizes = {1, 2, 3, 7, 10, 100, 9999, 10000};
    const auto                  typed_test = [&]<typename CoordT>() {
        for (int test_size : test_sizes) {
            const auto inner_test = [&]<typename TreeType>(auto &tree) {
                ASSERT_FALSE(std::distance(tree.begin(), tree.end()));

                for (int i = 0; i < test_size; ++i) {
                    ASSERT_TRUE(tree.emplace(i, i, i).second);
                }
                ASSERT_EQ(tree.size(), test_size);
                ASSERT_EQ(std::distance(tree.begin(), tree.end()), test_size);
            };
            st::KDTree<int, CoordT> tree;
            inner_test.template     operator()<st::KDTree<int, CoordT>>(tree);
        }
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestKDTree, SingleFind) {
    const std::vector<uint64_t> test_sizes = {1, 2, 3, 7, 10, 100, 9999, 10000};
    const auto                  typed_test = [&]<typename CoordT>() {
        for (int test_size : test_sizes) {
            const auto inner_test = [&]<typename TreeType>(auto &tree) {
                ASSERT_EQ(tree.find(0, 0), tree.end());

                for (int i = 0; i < test_size; ++i) {
                    ASSERT_TRUE(tree.emplace(i, i, i).second);
                }
                ASSERT_EQ(tree.size(), test_size);

                for (int i = 0; i < test_size; ++i) {
                    ASSERT_NE(tree.find(i, i), tree.end());
                }
            };
            st::KDTree<int, CoordT> tree;
            inner_test.template     operator()<st::KDTree<int, CoordT>>(tree);
        }
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestKDTree, BBoxFind) {
    static constexpr int DISTRIBUTION_BEG = -10000;
    static constexpr int DISTRIBUTION_END = 10000;
    static constexpr int RANGE_SIZE = DISTRIBUTION_END - DISTRIBUTION_BEG + 1;
    static_assert(RANGE_SIZE * RANGE_SIZE <= std::numeric_limits<int>::max());
    std::random_device            device;
    std::uniform_int_distribution distribution(static_cast<int>(DISTRIBUTION_BEG),
                                               static_cast<int>(DISTRIBUTION_END));

    const std::vector<uint64_t> test_sizes = {1, 2, 3, 7, 10, 100, 9999, 10000};
    const auto                  typed_test = [&]<typename CoordT>() {
        const auto inner_test = [&]<typename TreeType>(auto &tree) {
            for (uint64_t test_size : test_sizes) {
                std::unordered_set<int> added_points;
                for (uint64_t i = 0; i < test_size; ++i) {
                    const int x = distribution(device);
                    const int y = distribution(device);
                    tree.emplace(x, y, i);
                }

                static constexpr int N_BBOX = 100;
                for (int i = 0; i < N_BBOX; ++i) {
                    const CoordT x_start = distribution(device);
                    const CoordT y_start = distribution(device);

                    std::uniform_int_distribution distribution_stop_x(int(x_start),
                                                                                       DISTRIBUTION_END);
                    std::uniform_int_distribution distribution_stop_y(int(y_start),
                                                                                       DISTRIBUTION_END);

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
        st::KDTree<int, CoordT> tree;
        inner_test.template     operator()<st::KDTree<int, CoordT>>(tree);
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestKDTree, Nearest) {
    static constexpr int DISTRIBUTION_BEG = -10000;
    static constexpr int DISTRIBUTION_END = 10000;
    static constexpr int RANGE_SIZE = DISTRIBUTION_END - DISTRIBUTION_BEG + 1;
    static_assert(RANGE_SIZE * RANGE_SIZE <= std::numeric_limits<int>::max());
    std::random_device            device;
    std::uniform_int_distribution distribution(static_cast<int>(DISTRIBUTION_BEG),
                                               static_cast<int>(DISTRIBUTION_END));

    const std::vector<uint64_t> test_sizes = {1, 2, 3, 7, 10, 100, 9999, 10000};
    const auto                  typed_test = [&]<typename CoordT>() {
        const auto inner_test = [&]<typename TreeType>(auto &tree) {
            for (uint64_t test_size : test_sizes) {
                std::unordered_set<int> added_points;
                for (uint64_t i = 0; i < test_size; ++i) {
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
                    ASSERT_TRUE(
                        std::all_of(nearest_points.begin(), nearest_points.end(), [&](auto it) {
                            const auto [nearest_x2, nearest_y2, nearest_val2] = *it;
                            const CoordT distance2 =
                                st::euclidean_distance_squared(x_, nearest_x2, y_, nearest_y2);
                            return distance2 == distance1;
                        }));
                }
            }
        };
        st::KDTree<int, CoordT> tree;
        inner_test.template     operator()<st::KDTree<int, CoordT>>(tree);
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}