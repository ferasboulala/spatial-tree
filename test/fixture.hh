#pragma once

#include <gtest/gtest.h>

#include <random>
#include <unordered_set>

template <class TreeType, typename... Args>
class TreeTest : public ::testing::Test {
protected:
    void test_capacity(Args &&...args) {
        TreeType                  tree(std::forward<Args>(args)...);
        static constexpr uint64_t TEST_SIZE = 10000;
        tree.reserve(TEST_SIZE);
        ASSERT_EQ(tree.capacity(), TEST_SIZE);
        ASSERT_TRUE(tree.empty());
    }

    void test_unique_insertions(Args &&...args) {
        TreeType             tree(std::forward<Args>(args)...);
        static constexpr int TEST_SIZE = 10000;
        for (int i = 1; i <= TEST_SIZE; ++i) {
            ASSERT_TRUE(tree.emplace(i, i, i).second);
            ASSERT_TRUE(tree.emplace(i, -i, i).second);
            ASSERT_TRUE(tree.emplace(-i, i, i).second);
            ASSERT_TRUE(tree.emplace(-i, -i, i).second);
        }
        ASSERT_EQ(tree.size(), 4 * TEST_SIZE);
        tree.clear();
        ASSERT_TRUE(tree.empty());
        ASSERT_EQ(tree.size(), std::distance(tree.begin(), tree.end()));
    }

    void test_random_insertions(Args &&...args) {
        using CoordT = typename TreeType::__CoordinateType;
        if constexpr (!std::is_integral_v<CoordT>) {
            GTEST_SKIP();
        } else {
            static constexpr int    TEST_SIZE = 10000;
            static constexpr CoordT DISTRIBUTION_BEG = -100;
            static constexpr CoordT DISTRIBUTION_END = 100;
            static_assert((DISTRIBUTION_END - DISTRIBUTION_BEG + 1) *
                              (DISTRIBUTION_END - DISTRIBUTION_BEG + 1) <=
                          std::numeric_limits<CoordT>::max());

            st::KDTree<int>               tree;
            std::random_device            device;
            std::uniform_int_distribution distribution(DISTRIBUTION_BEG, DISTRIBUTION_END);
            std::unordered_set<CoordT>    added_points;
            for (int i = 0; i < TEST_SIZE; ++i) {
                const CoordT x = distribution(device);
                const CoordT y = distribution(device);
                const bool   was_inserted_in_tree = tree.emplace(x, y, i).second;
                const CoordT unique_hash = x + (DISTRIBUTION_END - DISTRIBUTION_BEG + 1) * y;
                const bool   was_inserted_in_map = added_points.insert(unique_hash).second;
                ASSERT_EQ(was_inserted_in_map, was_inserted_in_tree);
            }

            ASSERT_EQ(added_points.size(), tree.size());
            ASSERT_EQ(tree.size(), std::distance(tree.begin(), tree.end()));
        }
    }

    void test_iterators_coverage(Args &&...args) {
        TreeType tree(std::forward<Args>(args)...);
        ASSERT_FALSE(std::distance(tree.begin(), tree.end()));

        static constexpr int TEST_SIZE = 10000;
        for (int i = 0; i < TEST_SIZE; ++i) {
            ASSERT_TRUE(tree.emplace(i, i, i).second);
        }
        ASSERT_EQ(tree.size(), TEST_SIZE);
        ASSERT_EQ(std::distance(tree.begin(), tree.end()), TEST_SIZE);
    }

    void test_walk(Args &&...args) {
        TreeType   tree(std::forward<Args>(args)...);
        uint64_t   called = 0;
        const auto callback = [&](auto) { ++called; };
        tree.walk(callback);
        ASSERT_FALSE(called);

        static constexpr int TEST_SIZE = 10000;
        for (int i = 0; i < TEST_SIZE; ++i) {
            ASSERT_TRUE(tree.emplace(i, i, i).second);
        }
        ASSERT_EQ(tree.size(), TEST_SIZE);

        tree.walk(callback);
        ASSERT_EQ(called, TEST_SIZE);
    }

    void test_single_find(Args &&...args) {
        TreeType tree(std::forward<Args>(args)...);
        ASSERT_EQ(tree.find(0, 0), tree.end());

        static constexpr int TEST_SIZE = 10000;
        for (int i = 0; i < TEST_SIZE; ++i) {
            ASSERT_TRUE(tree.emplace(i, i, i).second);
        }
        ASSERT_EQ(tree.size(), TEST_SIZE);

        for (int i = 0; i < TEST_SIZE; ++i) {
            ASSERT_NE(tree.find(i, i), tree.end());
        }
    }

    void test_bbox_find(Args &&...args) {
        using CoordT = typename TreeType::__CoordinateType;
        if constexpr (!std::is_integral_v<CoordT>) {
            GTEST_SKIP();
        } else {
            static constexpr int    TEST_SIZE = 10000;
            static constexpr CoordT DISTRIBUTION_BEG = -100;
            static constexpr CoordT DISTRIBUTION_END = 100;
            static_assert((DISTRIBUTION_END - DISTRIBUTION_BEG + 1) *
                              (DISTRIBUTION_END - DISTRIBUTION_BEG + 1) <=
                          std::numeric_limits<CoordT>::max());

            st::KDTree<int>               tree;
            std::random_device            device;
            std::uniform_int_distribution distribution(DISTRIBUTION_BEG, DISTRIBUTION_END);
            for (int i = 0; i < TEST_SIZE; ++i) {
                const CoordT x = distribution(device);
                const CoordT y = distribution(device);
                tree.emplace(x, y, i);
            }

            static constexpr int N_BBOX = 100;
            for (int i = 0; i < N_BBOX; ++i) {
                const CoordT x_start = distribution(device);
                const CoordT y_start = distribution(device);

                std::uniform_int_distribution distribution_stop_x(x_start, DISTRIBUTION_END);
                std::uniform_int_distribution distribution_stop_y(y_start, DISTRIBUTION_END);

                const CoordT x_stop = distribution_stop_x(device);
                const CoordT y_stop = distribution_stop_y(device);

                const st::BoundingBox bbox = {x_start, y_stop, x_stop, y_start};

                int counter = std::count_if(tree.begin(), tree.end(), [&](auto it) {
                    const auto [x, y, val] = it;
                    return st::is_inside_bounding_box(x, y, bbox);
                });

                tree.find(bbox, [&](auto) { --counter; });

                ASSERT_EQ(counter, 0);
            }
        }
    }

    void test_nearest(Args &&...args) {
        using CoordT = typename TreeType::__CoordinateType;
        if constexpr (!std::is_integral_v<CoordT>) {
            GTEST_SKIP();
        } else {
            static constexpr int    TEST_SIZE = 10000;
            static constexpr CoordT DISTRIBUTION_BEG = -100;
            static constexpr CoordT DISTRIBUTION_END = 100;
            static_assert((DISTRIBUTION_END - DISTRIBUTION_BEG + 1) *
                              (DISTRIBUTION_END - DISTRIBUTION_BEG + 1) <=
                          std::numeric_limits<CoordT>::max());

            st::KDTree<int> tree;
            ASSERT_TRUE(tree.nearest(0, 0).empty());

            std::random_device            device;
            std::uniform_int_distribution distribution(DISTRIBUTION_BEG, DISTRIBUTION_END);
            for (int i = 0; i < TEST_SIZE; ++i) {
                const CoordT x = distribution(device);
                const CoordT y = distribution(device);
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

                const auto nearest_points = tree.nearest(x_, y_);
                ASSERT_TRUE(std::all_of(nearest_points.begin(), nearest_points.end(), [&](auto it) {
                    const auto [nearest_x2, nearest_y2, nearest_val2] = *it;
                    const CoordT distance2 =
                        st::euclidean_distance_squared(x_, nearest_x2, y_, nearest_y2);
                    return distance2 == distance1;
                }));
            }
        }
    }
};