
#include "spatial-tree/quadtree.hh"

#include <gtest/gtest.h>

#include <random>
#include <string>
#include <unordered_set>

#include "test/fixture.hh"

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
        // ASSERT_EQ(tree.size(), std::distance(tree.begin(), tree.end()));
    }
}

TEST(TestQuadTree, RandomInsertions) {
    const auto test = [&]<typename CoordT>() {
        static constexpr int    TEST_SIZE = 10000;
        static constexpr CoordT DISTRIBUTION_BEG = 0;
        static constexpr CoordT DISTRIBUTION_END = 1000;
        static_assert((DISTRIBUTION_END - DISTRIBUTION_BEG + 1) *
                          (DISTRIBUTION_END - DISTRIBUTION_BEG + 1) <=
                      std::numeric_limits<CoordT>::max());

        st::QuadTree<int, CoordT> tree(
            {DISTRIBUTION_BEG, DISTRIBUTION_END, DISTRIBUTION_END, DISTRIBUTION_BEG});
        std::random_device            device;
        std::uniform_int_distribution distribution(static_cast<int>(DISTRIBUTION_BEG),
                                                   static_cast<int>(DISTRIBUTION_END));
        std::unordered_set<int>       added_points;
        for (int i = 0; i < TEST_SIZE; ++i) {
            const CoordT x = distribution(device);
            const CoordT y = distribution(device);
            const bool   was_inserted_in_tree = tree.emplace(x, y, i).second;
            const CoordT unique_hash = x + (DISTRIBUTION_END - DISTRIBUTION_BEG + 1) * y;
            const bool   was_inserted_in_map = added_points.insert(unique_hash).second;
            ASSERT_EQ(was_inserted_in_map, was_inserted_in_tree);
        }

        ASSERT_EQ(added_points.size(), tree.size());
        // ASSERT_EQ(tree.size(), std::distance(tree.begin(), tree.end()));
    };

    test.operator()<float>();
    test.operator()<double>();
    test.operator()<int>();
    test.operator()<unsigned>();
}

TEST(TestQuadTree, IteratorCoverage) {
    static constexpr int   TEST_SIZE = 10000;
    st::QuadTree<int, int> tree({-TEST_SIZE, TEST_SIZE, TEST_SIZE, -TEST_SIZE});
    // ASSERT_FALSE(std::distance(tree.begin(), tree.end()));

    for (int i = 0; i < TEST_SIZE; ++i) {
        ASSERT_TRUE(tree.emplace(i, i, i).second);
    }
    ASSERT_EQ(tree.size(), TEST_SIZE);
    // ASSERT_EQ(std::distance(tree.begin(), tree.end()), TEST_SIZE);
}

// TEST(TestQuadTree, Walk) {
//     static constexpr int   TEST_SIZE = 10000;
//     st::QuadTree<int, int> tree({-TEST_SIZE, TEST_SIZE, TEST_SIZE, -TEST_SIZE});
//     uint64_t               called = 0;
//     const auto             callback = [&](auto) { ++called; };
//     tree.walk(callback);
//     ASSERT_FALSE(called);

//     for (int i = 0; i < TEST_SIZE; ++i) {
//         ASSERT_TRUE(tree.emplace(i, i, i).second);
//     }
//     ASSERT_EQ(tree.size(), TEST_SIZE);

//     tree.walk(callback);
//     ASSERT_EQ(called, TEST_SIZE);
// }

// TEST(TestQuadTree, SingleFind) {
//     static constexpr int   TEST_SIZE = 10000;
//     st::QuadTree<int, int> tree({-TEST_SIZE, TEST_SIZE, TEST_SIZE, -TEST_SIZE});
//     ASSERT_EQ(tree.find(0, 0), tree.end());

//     for (int i = 0; i < TEST_SIZE; ++i) {
//         ASSERT_TRUE(tree.emplace(i, i, i).second);
//     }
//     ASSERT_EQ(tree.size(), TEST_SIZE);

//     for (int i = 0; i < TEST_SIZE; ++i) {
//         ASSERT_NE(tree.find(i, i), tree.end());
//     }
// }

// TEST(TestQuadTree, BBoxFind) {
//     using CoordT = int;
//     static constexpr int    TEST_SIZE = 10000;
//     static constexpr CoordT DISTRIBUTION_BEG = 0;
//     static constexpr CoordT DISTRIBUTION_END = 100;
//     static_assert((DISTRIBUTION_END - DISTRIBUTION_BEG + 1) *
//                       (DISTRIBUTION_END - DISTRIBUTION_BEG + 1) <=
//                   std::numeric_limits<CoordT>::max());

//     st::QuadTree<int, CoordT> tree(
//         {DISTRIBUTION_BEG, DISTRIBUTION_END, DISTRIBUTION_END, DISTRIBUTION_BEG});
//     std::random_device            device;
//     std::uniform_int_distribution distribution(DISTRIBUTION_BEG, DISTRIBUTION_END);
//     for (int i = 0; i < TEST_SIZE; ++i) {
//         const CoordT x = distribution(device);
//         const CoordT y = distribution(device);
//         tree.emplace(x, y, i);
//     }

//     static constexpr int N_BBOX = 100;
//     for (int i = 0; i < N_BBOX; ++i) {
//         const CoordT x_start = distribution(device);
//         const CoordT y_start = distribution(device);

//         std::uniform_int_distribution distribution_stop_x(x_start, DISTRIBUTION_END);
//         std::uniform_int_distribution distribution_stop_y(y_start, DISTRIBUTION_END);

//         const CoordT x_stop = distribution_stop_x(device);
//         const CoordT y_stop = distribution_stop_y(device);

//         const st::BoundingBox bbox = {x_start, y_stop, x_stop, y_start};

//         int counter = std::count_if(tree.begin(), tree.end(), [&](auto it) {
//             const auto [x, y, val] = it;
//             return st::is_inside_bounding_box(x, y, bbox);
//         });

//         tree.find(bbox, [&](auto) { --counter; });

//         ASSERT_EQ(counter, 0);
//     }
// }

// TEST(TestQuadTree, Nearest) {
//     using CoordT = int;
//     if constexpr (!std::is_integral_v<CoordT>) {
//         GTEST_SKIP();
//     } else {
//         static constexpr int    TEST_SIZE = 10000;
//         static constexpr CoordT DISTRIBUTION_BEG = -100;
//         static constexpr CoordT DISTRIBUTION_END = 100;
//         static_assert((DISTRIBUTION_END - DISTRIBUTION_BEG + 1) *
//                           (DISTRIBUTION_END - DISTRIBUTION_BEG + 1) <=
//                       std::numeric_limits<CoordT>::max());

//         st::QuadTree<int, CoordT> tree(
//             {DISTRIBUTION_BEG, DISTRIBUTION_END, DISTRIBUTION_END, DISTRIBUTION_BEG});
//         ASSERT_TRUE(tree.nearest(0, 0).empty());

//         std::random_device            device;
//         std::uniform_int_distribution distribution(DISTRIBUTION_BEG, DISTRIBUTION_END);
//         for (int i = 0; i < TEST_SIZE; ++i) {
//             const CoordT x = distribution(device);
//             const CoordT y = distribution(device);
//             tree.emplace(x, y, i);

//             const CoordT x_ = distribution(device);
//             const CoordT y_ = distribution(device);

//             const auto [nearest_x1, nearest_y1, nearest_val1] =
//                 *std::min_element(tree.begin(), tree.end(), [&](auto lhs, auto rhs) {
//                     const auto [x1, y1, val1] = lhs;
//                     const auto [x2, y2, val2] = rhs;

//                     return st::euclidean_distance_squared(x1, x_, y1, y_) <
//                            st::euclidean_distance_squared(x2, x_, y2, y_);
//                 });
//             const CoordT distance1 = st::euclidean_distance_squared(x_, nearest_x1, y_,
//             nearest_y1);

//             const auto nearest_points = tree.nearest(x_, y_);
//             ASSERT_TRUE(std::all_of(nearest_points.begin(), nearest_points.end(), [&](auto it) {
//                 const auto [nearest_x2, nearest_y2, nearest_val2] = *it;
//                 const CoordT distance2 =
//                     st::euclidean_distance_squared(x_, nearest_x2, y_, nearest_y2);
//                 return distance2 == distance1;
//             }));
//         }
//     }
// }