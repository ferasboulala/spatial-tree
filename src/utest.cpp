#ifdef NDEBUG
#undef NDEBUG
#include "spatial-tree.h"
#define NDEBUG
#else
#include "spatial-tree.h"
#endif

#include <gtest/gtest.h>

#include <random>
#include <string>
#include <unordered_set>
#include <vector>

TEST(TestUnroll, TestUnroll) {
    const auto typed_test = [&]<typename CoordT>() {
        CoordT start = 1;
        CoordT stop = 5;
        CoordT sum = 0;

        st::internal::unroll_for<1>(start, stop, [&](CoordT i) { sum += i; });
        ASSERT_EQ(sum, 10);

        sum = 0;
        st::internal::unroll_for<2>(start, stop, [&](CoordT i) { sum += i; });
        ASSERT_EQ(sum, 10);

        sum = 0;
        st::internal::unroll_for<3>(start, stop, [&](CoordT i) { sum += i; });
        ASSERT_EQ(sum, 10);

        sum = 0;
        st::internal::unroll_for<4>(start, stop, [&](CoordT i) { sum += i; });
        ASSERT_EQ(sum, 10);
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestBoundingBox, TestDefaultConstructor) {
    const auto typed_test = [&]<typename CoordT>() {
        const st::bounding_box<CoordT, 3> bbox;
        ASSERT_EQ(bbox.starts.size(), CoordT(3));
        ASSERT_EQ(bbox.stops.size(), CoordT(3));
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestBoundingBox, TestConstructor) {
    const auto typed_test = [&]<typename CoordT>() {
        const st::bounding_box<CoordT, 3> bbox = {0, 0, 0, 3, 4, 5};
        ASSERT_EQ(bbox.starts.size(), CoordT(3));
        ASSERT_EQ(bbox.stops.size(), CoordT(3));
        ASSERT_EQ(bbox.area(), CoordT(3 * 4 * 5));
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestBoundingBox, Contains) {
    const auto typed_test = [&]<typename CoordT>() {
        const st::bounding_box<CoordT, 1> vec = {0, 10};
        ASSERT_TRUE(vec.contains({0}));
        ASSERT_TRUE(vec.contains({5}));
        ASSERT_TRUE(vec.contains({10}));
        ASSERT_FALSE(vec.contains({11}));

        const st::bounding_box<CoordT, 2> mat = {0, 0, 10, 10};
        ASSERT_TRUE(mat.contains({0, 0}));
        ASSERT_TRUE(mat.contains({0, 5}));
        ASSERT_TRUE(mat.contains({5, 0}));
        ASSERT_TRUE(mat.contains({10, 10}));
        ASSERT_FALSE(mat.contains({11, 0}));
        ASSERT_FALSE(mat.contains({0, 11}));
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestBoundingBox, Overlaps) {
    const auto typed_test = [&]<typename CoordT>() {
        const st::bounding_box<CoordT, 2> lhs = {0, 0, 10, 10};
        ASSERT_TRUE(lhs.overlaps(lhs));

        for (CoordT i = -10; i < -5; ++i) {
            const st::bounding_box<CoordT, 2> rhs = {i, i, i + 5, i + 5};
            ASSERT_TRUE(rhs.overlaps(rhs));
            ASSERT_FALSE(lhs.overlaps(rhs));
            ASSERT_FALSE(rhs.overlaps(lhs));
        }

        for (CoordT i = -5; i < 10; ++i) {
            const st::bounding_box<CoordT, 2> rhs = {i, i, i + 5, i + 5};
            ASSERT_TRUE(rhs.overlaps(rhs));
            ASSERT_TRUE(lhs.overlaps(rhs));
            ASSERT_TRUE(rhs.overlaps(lhs));
        }

        for (CoordT i = 11; i < 15; ++i) {
            const st::bounding_box<CoordT, 2> rhs = {i, i, i + 5, i + 5};
            ASSERT_TRUE(rhs.overlaps(rhs));
            ASSERT_FALSE(lhs.overlaps(rhs));
            ASSERT_FALSE(rhs.overlaps(lhs));
        }
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestBoundingBox, BelongsToQuadrant) {
    const auto typed_test = [&]<typename CoordT>() {
        const st::bounding_box<CoordT, 1> bbox1 = {0, 10};
        for (CoordT i = 0; i <= 5; ++i) {
            ASSERT_EQ(bbox1.quadrant({i}), 0);
        }

        for (CoordT i = 6; i <= 10; ++i) {
            ASSERT_EQ(bbox1.quadrant({i}), 1);
        }

        const st::bounding_box<CoordT, 2> bbox2 = {0, 0, 10, 10};
        for (CoordT i = 0; i <= 5; ++i) {
            ASSERT_EQ(bbox2.quadrant({i, i}), 0);
        }

        for (CoordT i = 6; i <= 10; ++i) {
            ASSERT_EQ(bbox2.quadrant({i, i}), 3);
        }
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestBoundingBox, QRecurse) {
    const auto typed_test = [&]<typename CoordT>() {
        const st::bounding_box<CoordT, 1> bbox1 = {0, 10};
        ASSERT_EQ(bbox1.qrecurse(0), (st::bounding_box<CoordT, 1>({0, 5})));
        ASSERT_EQ(bbox1.qrecurse(1), (st::bounding_box<CoordT, 1>({5, 10})));

        const st::bounding_box<CoordT, 2> bbox2 = {0, 0, 10, 10};
        ASSERT_EQ(bbox2.qrecurse(0), (st::bounding_box<CoordT, 2>({0, 0, 5, 5})));
        ASSERT_EQ(bbox2.qrecurse(1), (st::bounding_box<CoordT, 2>({5, 0, 10, 5})));
        ASSERT_EQ(bbox2.qrecurse(2), (st::bounding_box<CoordT, 2>({0, 5, 5, 10})));
        ASSERT_EQ(bbox2.qrecurse(3), (st::bounding_box<CoordT, 2>({5, 5, 10, 10})));
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestBoundingBox, Recurse) {
    const auto typed_test = [&]<typename CoordT>() {
        const st::bounding_box<CoordT, 1> bbox1 = {0, 10};
        ASSERT_EQ(std::get<0>(bbox1.recurse({2})), (st::bounding_box<CoordT, 1>{0, 5}));
        ASSERT_EQ(std::get<0>(bbox1.recurse({7})), (st::bounding_box<CoordT, 1>{5, 10}));

        ASSERT_EQ(std::get<1>(bbox1.recurse({2})), 0);
        ASSERT_EQ(std::get<1>(bbox1.recurse({7})), 1);

        const st::bounding_box<CoordT, 2> bbox2 = {0, 0, 10, 10};
        ASSERT_EQ(std::get<0>(bbox2.recurse({2, 2})), (st::bounding_box<CoordT, 2>{0, 0, 5, 5}));
        ASSERT_EQ(std::get<0>(bbox2.recurse({7, 2})), (st::bounding_box<CoordT, 2>{5, 0, 10, 5}));
        ASSERT_EQ(std::get<0>(bbox2.recurse({2, 7})), (st::bounding_box<CoordT, 2>{0, 5, 5, 10}));
        ASSERT_EQ(std::get<0>(bbox2.recurse({7, 7})), (st::bounding_box<CoordT, 2>{5, 5, 10, 10}));

        ASSERT_EQ(std::get<1>(bbox2.recurse({2, 2})), 0);
        ASSERT_EQ(std::get<1>(bbox2.recurse({7, 2})), 1);
        ASSERT_EQ(std::get<1>(bbox2.recurse({2, 7})), 2);
        ASSERT_EQ(std::get<1>(bbox2.recurse({7, 7})), 3);
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

// TODO: Test when a point is inside the bbox.
TEST(TestBoundingBox, SDistance) {
    const auto typed_test = [&]<typename CoordT>() {
        const st::bounding_box<CoordT, 2> bbox = {0, 0, 10, 10};
        ASSERT_EQ(bbox.sdistance({0, 0}), 0);
        ASSERT_EQ(bbox.sdistance({-2, 0}), 4);
        ASSERT_EQ(bbox.sdistance({0, -2}), 4);
        ASSERT_EQ(bbox.sdistance({2, 0}), 0);
        ASSERT_EQ(bbox.sdistance({0, 2}), 0);
        ASSERT_EQ(bbox.sdistance({13, 14}), 25);
        ASSERT_EQ(bbox.sdistance({13, 14}), 25);
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestEuclide, TestEuclideanDistanceDouble) {
    using namespace st::internal;
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(0.0, 0.0), 0.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(0.0, 1.0), 1.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(0.0, 2.0), 4.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(0.0, 0.0, 0.0, 0.0), 0.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(0.0, 0.0, 0.0, 1.0), 1.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(0.0, 3.0, 0.0, 4.0), 25.0);
    ASSERT_DOUBLE_EQ(euclidean_distance_squared<double>(1.0, 4.0, 10.0, 14.0), 25.0);
}

TEST(TestEuclide, TestEuclideanDistanceInt) {
    using namespace st::internal;
    ASSERT_EQ(euclidean_distance_squared<int>(0, 0), 0);
    ASSERT_EQ(euclidean_distance_squared<int>(0, 1), 1);
    ASSERT_EQ(euclidean_distance_squared<int>(0, 2), 4);
    ASSERT_EQ(euclidean_distance_squared<int>(0, 0, 0, 0), 0);
    ASSERT_EQ(euclidean_distance_squared<int>(0, 0, 0, 1), 1);
    ASSERT_EQ(euclidean_distance_squared<int>(0, 3, 0, 4), 25);
    ASSERT_EQ(euclidean_distance_squared<int>(1, 4, 10, 14), 25);
}

TEST(TestEuclide, TestEuclideanDistanceUnsigned) {
    using namespace st::internal;
    ASSERT_EQ(euclidean_distance_squared<unsigned>(0, 0), 0);
    ASSERT_EQ(euclidean_distance_squared<unsigned>(0, 1), 1);
    ASSERT_EQ(euclidean_distance_squared<unsigned>(0, 2), 4);
    ASSERT_EQ(euclidean_distance_squared<unsigned>(0, 0, 0, 0), 0);
    ASSERT_EQ(euclidean_distance_squared<unsigned>(0, 0, 0, 1), 1);
    ASSERT_EQ(euclidean_distance_squared<unsigned>(0, 3, 0, 4), 25);
    ASSERT_EQ(euclidean_distance_squared<unsigned>(1, 4, 10, 14), 25);
}

TEST(TestEuclide, TestAbsDiff) {
    using namespace st::internal;
    ASSERT_EQ(absdiff<unsigned>(0, 2), 2);
    ASSERT_EQ(absdiff<unsigned>(2, 0), 2);
    ASSERT_EQ(absdiff<int>(0, 2), 2);
    ASSERT_EQ(absdiff<int>(0, -2), 2);
    ASSERT_EQ(absdiff<int>(-2, 0), 2);
    ASSERT_EQ(absdiff<double>(0, 2), 2);
    ASSERT_EQ(absdiff<double>(0, -2), 2);
    ASSERT_EQ(absdiff<double>(-2, 0), 2);
}

TEST(TestSpatialTree, ManyValueTypes) {
    st::internal::spatial_tree<float, float, 2>                    floating;
    st::internal::spatial_tree<float, std::string, 2>              str;
    st::internal::spatial_tree<float, std::vector<std::string>, 2> vec;
    st::internal::spatial_tree<float, void, 2>                     nothing;
}

TEST(TestSpatialTree, MutableAccesses) {
    st::internal::spatial_tree<int, int, 2> tree;
    auto                                    it = tree.emplace({0, 0}, 0).first;
    auto [coordinates, data1] = *it;
    ++data1;
    ++tree[{0, 0}];
    auto it2 = tree.find({0, 0});
    auto [_, data2] = *it2;
    ASSERT_EQ(data2, 2);
}

TEST(TestSpatialTree, Destructors) {
    int counter = 0;
    class MyObject {
    public:
        inline MyObject(int &c) : counter(&c) {
            assert(enabled);
            enabled = true;
            *counter = *counter + 1;
        }
        inline MyObject(MyObject &&other) : enabled(other.enabled), counter(other.counter) {
            assert(enabled);
            *counter = *counter + 1;
        }
        inline MyObject(const MyObject &other) : enabled(other.enabled), counter(other.counter) {
            assert(enabled);
            *counter = *counter + 1;
        }
        inline ~MyObject() {
            assert(enabled);
            enabled = false;
            *counter = *counter - 1;
        }

    private:
        bool enabled;
        int *counter;
    };

    st::internal::spatial_tree<int, MyObject, 2> tree;
    ASSERT_TRUE(tree.emplace({0, 0}, counter).second);
    ASSERT_TRUE(tree.emplace({0, 1}, counter).second);
    ASSERT_TRUE(tree.emplace({1, 0}, counter).second);
    ASSERT_TRUE(tree.emplace({1, 1}, counter).second);
    ASSERT_EQ(counter, 4);
    ASSERT_TRUE(tree.erase({0, 0}));
    ASSERT_EQ(counter, 3);
    ASSERT_TRUE(tree.erase({1, 0}));
    ASSERT_EQ(counter, 2);
    ASSERT_TRUE(tree.erase({0, 1}));
    ASSERT_EQ(counter, 1);
    ASSERT_TRUE(tree.erase({1, 1}));
    ASSERT_EQ(counter, 0);
    ASSERT_TRUE(tree.emplace({0, 0}, counter).second);
    ASSERT_TRUE(tree.emplace({0, 1}, counter).second);
    ASSERT_TRUE(tree.emplace({1, 0}, counter).second);
    ASSERT_TRUE(tree.emplace({1, 1}, counter).second);
    ASSERT_EQ(counter, 4);
    tree.clear();
    ASSERT_EQ(counter, 0);
}

TEST(TestSpatialTree, BinarySearch) {
    st::internal::spatial_tree<int, int, 1> tree;
    tree.emplace({0});
    tree.emplace({2});
    tree.emplace({4});
    auto nearest = tree.nearest({3});
    ASSERT_EQ(nearest.size(), 2);
    tree.erase({0});
    nearest = tree.nearest({3});
    ASSERT_EQ(nearest.size(), 2);
    tree.erase({2});
    nearest = tree.nearest({3});
    ASSERT_EQ(nearest.size(), 1);
}

TEST(TestSpatialTree, Octree) {
    st::internal::spatial_tree<int, int, 3> tree;
    tree.emplace({0, 0, 0});
    tree.emplace({2, 2, 2});
    tree.emplace({4, 4, 4});
    auto nearest = tree.nearest({3, 3, 3});
    ASSERT_EQ(nearest.size(), 2);
    tree.erase({0, 0, 0});
    nearest = tree.nearest({3, 3, 3});
    ASSERT_EQ(nearest.size(), 2);
    tree.erase({2, 2, 2});
    nearest = tree.nearest({3, 3, 3});
    ASSERT_EQ(nearest.size(), 1);
}

TEST(TestSpatialTree, VolumeAndBSize) {
    st::internal::spatial_tree<int, int, 3> tree;

    uint64_t volume = 0;
    uint64_t bsize = 0;
    for (int i = 0; i < 10000; ++i) {
        ASSERT_GE(tree.volume(), volume);
        ASSERT_GE(tree.bsize(), bsize);
        volume = tree.volume();
        tree.emplace({i, i, -i}, i);
        ASSERT_GE(tree.volume(), volume);
        ASSERT_GE(tree.bsize(), bsize);
        volume = tree.volume();
        tree.emplace({-i, i, i}, i);
        ASSERT_GE(tree.volume(), volume);
        ASSERT_GE(tree.bsize(), bsize);
        volume = tree.volume();
        tree.emplace({i, -i, i}, i);
        ASSERT_GE(tree.volume(), volume);
        ASSERT_GE(tree.bsize(), bsize);
        volume = tree.volume();
    }

    for (int i = 0; i < 10000; ++i) {
        ASSERT_LE(tree.volume(), volume);
        ASSERT_GE(tree.bsize(), bsize);
        volume = tree.volume();
        tree.erase({i, i, -i});
        ASSERT_LE(tree.volume(), volume);
        ASSERT_GE(tree.bsize(), bsize);
        volume = tree.volume();
        tree.erase({-i, i, i});
        ASSERT_LE(tree.volume(), volume);
        ASSERT_GE(tree.bsize(), bsize);
        volume = tree.volume();
        tree.erase({i, -i, i});
        ASSERT_LE(tree.volume(), volume);
        ASSERT_GE(tree.bsize(), bsize);
        volume = tree.volume();
    }
}

TEST(TestSpatialTree, UniqueInsertions) {
    const auto typed_test = [&]<typename CoordT>() {
        const std::vector<uint64_t> test_sizes = {1, 2, 3, 7, 10, 100, 1000};
        for (uint64_t test_size : test_sizes) {
            const auto inner_test = [&]<typename TreeType>(auto &tree) {
                for (uint64_t i = 1; i <= test_size; ++i) {
                    CoordT coord = i;
                    ASSERT_TRUE(tree.emplace({coord, coord}, i).second);
                    ASSERT_TRUE(tree.emplace({coord, -coord}, i).second);
                    ASSERT_TRUE(tree.emplace({-coord, coord}, i).second);
                    ASSERT_TRUE(tree.emplace({-coord, -coord}, i).second);
                }
                ASSERT_EQ(tree.size(), 4 * test_size);
                ASSERT_EQ(tree.size(), std::distance(tree.begin(), tree.end()));
                tree.clear();
            };
            CoordT                            bound = test_size;
            const st::bounding_box<CoordT, 2> bounds = {-bound, -bound, bound, bound};

            st::internal::spatial_tree<CoordT, int, 2> bounded(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2>>(bounded);

            st::internal::spatial_tree<CoordT, int, 2> unbounded;
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2>>(unbounded);

            st::internal::spatial_tree<CoordT, int, 2, 1> small_size(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 1>>(
                small_size);

            st::internal::spatial_tree<CoordT, int, 2, 16> large_size(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 16>>(
                large_size);

            st::internal::spatial_tree<CoordT, int, 2, 7> odd_size(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 7>>(odd_size);
        }
    };
    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestSpatialTree, UniqueDeletions) {
    const auto typed_test = [&]<typename CoordT>() {
        const std::vector<uint64_t> test_sizes = {1, 2, 3, 7, 10, 100, 1000};
        for (uint64_t test_size : test_sizes) {
            const auto inner_test = [&]<typename TreeType>(auto &tree) {
                for (uint64_t i = 1; i <= test_size; ++i) {
                    CoordT coord = i;
                    ASSERT_TRUE(tree.emplace({coord, coord}, i).second);
                    ASSERT_TRUE(tree.emplace({coord, -coord}, i).second);
                    ASSERT_TRUE(tree.emplace({-coord, coord}, i).second);
                    ASSERT_TRUE(tree.emplace({-coord, -coord}, i).second);
                }
                ASSERT_TRUE(!tree.empty());
                for (uint64_t i = 1; i <= test_size; ++i) {
                    CoordT coord = i;
                    ASSERT_TRUE(tree.erase({coord, coord}));
                    ASSERT_TRUE(tree.erase({coord, -coord}));
                    ASSERT_TRUE(tree.erase({-coord, coord}));
                    ASSERT_TRUE(tree.erase({-coord, -coord}));
                }
                ASSERT_TRUE(tree.empty());
                ASSERT_EQ(tree.size(), std::distance(tree.begin(), tree.end()));
                for (uint64_t i = 1; i <= test_size; ++i) {
                    CoordT coord = i;
                    ASSERT_TRUE(tree.emplace({coord, coord}, i).second);
                    ASSERT_TRUE(tree.emplace({coord, -coord}, i).second);
                    ASSERT_TRUE(tree.emplace({-coord, coord}, i).second);
                    ASSERT_TRUE(tree.emplace({-coord, -coord}, i).second);
                }
                ASSERT_TRUE(!tree.empty());
                for (uint64_t i = 1; i <= test_size; ++i) {
                    CoordT coord = i;
                    ASSERT_TRUE(tree.erase({coord, coord}));
                    ASSERT_TRUE(tree.erase({coord, -coord}));
                    ASSERT_TRUE(tree.erase({-coord, coord}));
                    ASSERT_TRUE(tree.erase({-coord, -coord}));
                }
                ASSERT_TRUE(tree.empty());
                ASSERT_EQ(tree.size(), std::distance(tree.begin(), tree.end()));
            };
            CoordT                            bound = test_size;
            const st::bounding_box<CoordT, 2> bounds = {-bound, -bound, bound, bound};

            st::internal::spatial_tree<CoordT, int, 2> bounded(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2>>(bounded);

            st::internal::spatial_tree<CoordT, int, 2> unbounded;
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2>>(unbounded);

            st::internal::spatial_tree<CoordT, int, 2, 1> small_size(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 1>>(
                small_size);

            st::internal::spatial_tree<CoordT, int, 2, 16> large_size(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 16>>(
                large_size);

            st::internal::spatial_tree<CoordT, int, 2, 7> odd_size(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 7>>(odd_size);
        }
    };
    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestSpatialTree, RandomInsertionsAndDeletions) {
    static constexpr int DISTRIBUTION_BEG = -10000;
    static constexpr int DISTRIBUTION_END = 10000;
    static constexpr int RANGE_SIZE = DISTRIBUTION_END - DISTRIBUTION_BEG + 1;
    static_assert(RANGE_SIZE * RANGE_SIZE <= std::numeric_limits<int>::max());
    std::random_device            device;
    std::uniform_int_distribution distribution(static_cast<int>(DISTRIBUTION_BEG),
                                               static_cast<int>(DISTRIBUTION_END));

    const std::vector<uint64_t> test_sizes = {1, 2, 3, 7, 10, 100, 1000};
    const auto                  typed_test = [&]<typename CoordT>() {
        const auto inner_test = [&]<typename TreeType>(auto &tree) {
            for (uint64_t test_size : test_sizes) {
                std::unordered_set<int> added_points;
                for (uint64_t i = 0; i < test_size; ++i) {
                    const int  x = distribution(device);
                    const int  y = distribution(device);
                    const bool was_inserted_in_tree =
                        tree.emplace({CoordT(x), CoordT(y)}, i).second;
                    const int64_t unique_hash = x + RANGE_SIZE * y;
                    const bool    was_inserted_in_map = added_points.insert(unique_hash).second;
                    ASSERT_EQ(was_inserted_in_map, was_inserted_in_tree);
                }

                ASSERT_EQ(tree.size(), std::distance(tree.begin(), tree.end()));
                ASSERT_EQ(added_points.size(), tree.size());

                for (uint64_t i = 0; i < 100 * test_size; ++i) {
                    const int     x = distribution(device);
                    const int     y = distribution(device);
                    const bool    was_deleted_from_tree = tree.erase({CoordT(x), CoordT(y)});
                    const int64_t unique_hash = x + RANGE_SIZE * y;
                    const bool    was_deleted_from_map = added_points.erase(unique_hash);
                    ASSERT_EQ(was_deleted_from_map, was_deleted_from_tree);
                }
                ASSERT_EQ(added_points.size(), tree.size());
                tree.clear();
            }
        };
        const st::bounding_box<CoordT, 2> bounds = {
            DISTRIBUTION_BEG, DISTRIBUTION_BEG, DISTRIBUTION_END, DISTRIBUTION_END};

        st::internal::spatial_tree<CoordT, int, 2> bounded(bounds);
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2>>(bounded);

        st::internal::spatial_tree<CoordT, int, 2, 1> small_size(bounds);
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 1>>(small_size);

        st::internal::spatial_tree<CoordT, int, 2, 16> large_size(bounds);
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 16>>(large_size);

        st::internal::spatial_tree<CoordT, int, 2, 7> odd_size(bounds);
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 7>>(odd_size);
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestSpatialTree, IteratorCoverage) {
    const std::vector<uint64_t> test_sizes = {1, 2, 3, 7, 10, 100, 1000};
    const auto                  typed_test = [&]<typename CoordT>() {
        for (int test_size : test_sizes) {
            const auto inner_test = [&]<typename TreeType>(auto &tree) {
                ASSERT_FALSE(std::distance(tree.begin(), tree.end()));

                for (int i = 0; i < test_size; ++i) {
                    ASSERT_TRUE(tree.emplace({CoordT(i), CoordT(i)}, i).second);
                }
                ASSERT_EQ(tree.size(), test_size);
                ASSERT_EQ(std::distance(tree.begin(), tree.end()), test_size);
            };
            CoordT                            boundary = test_size;
            const st::bounding_box<CoordT, 2> bounds = {-boundary, -boundary, boundary, boundary};

            st::internal::spatial_tree<CoordT, int, 2> bounded(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2>>(bounded);

            st::internal::spatial_tree<CoordT, int, 2> unbounded;
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2>>(unbounded);

            st::internal::spatial_tree<CoordT, int, 2, 1> small_size(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 1>>(
                small_size);

            st::internal::spatial_tree<CoordT, int, 2, 16> large_size(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 16>>(
                large_size);

            st::internal::spatial_tree<CoordT, int, 2, 7> odd_size(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 7>>(odd_size);
        }
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestSpatialTree, SingleFind) {
    const std::vector<uint64_t> test_sizes = {1, 2, 3, 7, 10, 100, 1000};
    const auto                  typed_test = [&]<typename CoordT>() {
        for (int test_size : test_sizes) {
            const auto inner_test = [&]<typename TreeType>(auto &tree) {
                ASSERT_EQ(tree.find({0, 0}), tree.end());

                for (int i = 0; i < test_size; ++i) {
                    ASSERT_TRUE(tree.emplace({CoordT(i), CoordT(i)}, i).second);
                }
                ASSERT_EQ(tree.size(), test_size);

                for (int i = 0; i < test_size; ++i) {
                    ASSERT_NE(tree.find({CoordT(i), CoordT(i)}), tree.end());
                }
            };
            CoordT                            boundary = test_size;
            const st::bounding_box<CoordT, 2> bounds = {-boundary, -boundary, boundary, boundary};

            st::internal::spatial_tree<CoordT, int, 2> bounded(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2>>(bounded);

            st::internal::spatial_tree<CoordT, int, 2> unbounded;
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2>>(unbounded);

            st::internal::spatial_tree<CoordT, int, 2, 1> small_size(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 1>>(
                small_size);

            st::internal::spatial_tree<CoordT, int, 2, 16> large_size(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 16>>(
                large_size);

            st::internal::spatial_tree<CoordT, int, 2, 7> odd_size(bounds);
            inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 7>>(odd_size);
        }
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestSpatialTree, BBoxFind) {
    static constexpr int DISTRIBUTION_BEG = -10000;
    static constexpr int DISTRIBUTION_END = 10000;
    static constexpr int RANGE_SIZE = DISTRIBUTION_END - DISTRIBUTION_BEG + 1;
    static_assert(RANGE_SIZE * RANGE_SIZE <= std::numeric_limits<int>::max());
    std::random_device            device;
    std::uniform_int_distribution distribution(static_cast<int>(DISTRIBUTION_BEG),
                                               static_cast<int>(DISTRIBUTION_END));

    const std::vector<uint64_t> test_sizes = {1, 2, 3, 7, 10, 100, 1000};
    const auto                  typed_test = [&]<typename CoordT>() {
        const auto inner_test = [&]<typename TreeType>(auto &tree) {
            for (uint64_t test_size : test_sizes) {
                std::unordered_set<int> added_points;
                for (uint64_t i = 0; i < test_size; ++i) {
                    const int x = distribution(device);
                    const int y = distribution(device);
                    tree.emplace({CoordT(x), CoordT(y)}, i);
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

                    const st::bounding_box<CoordT, 2> bbox = {x_start, y_start, x_stop, y_stop};

                    int counter = std::count_if(tree.begin(), tree.end(), [&](auto it) {
                        const auto [coordinates, val] = it;
                        return bbox.contains(coordinates);
                    });

                    int rcounter = 0;
                    tree.find(bbox, [&](auto) { ++rcounter; });

                    ASSERT_EQ(counter, rcounter);
                }
            }
        };
        const st::bounding_box<CoordT, 2> bounds = {
            DISTRIBUTION_BEG, DISTRIBUTION_BEG, DISTRIBUTION_END, DISTRIBUTION_END};

        st::internal::spatial_tree<CoordT, int, 2> bounded(bounds);
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2>>(bounded);

        st::internal::spatial_tree<CoordT, int, 2> unbounded;
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2>>(unbounded);

        st::internal::spatial_tree<CoordT, int, 2, 1> small_size(bounds);
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 1>>(small_size);

        st::internal::spatial_tree<CoordT, int, 2, 16> large_size(bounds);
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 16>>(large_size);

        st::internal::spatial_tree<CoordT, int, 2, 7> odd_size(bounds);
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 7>>(odd_size);
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}

TEST(TestSpatialTree, Nearest) {
    static constexpr int DISTRIBUTION_BEG = -10000;
    static constexpr int DISTRIBUTION_END = 10000;
    static constexpr int RANGE_SIZE = DISTRIBUTION_END - DISTRIBUTION_BEG + 1;
    static_assert(RANGE_SIZE * RANGE_SIZE <= std::numeric_limits<int>::max());
    std::random_device            device;
    std::uniform_int_distribution distribution(static_cast<int>(DISTRIBUTION_BEG),
                                               static_cast<int>(DISTRIBUTION_END));

    const std::vector<uint64_t> test_sizes = {1, 2, 3, 7, 10, 100, 1000};
    const auto                  typed_test = [&]<typename CoordT>() {
        const auto inner_test = [&]<typename TreeType>(auto &tree) {
            for (uint64_t test_size : test_sizes) {
                std::unordered_set<int> added_points;
                for (uint64_t i = 0; i < test_size; ++i) {
                    const int x = distribution(device);
                    const int y = distribution(device);
                    tree.emplace({CoordT(x), CoordT(y)}, i);

                    const CoordT x_ = distribution(device);
                    const CoordT y_ = distribution(device);

                    const auto [pt, nearest_val1] =
                        *std::min_element(tree.begin(), tree.end(), [&](auto lhs, auto rhs) {
                            const auto [p1, val1] = lhs;
                            const auto [p2, val2] = rhs;

                            return st::internal::euclidean_distance_squared(p1[0], x_, p1[1], y_) <
                                   st::internal::euclidean_distance_squared(p2[0], x_, p2[1], y_);
                        });
                    const CoordT distance1 =
                        st::internal::euclidean_distance_squared(x_, pt[0], y_, pt[1]);

                    const auto nearest_points = tree.nearest({CoordT(x_), CoordT(y_)});
                    ASSERT_FALSE(nearest_points.empty());
                    ASSERT_TRUE(
                        std::all_of(nearest_points.begin(), nearest_points.end(), [&](auto it) {
                            const auto [p2, nearest_val2] = it;
                            const CoordT distance2 =
                                st::internal::euclidean_distance_squared(x_, p2[0], y_, p2[1]);
                            return distance2 == distance1;
                        }));
                }
            }
        };
        const st::bounding_box<CoordT, 2> bounds = {
            DISTRIBUTION_BEG, DISTRIBUTION_BEG, DISTRIBUTION_END, DISTRIBUTION_END};

        st::internal::spatial_tree<CoordT, int, 2> bounded(bounds);
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2>>(bounded);

        st::internal::spatial_tree<CoordT, int, 2> unbounded;
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2>>(unbounded);

        st::internal::spatial_tree<CoordT, int, 2, 1> small_size(bounds);
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 1>>(small_size);

        st::internal::spatial_tree<CoordT, int, 2, 16> large_size(bounds);
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 16>>(large_size);

        st::internal::spatial_tree<CoordT, int, 2, 7> odd_size(bounds);
        inner_test.template operator()<st::internal::spatial_tree<CoordT, int, 2, 7>>(odd_size);
    };

    typed_test.operator()<float>();
    typed_test.operator()<double>();
    typed_test.operator()<int>();
}
