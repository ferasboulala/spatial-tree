#include <benchmark/benchmark.h>

#include <limits>
#include <optional>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

#include "spatial.h"

#define BEG_TYPELESS -10000
#define END_TYPELESS 10000
static_assert(END_TYPELESS > BEG_TYPELESS);

static constexpr uint64_t LOW = 1 << 10;
static constexpr uint64_t HIGH = 1 << 20;

struct opaque_data {
    char opaque[1];
    opaque_data() {}
};
using VALUE_TYPE = opaque_data;
static constexpr uint64_t MAX_NODE_SIZE = 32;

template <typename CoordT>
std::vector<std::pair<CoordT, CoordT>> generate_points(CoordT BEG, CoordT END, uint64_t test_size) {
    using DistributionType = std::conditional<std::is_floating_point_v<CoordT>,
                                              std::uniform_real_distribution<CoordT>,
                                              std::uniform_int_distribution<CoordT>>::type;
    DistributionType                       distribution(BEG, END);
    std::default_random_engine             device;
    std::vector<std::pair<CoordT, CoordT>> points;
    for (uint64_t i = 0; i < test_size; ++i) {
        CoordT x = distribution(device);
        CoordT y = distribution(device);
        points.emplace_back(x, y);
    }

    return points;
}

template <typename CoordT, typename TreeType>
void benchmark_insertions(benchmark::State &state, TreeType &tree) {
    static_assert(std::is_signed_v<CoordT>);

    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(BEG, END, test_size);

    uint64_t n_points_inserted = 0;
    uint64_t n_collisions = 0;
    for (auto _ : state) {
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            bool inserted = tree.emplace(x, y).second;
            n_points_inserted += inserted;
            n_collisions += (1 - inserted);
        }
        tree.clear();
    }
    state.counters["inserted"] = n_points_inserted / state.iterations();
    state.counters["collisions"] = n_collisions / state.iterations();
    state.SetItemsProcessed(test_size * state.iterations());
}

template <typename CoordT, typename TreeType>
void benchmark_find(benchmark::State &state, TreeType &tree) {
    static_assert(std::is_signed_v<CoordT>);

    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;
    static constexpr CoordT AREA = (END - BEG) * (END - BEG);
    static_assert(AREA > BEG && AREA > END);

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(BEG, END, test_size);

    for (auto _ : state) {
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            tree.emplace(x, y);
        }
    }

    static constexpr uint64_t n_bounding_boxes = 1 << 7;
    const auto bounding_boxes_points = generate_points(BEG, END, n_bounding_boxes * 2);
    std::vector<st::bounding_box<CoordT>> bounding_boxes;
    bounding_boxes.reserve(n_bounding_boxes);
    CoordT average_area = 0;
    for (uint64_t i = 0; i < 2 * n_bounding_boxes; i += 2) {
        auto [x1, x2] = bounding_boxes_points[i];
        auto [y1, y2] = bounding_boxes_points[i + 1];
        auto [top_x, bottom_x] = std::minmax(x1, x2);
        auto [bottom_y, top_y] = std::minmax(y1, y2);
        bounding_boxes.emplace_back(top_x, top_y, bottom_x, bottom_y);
        average_area += bounding_boxes.back().area();
    }
    average_area /= n_bounding_boxes;

    uint64_t n_points_found = 0;
    for (auto _ : state) {
        for (const auto &bbox : bounding_boxes) {
            tree.find(bbox, [&](auto) { ++n_points_found; });
        }
    }
    state.counters["avg_area (%)"] = average_area / AREA;
    state.counters["found"] = n_points_found / state.iterations() / n_bounding_boxes;
    state.SetItemsProcessed(n_bounding_boxes * state.iterations());
}

template <typename CoordT, typename TreeType>
void benchmark_find_single(benchmark::State &state, TreeType &tree) {
    static_assert(std::is_signed_v<CoordT>);

    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;
    static constexpr CoordT AREA = (END - BEG) * (END - BEG);
    static_assert(AREA > BEG && AREA > END);

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(BEG, END, test_size);

    for (auto _ : state) {
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            tree.emplace(x, y);
        }
    }

    static constexpr uint64_t n_points_to_find = 1 << 7;
    const auto                points_to_find = generate_points(BEG, END, n_points_to_find);

    for (auto _ : state) {
        for (auto [x, y] : points_to_find) {
            tree.find(x, y);
        }
    }
    state.SetItemsProcessed(n_points_to_find * state.iterations());
}

template <typename CoordT, typename TreeType>
void benchmark_nearest(benchmark::State &state, TreeType &tree) {
    static_assert(std::is_signed_v<CoordT>);

    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(BEG, END, test_size);

    for (auto _ : state) {
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            tree.emplace(x, y);
        }
    }

    static constexpr uint64_t n_queries = 1 << 10;
    const auto                query_points = generate_points(BEG, END, n_queries);
    uint64_t                  n_points_found = 0;
    for (auto _ : state) {
        for (auto [x, y] : query_points) {
            n_points_found += tree.nearest(x, y).size();
        }
    }
    state.counters["found"] = n_points_found / state.iterations();
    state.SetItemsProcessed(n_queries * state.iterations());
}

template <typename StorageType, typename CoordinateType>
class hash_table_oracle {
private:
    struct HashFunc {
        __always_inline size_t
        operator()(const std::pair<CoordinateType, CoordinateType> &coords) const {
            auto [x, y] = coords;

            size_t seed = 0;
            seed ^= std::hash<double>()(x) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= std::hash<double>()(y) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

            return seed;
        }
    };

    using TableType =
        std::unordered_map<std::pair<CoordinateType, CoordinateType>, StorageType, HashFunc>;

public:
    __always_inline void clear() { table_.clear(); }

    template <typename... Args>
    __always_inline auto emplace(CoordinateType x, CoordinateType y, Args &&...args) {
        return table_.emplace(std::pair<CoordinateType, CoordinateType>{x, y},
                              StorageType(std::forward<Args>(args)...));
    }

    void find(
        const st::bounding_box<CoordinateType> &bbox,
        std::function<void(std::pair<const std::pair<CoordinateType, CoordinateType>, StorageType>)>
            func) const {
        for (const auto &it : table_) {
            const auto [x, y] = it.first;
            if (st::internal::is_inside_bounding_box(x, y, bbox)) {
                func(it);
            }
        }
    }

    __always_inline auto find(CoordinateType x, CoordinateType y) const {
        return table_.find({x, y});
    }

    auto nearest(CoordinateType x, CoordinateType y) const {
        std::vector<std::pair<const std::pair<CoordinateType, CoordinateType>, StorageType>>
             nearest_points;
        auto shortest_distance_squared = std::numeric_limits<CoordinateType>::max();
        for (const auto &it : table_) {
            const auto &[x_, y_] = it.first;
            const auto dist = st::internal::euclidean_distance_squared(x, x_, y, y_);
            if (dist < shortest_distance_squared) {
                shortest_distance_squared = dist;
                nearest_points.clear();
                nearest_points.push_back(it);
            } else if (dist == shortest_distance_squared) {
                nearest_points.emplace_back(it);
            }
        }

        return nearest_points;
    }

    static void insertions(benchmark::State &state) {
        auto tree = hash_table_oracle<VALUE_TYPE, CoordinateType>();
        benchmark_insertions<CoordinateType>(state, tree);
    }

    static void find(benchmark::State &state) {
        auto tree = hash_table_oracle<VALUE_TYPE, CoordinateType>();
        benchmark_find<CoordinateType>(state, tree);
    }

    static void find_single(benchmark::State &state) {
        auto tree = hash_table_oracle<VALUE_TYPE, CoordinateType>();
        benchmark_find_single<CoordinateType>(state, tree);
    }

    static void nearest(benchmark::State &state) {
        auto tree = hash_table_oracle<VALUE_TYPE, CoordinateType>();
        benchmark_nearest<CoordinateType>(state, tree);
    }

private:
    TableType table_;
};

template <typename CoordT>
void insertions(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::internal::spatial_tree<VALUE_TYPE, CoordT, MAX_NODE_SIZE>({BEG, END, END, BEG});
    benchmark_insertions<CoordT>(state, tree);
}

template <typename CoordT>
void find(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::internal::spatial_tree<VALUE_TYPE, CoordT, MAX_NODE_SIZE>({BEG, END, END, BEG});
    benchmark_find<CoordT>(state, tree);
}

template <typename CoordT>
void find_single(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::internal::spatial_tree<VALUE_TYPE, CoordT, MAX_NODE_SIZE>({BEG, END, END, BEG});
    benchmark_find_single<CoordT>(state, tree);
}

template <typename CoordT>
void nearest(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::internal::spatial_tree<VALUE_TYPE, CoordT, MAX_NODE_SIZE>({BEG, END, END, BEG});
    benchmark_nearest<CoordT>(state, tree);
}

using COORD_TYPE = float;

BENCHMARK(insertions<COORD_TYPE>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(find<COORD_TYPE>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(find_single<COORD_TYPE>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(nearest<COORD_TYPE>)->RangeMultiplier(2)->Range(LOW, HIGH);

BENCHMARK(hash_table_oracle<int, COORD_TYPE>::insertions)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(hash_table_oracle<int, COORD_TYPE>::find)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(hash_table_oracle<int, COORD_TYPE>::find_single)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(hash_table_oracle<int, COORD_TYPE>::nearest)->RangeMultiplier(2)->Range(LOW, HIGH);

BENCHMARK_MAIN();