#include <benchmark/benchmark.h>

#include <limits>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

#include "benchmark/common.hh"
#include "spatial-tree/euclide.hh"
#include "spatial-tree/intervals.hh"

using VALUE_TYPE = char;

template <typename StorageType, typename CoordinateType>
class HashTable {
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
        return table_.emplace(std::pair<CoordinateType, CoordinateType>{x, y}, StorageType(std::forward<Args>(args)...));
    }

    void find(
        const st::BoundingBox<CoordinateType> &bbox,
        std::function<void(std::pair<const std::pair<CoordinateType, CoordinateType>, StorageType>)>
            func) const {
        for (const auto &it : table_) {
            const auto [x, y] = it.first;
            if (st::is_inside_bounding_box(x, y, bbox)) {
                func(it);
            }
        }
    }

    __always_inline auto find(CoordinateType x, CoordinateType y) const { return table_.find({x, y}); }

    auto nearest(CoordinateType x, CoordinateType y) const {
        std::vector<std::pair<const std::pair<CoordinateType, CoordinateType>, StorageType>>
             nearest_points;
        auto shortest_distance_squared = std::numeric_limits<CoordinateType>::max();
        for (const auto &it : table_) {
            const auto &[x_, y_] = it.first;
            const auto dist = st::euclidean_distance_squared(x, x_, y, y_);
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

private:
    TableType table_;
};

template <typename CoordT>
void insertions(benchmark::State &state) {
    auto tree = HashTable<VALUE_TYPE, CoordT>();
    benchmark_insertions<CoordT>(state, tree);
}

template <typename CoordT>
void find(benchmark::State &state) {
    auto tree = HashTable<VALUE_TYPE, CoordT>();
    benchmark_find<CoordT>(state, tree);
}

template <typename CoordT>
void find_single(benchmark::State &state) {
    auto tree = HashTable<VALUE_TYPE, CoordT>();
    benchmark_find_single<CoordT>(state, tree);
}

template <typename CoordT>
void nearest(benchmark::State &state) {
    auto tree = HashTable<VALUE_TYPE, CoordT>();
    benchmark_nearest<CoordT>(state, tree);
}

static constexpr uint64_t LOW = 1 << 10;
static constexpr uint64_t HIGH = 1 << 20;

BENCHMARK(insertions<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(find<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(find_single<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(nearest<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK_MAIN();