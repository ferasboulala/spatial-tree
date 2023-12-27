#include "spatial-tree/quadtree.hh"

#include <benchmark/benchmark.h>

#include <vector>

#include "benchmark/common.hh"

using VALUE_TYPE = char;
static constexpr uint64_t MAX_NODE_SIZE = 4;

template <typename CoordT>
void insertions(benchmark::State &state) {
    static_assert(std::is_signed_v<CoordT>);

    static constexpr CoordT BEG = -10000;
    static constexpr CoordT END = 10000;

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(BEG, END, test_size);

    st::QuadTree<VALUE_TYPE, CoordT, MAX_NODE_SIZE> tree({BEG, END, END, BEG});
    for (auto _ : state) {
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            tree.emplace(x, y, i);
        }
    }
    // state.SetItemsProcessed(test_size * state.iterations());
}

template <typename CoordT>
void find(benchmark::State &state) {
    static_assert(std::is_signed_v<CoordT>);

    static constexpr CoordT BEG = -10000;
    static constexpr CoordT END = 10000;

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(BEG, END, test_size);

    st::QuadTree<VALUE_TYPE, CoordT, MAX_NODE_SIZE> tree({BEG, END, END, BEG});
    for (auto _ : state) {
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            tree.emplace(x, y, i);
        }
    }

    const uint64_t n_bounding_boxes =
        std::min<uint64_t>(100, std::max<uint64_t>(100, test_size / 1024));
    const auto bounding_boxes_points = generate_points(BEG, END, n_bounding_boxes * 2);
    std::vector<st::BoundingBox<CoordT>> bounding_boxes;
    for (uint64_t i = 0; i < 2 * n_bounding_boxes; i += 2) {
        auto [x1, x2] = bounding_boxes_points[i];
        auto [y1, y2] = bounding_boxes_points[i + 1];
        auto [top_x, bottom_x] = std::minmax(x1, x2);
        auto [bottom_y, top_y] = std::minmax(y1, y2);
        bounding_boxes.emplace_back(top_x, top_y, bottom_x, bottom_y);
    }

    uint64_t n_points_found = 0;
    for (auto _ : state) {
        for (const auto &bbox : bounding_boxes) {
            tree.find(bbox, [&](auto) { ++n_points_found; });
        }
    }
    // state.SetItemsProcessed(n_points_found);
}

template <typename CoordT>
void nearest(benchmark::State &state) {
    static_assert(std::is_signed_v<CoordT>);

    static constexpr CoordT BEG = -10000;
    static constexpr CoordT END = 10000;

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(BEG, END, test_size);

    st::QuadTree<VALUE_TYPE, CoordT, MAX_NODE_SIZE> tree({BEG, END, END, BEG});
    for (auto _ : state) {
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            tree.emplace(x, y, i);
        }
    }

    const uint64_t n_queries = std::min<uint64_t>(100, std::max<uint64_t>(100, test_size / 1024));
    const auto     query_points = generate_points(BEG, END, n_queries);
    uint64_t       n_points_found = 0;
    for (auto _ : state) {
        for (auto [x, y] : query_points) {
            n_points_found += tree.nearest(x, y).size();
        }
    }
    // state.SetItemsProcessed(n_points_found);
}

static constexpr uint64_t LOW = 1 << 10;
static constexpr uint64_t HIGH = 1 << 20;

BENCHMARK(insertions<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(find<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(nearest<double>)->RangeMultiplier(2)->Range(LOW, HIGH);

// BENCHMARK(insertions<int>)->RangeMultiplier(2)->Range(LOW, HIGH);
// BENCHMARK(find<int>)->RangeMultiplier(2)->Range(LOW, HIGH);
// BENCHMARK(nearest<int>)->RangeMultiplier(2)->Range(LOW, HIGH);

BENCHMARK_MAIN();