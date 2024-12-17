#define NDEBUG

#include <benchmark/benchmark.h>

#include <random>
#include <utility>
#include <vector>

#include "spatial-tree.h"

#define BegTypeless -10000
#define EndTypeless 10000
static_assert(EndTypeless > BegTypeless);

static constexpr uint64_t BoundaryLow = 1 << 10;
static constexpr uint64_t BoundaryHigh = 1 << 22;

struct opaque_data {
    char opaque[sizeof(void *)];
    opaque_data() {}
};
using coordinate_type = float;
using tree_type = st::internal::spatial_tree<coordinate_type, void, 2, 1 << 15, 32, true>;
auto create_tree() {
    return tree_type(
        st::bounding_box<coordinate_type, 2>({BegTypeless, BegTypeless, EndTypeless, EndTypeless}));
}

std::vector<std::pair<coordinate_type, coordinate_type>> generate_points(uint64_t test_size) {
    using DistributionType =
        typename std::conditional<std::is_floating_point_v<coordinate_type>,
                                  std::uniform_real_distribution<coordinate_type>,
                                  std::uniform_int_distribution<coordinate_type>>::type;
    DistributionType                                         distribution(BegTypeless, EndTypeless);
    std::default_random_engine                               device;
    std::vector<std::pair<coordinate_type, coordinate_type>> points;
    for (uint64_t i = 0; i < test_size; ++i) {
        coordinate_type x = distribution(device);
        coordinate_type y = distribution(device);
        points.emplace_back(x, y);
    }

    return points;
}

void insertions(benchmark::State &state) {
    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(test_size);

    auto tree = create_tree();
    tree.reserve(test_size);
    for (auto _ : state) {
        state.PauseTiming();
        tree.clear();
        state.ResumeTiming();
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            benchmark::DoNotOptimize(tree.emplace({x, y}));
        }
        state.PauseTiming();
        benchmark::ClobberMemory();
        state.ResumeTiming();
    }
    state.counters["volume"] = tree.volume();
    state.SetItemsProcessed(state.iterations());
}

void insertions_duplicate(benchmark::State &state) {
    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(test_size);
    auto           tree = create_tree();
    tree.reserve(test_size);
    for (uint64_t i = 0; i < test_size; ++i) {
        auto [x, y] = points[i];
        tree.emplace({x, y});
    }

    for (auto _ : state) {
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            benchmark::DoNotOptimize(tree.emplace({x, y}));
        }
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(test_size * state.iterations());
}

void deletions(benchmark::State &state) {
    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(test_size);

    auto tree = create_tree();
    tree.reserve(test_size);
    for (auto _ : state) {
        state.PauseTiming();
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            tree.emplace({x, y});
        }
        state.ResumeTiming();
        for (int64_t i = test_size; i >= 0; --i) {
            auto [x, y] = points[i];
            tree.erase({x, y});
        }
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

void deletions_non_existent(benchmark::State &state) {
    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(test_size);

    auto tree = create_tree();
    tree.reserve(test_size);
    for (uint64_t i = 0; i < test_size; ++i) {
        auto [x, y] = points[i];
        tree.emplace({x, y});
    }

    const auto other_points = generate_points(test_size);
    for (auto _ : state) {
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = other_points[i];
            tree.erase({x, y});
        }
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(test_size * state.iterations());
}

void find(benchmark::State &state) {
    static constexpr coordinate_type BEG = BegTypeless;
    static constexpr coordinate_type END = EndTypeless;
    static constexpr coordinate_type RANGE = END - BEG;
    static constexpr double          AREA = RANGE * RANGE;
    static constexpr double          BBOX_DIM_RATIO = 0.1;
    static constexpr double          BBOX_DIM_SIZE = BBOX_DIM_RATIO * RANGE;
    static_assert(AREA > BEG && AREA > END);
    static_assert(BBOX_DIM_SIZE);

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(test_size);

    auto tree = create_tree();
    tree.reserve(test_size);
    for (uint64_t i = 0; i < test_size; ++i) {
        auto [x, y] = points[i];
        tree.emplace({x, y});
    }

    static constexpr uint64_t n_bounding_boxes = 1 << 16;
    const auto                bounding_boxes_points = generate_points(n_bounding_boxes);
    std::vector<st::bounding_box<coordinate_type, 2>> bounding_boxes;
    bounding_boxes.reserve(n_bounding_boxes);
    for (uint64_t i = 0; i < n_bounding_boxes; ++i) {
        auto [x1, y1] = bounding_boxes_points[i];
        coordinate_type x2 = x1 + BBOX_DIM_SIZE;
        coordinate_type y2 = y1 + BBOX_DIM_SIZE;
        bounding_boxes.push_back(st::bounding_box<coordinate_type, 2>({x1, y1, x2, y2}));
    }

    uint64_t n_points_found = 0;
    for (auto _ : state) {
        for (const auto &bbox : bounding_boxes) {
            tree.find(bbox, [&](const auto) { ++n_points_found; });
        }
        benchmark::ClobberMemory();
    }
    state.counters["found"] = n_points_found / state.iterations() / n_bounding_boxes;
    state.SetItemsProcessed(n_bounding_boxes * state.iterations());
}

void find_single(benchmark::State &state) {
    static constexpr coordinate_type BEG = BegTypeless;
    static constexpr coordinate_type END = EndTypeless;
    static constexpr double          AREA = (END - BEG) * (END - BEG);
    static_assert(AREA > BEG && AREA > END);

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(test_size);

    auto tree = create_tree();
    tree.reserve(test_size);
    for (uint64_t i = 0; i < test_size; ++i) {
        auto [x, y] = points[i];
        tree.emplace({x, y});
    }

    static constexpr uint64_t n_points_to_find = 1 << 24;
    const auto                points_to_find = generate_points(n_points_to_find);

    for (auto _ : state) {
        for (auto [x, y] : points_to_find) {
            benchmark::DoNotOptimize(tree.find({x, y}));
        }
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(n_points_to_find * state.iterations());
}

void contains(benchmark::State &state) {
    static constexpr coordinate_type BEG = BegTypeless;
    static constexpr coordinate_type END = EndTypeless;
    static constexpr double          AREA = (END - BEG) * (END - BEG);
    static_assert(AREA > BEG && AREA > END);

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(test_size);

    auto tree = create_tree();
    tree.reserve(test_size);
    for (uint64_t i = 0; i < test_size; ++i) {
        auto [x, y] = points[i];
        tree.emplace({x, y});
    }

    static constexpr uint64_t n_points_to_find = 1 << 24;
    const auto                points_to_find = generate_points(n_points_to_find);

    for (auto _ : state) {
        for (auto [x, y] : points_to_find) {
            benchmark::DoNotOptimize(tree.contains({x, y}));
        }
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(n_points_to_find * state.iterations());
}

void nearest(benchmark::State &state) {
    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(test_size);

    auto tree = create_tree();
    tree.reserve(test_size);
    for (uint64_t i = 0; i < test_size; ++i) {
        auto [x, y] = points[i];
        tree.emplace({x, y});
    }

    static constexpr uint64_t n_queries = 1 << 24;
    const auto                query_points = generate_points(n_queries);
    uint64_t                  n_points_found = 0;
    for (auto _ : state) {
        for (auto [x, y] : query_points) {
            n_points_found += tree.nearest({x, y}).size();
        }
        benchmark::ClobberMemory();
    }
    state.counters["found"] = n_points_found / state.iterations();
    state.SetItemsProcessed(n_queries * state.iterations());
}

void iteration(benchmark::State &state) {
    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(test_size);

    auto tree = create_tree();
    tree.reserve(test_size);
    for (uint64_t i = 0; i < test_size; ++i) {
        auto [x, y] = points[i];
        tree.emplace({x, y});
    }

    for (uint64_t i = 0; i < test_size; i += 4) {
        auto [x, y] = points[i];
        tree.erase({x, y});
    }

    for (auto _ : state) {
        int i = 0;
        for (const auto &entry : tree) {
            benchmark::DoNotOptimize(++i);
        }
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * test_size);
}

BENCHMARK(insertions)->RangeMultiplier(2)->Range(BoundaryLow, BoundaryHigh);
BENCHMARK(insertions_duplicate)->RangeMultiplier(2)->Range(BoundaryLow, BoundaryHigh);
BENCHMARK(deletions)->RangeMultiplier(2)->Range(BoundaryLow, BoundaryHigh);
BENCHMARK(deletions_non_existent)->RangeMultiplier(2)->Range(BoundaryLow, BoundaryHigh);
BENCHMARK(find)->RangeMultiplier(2)->Range(BoundaryLow, BoundaryHigh);
BENCHMARK(find_single)->RangeMultiplier(2)->Range(BoundaryLow, BoundaryHigh);
BENCHMARK(contains)->RangeMultiplier(2)->Range(BoundaryLow, BoundaryHigh);
BENCHMARK(nearest)->RangeMultiplier(2)->Range(BoundaryLow, BoundaryHigh);
BENCHMARK(iteration)->RangeMultiplier(2)->Range(BoundaryLow, BoundaryHigh);

BENCHMARK_MAIN();
