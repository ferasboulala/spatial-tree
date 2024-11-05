#define NDEBUG

#include <benchmark/benchmark.h>

#include <limits>
#include <random>
#include <utility>
#include <vector>

#include "spatial-tree.h"

#define BEG_TYPELESS -10000
#define END_TYPELESS 10000
static_assert(END_TYPELESS > BEG_TYPELESS);

static constexpr uint64_t LOW = 1 << 10;
static constexpr uint64_t HIGH = 1 << 22;

struct opaque_data {
    char opaque[sizeof(void *)];
    opaque_data() {}
};
static constexpr uint64_t MAX_NODE_SIZE = 256;

template <typename CoordT>
std::vector<std::pair<CoordT, CoordT>> generate_points(CoordT BEG, CoordT END, uint64_t test_size) {
    using DistributionType = typename std::conditional<std::is_floating_point_v<CoordT>,
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

    tree.reserve(test_size);
    for (auto _ : state) {
        state.PauseTiming();
        tree.clear();
        state.ResumeTiming();
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            tree.emplace({x, y}).second;
        }
        benchmark::ClobberMemory();
    }
    state.counters["volume"] = tree.volume();
    state.SetItemsProcessed(test_size * state.iterations());
}

template <typename CoordT, typename TreeType>
void benchmark_insertions_duplicate(benchmark::State &state, TreeType &tree) {
    static_assert(std::is_signed_v<CoordT>);

    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(BEG, END, test_size);
    tree.clear();
    tree.reserve(test_size);
    for (uint64_t i = 0; i < test_size; ++i) {
        auto [x, y] = points[i];
        tree.emplace({x, y});
    }

    for (auto _ : state) {
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            tree.emplace({x, y});
        }
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(test_size * state.iterations());
}

template <typename CoordT, typename TreeType>
void benchmark_deletions(benchmark::State &state, TreeType &tree) {
    static_assert(std::is_signed_v<CoordT>);

    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(BEG, END, test_size);

    tree.reserve(test_size);
    tree.clear();
    for (uint64_t i = 0; i < test_size; ++i) {
        auto [x, y] = points[i];
        tree.emplace({x, y});
    }

    for (auto _ : state) {
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            tree.erase({x, y});
        }
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(test_size * state.iterations());
}

template <typename CoordT, typename TreeType>
void benchmark_deletions_non_existent(benchmark::State &state, TreeType &tree) {
    static_assert(std::is_signed_v<CoordT>);

    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(BEG, END, test_size);

    tree.reserve(test_size);
    tree.clear();
    for (uint64_t i = 0; i < test_size; ++i) {
        auto [x, y] = points[i];
        tree.emplace({x, y});
    }

    const auto other_points = generate_points(BEG, END, test_size);
    for (auto _ : state) {
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = other_points[i];
            tree.erase({x, y});
        }
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(test_size * state.iterations());
}

template <typename CoordT, typename TreeType>
void benchmark_find(benchmark::State &state, TreeType &tree) {
    static_assert(std::is_signed_v<CoordT>);

    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;
    static constexpr CoordT RANGE = END - BEG;
    static constexpr double AREA = RANGE * RANGE;
    static constexpr double BBOX_DIM_RATIO = 0.1;
    static constexpr double BBOX_DIM_SIZE = BBOX_DIM_RATIO * RANGE;
    static_assert(AREA > BEG && AREA > END);
    static_assert(BBOX_DIM_SIZE);

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(BEG, END, test_size);

    tree.reserve(test_size);
    tree.clear();
    for (uint64_t i = 0; i < test_size; ++i) {
        auto [x, y] = points[i];
        tree.emplace({x, y});
    }

    static constexpr uint64_t n_bounding_boxes = 1 << 16;
    const auto                bounding_boxes_points = generate_points(BEG, END, n_bounding_boxes);
    std::vector<st::bounding_box<CoordT, 2>> bounding_boxes;
    bounding_boxes.reserve(n_bounding_boxes);
    for (uint64_t i = 0; i < n_bounding_boxes; ++i) {
        auto [x1, y1] = bounding_boxes_points[i];
        CoordT x2 = x1 + BBOX_DIM_SIZE;
        CoordT y2 = y1 + BBOX_DIM_SIZE;
        bounding_boxes.push_back(st::bounding_box<CoordT, 2>({x1, y1, x2, y2}));
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

template <typename CoordT, typename TreeType>
void benchmark_find_single(benchmark::State &state, TreeType &tree) {
    static_assert(std::is_signed_v<CoordT>);

    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;
    static constexpr double AREA = (END - BEG) * (END - BEG);
    static_assert(AREA > BEG && AREA > END);

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(BEG, END, test_size);

    tree.reserve(test_size);
    tree.clear();
    for (uint64_t i = 0; i < test_size; ++i) {
        auto [x, y] = points[i];
        tree.emplace({x, y});
    }

    static constexpr uint64_t n_points_to_find = 1 << 24;
    const auto                points_to_find = generate_points(BEG, END, n_points_to_find);

    for (auto _ : state) {
        for (auto [x, y] : points_to_find) {
            benchmark::DoNotOptimize(tree.find({x, y}));
        }
        benchmark::ClobberMemory();
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

    tree.reserve(test_size);
    tree.clear();
    for (uint64_t i = 0; i < test_size; ++i) {
        auto [x, y] = points[i];
        tree.emplace({x, y});
    }

    static constexpr uint64_t n_queries = 1 << 24;
    const auto                query_points = generate_points(BEG, END, n_queries);
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

template <typename CoordT, typename TreeType>
void benchmark_iteration(benchmark::State &state, TreeType &tree) {
    static_assert(std::is_signed_v<CoordT>);

    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    const uint64_t test_size = state.range(0);
    const auto     points = generate_points(BEG, END, test_size);

    tree.reserve(test_size);
    for (auto _ : state) {
        for (uint64_t i = 0; i < test_size; ++i) {
            auto [x, y] = points[i];
            tree.emplace({x, y});
        }
    }

    for (auto _ : state) {
        int i = 0;
        for (const auto &_ : tree) {
            ++i;
        }
        assert(i == test_size);
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * test_size);
}

using VALUE_TYPE = void;

template <typename CoordT>
void insertions(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::internal::spatial_tree<CoordT, VALUE_TYPE, 2, MAX_NODE_SIZE>(
        st::bounding_box<CoordT, 2>({BEG, BEG, END, END}));
    benchmark_insertions<CoordT>(state, tree);
}

template <typename CoordT>
void insertions_duplicate(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::internal::spatial_tree<CoordT, VALUE_TYPE, 2, MAX_NODE_SIZE>(
        st::bounding_box<CoordT, 2>({BEG, BEG, END, END}));
    benchmark_insertions_duplicate<CoordT>(state, tree);
}

template <typename CoordT>
void deletions(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::internal::spatial_tree<CoordT, VALUE_TYPE, 2, MAX_NODE_SIZE>(
        st::bounding_box<CoordT, 2>({BEG, BEG, END, END}));
    benchmark_deletions<CoordT>(state, tree);
}

template <typename CoordT>
void deletions_non_existent(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::internal::spatial_tree<CoordT, VALUE_TYPE, 2, MAX_NODE_SIZE>(
        st::bounding_box<CoordT, 2>({BEG, BEG, END, END}));
    benchmark_deletions_non_existent<CoordT>(state, tree);
}

template <typename CoordT>
void find(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::internal::spatial_tree<CoordT, VALUE_TYPE, 2, MAX_NODE_SIZE>(
        st::bounding_box<CoordT, 2>({BEG, BEG, END, END}));
    benchmark_find<CoordT>(state, tree);
}

template <typename CoordT>
void find_single(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::internal::spatial_tree<CoordT, VALUE_TYPE, 2, MAX_NODE_SIZE>(
        st::bounding_box<CoordT, 2>({BEG, BEG, END, END}));
    benchmark_find_single<CoordT>(state, tree);
}

template <typename CoordT>
void nearest(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::internal::spatial_tree<CoordT, VALUE_TYPE, 2, MAX_NODE_SIZE>(
        st::bounding_box<CoordT, 2>({BEG, BEG, END, END}));
    benchmark_nearest<CoordT>(state, tree);
}

template <typename CoordT>
void iteration(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::internal::spatial_tree<CoordT, VALUE_TYPE, 2, MAX_NODE_SIZE>(
        st::bounding_box<CoordT, 2>({BEG, BEG, END, END}));
    benchmark_iteration<CoordT>(state, tree);
}

using COORD_TYPE = float;

BENCHMARK(insertions<COORD_TYPE>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(insertions_duplicate<COORD_TYPE>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(deletions<COORD_TYPE>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(deletions_non_existent<COORD_TYPE>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(find<COORD_TYPE>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(find_single<COORD_TYPE>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(nearest<COORD_TYPE>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(iteration<COORD_TYPE>)->RangeMultiplier(2)->Range(LOW, HIGH);

BENCHMARK_MAIN();
