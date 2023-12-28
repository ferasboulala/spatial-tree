#include "spatial-tree/quadtree.hh"

#include <benchmark/benchmark.h>

#include <vector>

#include "benchmark/common.hh"

struct Data {
    char opaque[1];
    Data() {}
};
using VALUE_TYPE = Data;
static constexpr uint64_t MAX_NODE_SIZE = 32;

template <typename CoordT>
void insertions(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::QuadTree<VALUE_TYPE, CoordT, MAX_NODE_SIZE>({BEG, END, END, BEG});
    benchmark_insertions<CoordT>(state, tree);
}

template <typename CoordT>
void find(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::QuadTree<VALUE_TYPE, CoordT, MAX_NODE_SIZE>({BEG, END, END, BEG});
    benchmark_find<CoordT>(state, tree);
}

template <typename CoordT>
void find_single(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::QuadTree<VALUE_TYPE, CoordT, MAX_NODE_SIZE>({BEG, END, END, BEG});
    benchmark_find_single<CoordT>(state, tree);
}

template <typename CoordT>
void nearest(benchmark::State &state) {
    static constexpr CoordT BEG = BEG_TYPELESS;
    static constexpr CoordT END = END_TYPELESS;

    auto tree = st::QuadTree<VALUE_TYPE, CoordT, MAX_NODE_SIZE>({BEG, END, END, BEG});
    benchmark_nearest<CoordT>(state, tree);
}

static constexpr uint64_t LOW = 1 << 10;
static constexpr uint64_t HIGH = 1 << 20;

BENCHMARK(insertions<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(find<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(find_single<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(nearest<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK_MAIN();