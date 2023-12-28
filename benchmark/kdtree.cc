#include "spatial-tree/kdtree.hh"

#include <benchmark/benchmark.h>

#include <vector>

#include "benchmark/common.hh"

using VALUE_TYPE = char;

template <typename CoordT>
void insertions(benchmark::State &state) {
    auto tree = st::KDTree<VALUE_TYPE, CoordT>();
    benchmark_insertions<CoordT>(state, tree);
}

template <typename CoordT>
void find(benchmark::State &state) {
    auto tree = st::KDTree<VALUE_TYPE, CoordT>();
    benchmark_find<CoordT>(state, tree);
}

template <typename CoordT>
void find_single(benchmark::State &state) {
    auto tree = st::KDTree<VALUE_TYPE, CoordT>();
    benchmark_find_single<CoordT>(state, tree);
}

template <typename CoordT>
void nearest(benchmark::State &state) {
    auto tree = st::KDTree<VALUE_TYPE, CoordT>();
    benchmark_nearest<CoordT>(state, tree);
}

static constexpr uint64_t LOW = 1 << 10;
static constexpr uint64_t HIGH = 1 << 20;

BENCHMARK(insertions<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(find<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(find_single<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK(nearest<double>)->RangeMultiplier(2)->Range(LOW, HIGH);
BENCHMARK_MAIN();