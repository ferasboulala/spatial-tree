#define NDEBUG

#include <arm_neon.h>
#include <benchmark/benchmark.h>

#include <iostream>
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
using tree_type = st::internal::spatial_tree<coordinate_type, void, 2, 64, 32, true>;
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

static constexpr unsigned SIMDRank = 3;
static constexpr unsigned SIMDSize = 64;
static constexpr uint64_t SIMDRepititions = (1ul << 34) / SIMDSize;

inline float reduce_sum(float32x4_t vec) {
    float32x2_t sum_pair = vadd_f32(vget_low_f32(vec), vget_high_f32(vec));
    return vget_lane_f32(vpadd_f32(sum_pair, sum_pair), 0);
}

void simd(benchmark::State &state) {
    struct StructOfArrays {
        std::array<std::array<float, SIMDSize>, SIMDRank> positions;
        std::array<float, SIMDSize>                       masses;
    };
    StructOfArrays                          arrays;
    std::array<float, SIMDRank>             position;
    std::array<float, SIMDRank>             acceleration;
    std::array<float, SIMDSize>             tmp;
    std::default_random_engine              gen;
    std::uniform_int_distribution<unsigned> uniform(SIMDSize / 3, SIMDSize / 2);
    std::vector<unsigned>                   sizes(SIMDSize);
    for (unsigned i = 0; i < SIMDSize; ++i) {
        sizes[i] = uniform(gen);
    }

    memset(&arrays, 0, sizeof(arrays));
    arrays.masses.fill(1);
    tmp.fill(0);
    position.fill(1);
    acceleration.fill(0);

    int      iter = 0;
    uint64_t total = 0;
    for (auto _ : state) {
        auto size = sizes[iter % sizes.size()];
        total += size;
        ++iter;
        arrays.masses.fill(1);
        std::fill(arrays.masses.begin() + size, arrays.masses.end(), 0);
        tmp.fill(1);
        for (uint64_t i = 0; i < SIMDRank; ++i) {
            float       x = position[i];
            float32x4_t x_splat = vdupq_n_f32(x);
            _Pragma("clang loop unroll_count(2)") for (int64_t j = 0; j < size; j += 4) {
                float32x4_t p = vld1q_f32(&arrays.positions[i][j]);
                float32x4_t delta = vsubq_f32(p, x_splat);
                float32x4_t delta_squared = vmulq_f32(delta, delta);
                float32x4_t old = vld1q_f32(&tmp[j]);
                float32x4_t s = vaddq_f32(delta_squared, old);
                vst1q_f32(&tmp[j], s);
            }
        }
        _Pragma("clang loop unroll_count(2)") for (int64_t j = 0; j < size; j += 4) {
            float32x4_t t = vld1q_f32(&tmp[j]);
            float32x4_t m = vld1q_f32(&arrays.masses[j]);
            float32x4_t reciprocal = vrsqrteq_f32(t);
            float32x4_t reciprocal_squared = vmulq_f32(reciprocal, reciprocal);
            float32x4_t reciprocal_denom = vmulq_f32(reciprocal_squared, reciprocal);
            float32x4_t coeff = vmulq_f32(m, reciprocal_denom);
            vst1q_f32(&tmp[j], coeff);
        }

        for (uint64_t i = 0; i < SIMDRank; ++i) {
            float       x = position[i];
            float32x4_t x_splat = vdupq_n_f32(x);
            _Pragma("clang loop unroll_count(2)") for (int64_t j = 0; j < size; j += 4) {
                float32x4_t t = vld1q_f32(&tmp[j]);
                float32x4_t p = vld1q_f32(&arrays.positions[i][j]);
                float32x4_t delta = vsubq_f32(p, x_splat);
                float32x4_t to_reduce = vmulq_f32(delta, t);
                float       contribution = reduce_sum(to_reduce);
                acceleration[i] += contribution;
            }
        };
    }
    state.SetItemsProcessed(total);
    st::internal::unroll_for<SIMDRank>([&](auto i) { std::cout << acceleration[i] << std::endl; });
}

void print_float32x4(float32x4_t vec) {
    float values[4];         // Array to store the vector elements
    vst1q_f32(values, vec);  // Store the vector into the array

    // Print the elements
    printf("[%f, %f, %f, %f]\n", values[0], values[1], values[2], values[3]);
    exit(0);
}

void sparse(benchmark::State &state) {
    struct ArrayOfStructs {
        std::array<std::array<float, SIMDRank>, SIMDSize> positions;
        std::array<float, SIMDSize>                       masses;
    };
    ArrayOfStructs                          structs;
    std::array<float, SIMDRank>             position;
    std::array<float, SIMDRank>             acceleration;
    std::array<float, SIMDSize>             tmp;
    std::default_random_engine              gen;
    std::uniform_int_distribution<unsigned> uniform(SIMDSize / 3, SIMDSize / 2);
    std::vector<unsigned>                   sizes(SIMDSize);
    for (unsigned i = 0; i < SIMDSize; ++i) {
        sizes[i] = uniform(gen);
    }

    memset(&structs, 0, sizeof(structs));
    structs.masses.fill(1);
    tmp.fill(0);
    position.fill(1);
    acceleration.fill(0);

    int      iter = 0;
    uint64_t total = 0;
    for (auto _ : state) {
        auto size = sizes[iter % sizes.size()];
        total += size;
        ++iter;
        structs.masses.fill(1);
        std::fill(structs.masses.begin() + size, structs.masses.end(), 0);
        for (uint64_t i = 0; i < SIMDRank; ++i) {
            float       x = position[i];
            float32x4_t x_splat = vdupq_n_f32(x);
            float32x4_t one = vdupq_n_f32(1.0);
            _Pragma("clang loop unroll_count(2)") for (int64_t j = 0; j < size; j += 4) {
                float32x4_t p;
                p = vld1q_lane_f32(&structs.positions[j][i], p, 0);
                p = vld1q_lane_f32(&structs.positions[j + 1][i], p, 1);
                p = vld1q_lane_f32(&structs.positions[j + 2][i], p, 2);
                p = vld1q_lane_f32(&structs.positions[j + 3][i], p, 3);

                float32x4_t delta = vsubq_f32(p, x_splat);
                float32x4_t delta_squared_plus_one = vmlaq_f32(one, delta, delta);
                vst1q_f32(&tmp[j], delta_squared_plus_one);
            }
        }
        _Pragma("clang loop unroll_count(2)") for (int64_t j = 0; j < size; j += 4) {
            float32x4_t t = vld1q_f32(&tmp[j]);
            float32x4_t m = vld1q_f32(&structs.masses[j]);
            float32x4_t reciprocal = vrsqrteq_f32(t);
            float32x4_t reciprocal_squared = vmulq_f32(reciprocal, reciprocal);
            float32x4_t reciprocal_denom = vmulq_f32(reciprocal_squared, reciprocal);
            float32x4_t coeff = vmulq_f32(m, reciprocal_denom);
            vst1q_f32(&tmp[j], coeff);
        }

        for (uint64_t i = 0; i < SIMDRank; ++i) {
            float       x = position[i];
            float32x4_t x_splat = vdupq_n_f32(x);
            _Pragma("clang loop unroll_count(2)") for (int64_t j = 0; j < size; j += 4) {
                float32x4_t p;
                p = vld1q_lane_f32(&structs.positions[j][i], p, 0);
                p = vld1q_lane_f32(&structs.positions[j + 1][i], p, 1);
                p = vld1q_lane_f32(&structs.positions[j + 2][i], p, 2);
                p = vld1q_lane_f32(&structs.positions[j + 3][i], p, 3);

                float32x4_t t = vld1q_f32(&tmp[j]);
                float32x4_t delta = vsubq_f32(p, x_splat);
                float32x4_t to_reduce = vmulq_f32(delta, t);
                float       contribution = reduce_sum(to_reduce);
                acceleration[i] += contribution;
            }
        };
    }
    state.SetItemsProcessed(total);
    st::internal::unroll_for<SIMDRank>([&](auto i) { std::cout << acceleration[i] << std::endl; });
}

void scalar(benchmark::State &state) {
    struct ArrayOfStructs {
        std::array<std::array<float, SIMDRank>, SIMDSize> positions;
        std::array<float, SIMDSize>                       masses;
    };
    ArrayOfStructs                          structs;
    std::array<float, SIMDRank>             position;
    std::array<float, SIMDRank>             acceleration;
    std::default_random_engine              gen;
    std::uniform_int_distribution<unsigned> uniform(SIMDSize / 3, SIMDSize / 2);
    std::vector<unsigned>                   sizes(SIMDSize);
    for (unsigned i = 0; i < SIMDSize; ++i) {
        sizes[i] = uniform(gen);
    }

    memset(&structs, 0, sizeof(structs));
    structs.masses.fill(1);
    position.fill(1);
    acceleration.fill(0);

    int      iter = 0;
    uint64_t total = 0;
    for (auto _ : state) {
        auto size = sizes[iter % sizes.size()];
        total += size;
        ++iter;
        for (uint64_t i = 0; i < size; ++i) {
            auto distance_squared = st::internal::euclidean_distance_squared_arr<float, SIMDRank>(
                position, structs.positions[i]);
            distance_squared += 1;
            auto reciprocal = 1.0 / std::sqrt(distance_squared);
            auto reciprocal_squared = reciprocal * reciprocal;
            auto coeff = structs.masses[i] * reciprocal * reciprocal_squared;
            st::internal::unroll_for<SIMDRank>([&](auto j) {
                auto contribution = coeff * (structs.positions[i][j] - position[j]);
                acceleration[j] += contribution;
            });
        }
    }
    state.SetItemsProcessed(total);
    st::internal::unroll_for<SIMDRank>([&](auto i) { std::cout << acceleration[i] << std::endl; });
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
BENCHMARK(scalar)->Iterations(SIMDRepititions);
BENCHMARK(simd)->Iterations(SIMDRepititions);
BENCHMARK(sparse)->Iterations(SIMDRepititions);

BENCHMARK_MAIN();
