#pragma once

#include <random>
#include <vector>

#include "spatial-tree/helpers.hh"
#include "spatial-tree/intervals.hh"

#define BEG_TYPELESS -10000
#define END_TYPELESS 10000

static_assert(END_TYPELESS > BEG_TYPELESS);

namespace {
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
}  // namespace

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
    std::vector<st::BoundingBox<CoordT>> bounding_boxes;
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
    const auto points_to_find = generate_points(BEG, END, n_points_to_find);

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