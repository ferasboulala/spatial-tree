#include <raylib.h>

#include <random>

#include "spatial-tree.h"

static constexpr int Rank = 2;
static constexpr int MaximumLeafSize = 32;
using CoordinateType = float;
using MassType = float;

struct n_body_data {
    std::array<CoordinateType, Rank> position;
    std::array<CoordinateType, Rank> velocity;
    std::array<CoordinateType, Rank> acceleration;
    MassType                         mass;
};

struct n_body_tree_data_external {
    std::array<CoordinateType, Rank> position;
    MassType                         mass;
};

static inline float isqrt(float x) {
    union {
        float    f;
        uint32_t i;
    } conv;

    float x2 = x * 0.5F;
    conv.f = x;
    conv.i = 0x5f3759df - (conv.i >> 1);
    conv.f = conv.f * (1.5F - (x2 * conv.f * conv.f));

    return conv.f;
}
using spatial_tree_type =
    st::internal::spatial_tree<CoordinateType, MassType, Rank, MaximumLeafSize, 32, true>;
class barnes_hut : public spatial_tree_type {
public:
    barnes_hut() = default;
    ~barnes_hut() = default;

    void build(const auto& points) {
        st::bounding_box<CoordinateType, Rank> boundary;
        st::bounding_box<CoordinateType, Rank> boundary_global;

        boundary_global.stops.fill(std::numeric_limits<CoordinateType>::lowest() / 2);
        boundary_global.starts.fill(std::numeric_limits<CoordinateType>::max() / 2);

#pragma omp parallel private(boundary)
        {
            boundary.stops.fill(std::numeric_limits<CoordinateType>::lowest() / 2);
            boundary.starts.fill(std::numeric_limits<CoordinateType>::max() / 2);
// TODO: parallelize this
#pragma omp parallel for
            for (const auto& point : points) {
                assert(point.mass);
                st::internal::unroll_for<Rank>([&](auto i) {
                    boundary.starts[i] = std::min(point.position[i], boundary.starts[i]);
                    boundary.stops[i] = std::max(point.position[i], boundary.stops[i]);
                });
            }
#pragma omp critical
            st::internal::unroll_for<Rank>([&](auto i) {
                boundary_global.starts[i] = std::min(boundary_global.starts[i], boundary.starts[i]);
                boundary_global.stops[i] = std::max(boundary_global.stops[i], boundary.stops[i]);
            });
        }

        this->reset(boundary_global);
        for (const auto& point : points) {
            this->emplace(point.position, point.mass);
        }
    }

    void propagate() {
        branch_data_.resize(this->branches_.size());
        for (int64_t i = this->branches_.size() - 1; i >= 0; --i) {
            const auto& branch = this->branches_[i];
            auto&       branch_data_entry = branch_data_[i];
            branch_data_entry.mass = 0;
            branch_data_entry.position = {0, 0};

            if (branch.is_terminal()) {
                auto& leaf = this->leaves_[branch.index()];
                std::for_each(leaf.items.begin(),
                              leaf.items.begin() + leaf.size,
                              [&](const auto& data) { branch_data_entry.mass += data.data; });
                if (branch_data_entry.mass == 0.0) {
                    continue;
                }

                const MassType reciprocal = 1.0 / branch_data_entry.mass;
                for (uint64_t i = 0; i < leaf.size; ++i) {
                    const MassType weight = leaf.items[i].data * reciprocal;
                    st::internal::unroll_for<Rank>([&](auto j) {
                        branch_data_entry.position[j] += leaf.coordinates[i][j] * weight;
                    });
                }
            } else {
                st::internal::unroll_for<spatial_tree_type::BranchingFactor>([&](auto i) {
                    const auto& child_branch_data_entry =
                        branch_data_[branch.index_of_first_child + i];
                    branch_data_entry.mass += child_branch_data_entry.mass;
                });
                if (branch_data_entry.mass == 0.0) continue;

                const MassType reciprocal = 1.0 / branch_data_entry.mass;
                st::internal::unroll_for<spatial_tree_type::BranchingFactor>([&](auto i) {
                    const auto& child_branch_data_entry =
                        branch_data_[branch.index_of_first_child + i];
                    st::internal::unroll_for<Rank>([&](auto j) {
                        branch_data_entry.position[j] += child_branch_data_entry.position[j] *
                                                         child_branch_data_entry.mass * reciprocal;
                    });
                });
            }
        }
    }

    void update(auto& points, float dt) {
        float tree_span = 0;
        st::internal::unroll_for<Rank>([&](auto i) {
            tree_span =
                std::max<float>(tree_span, this->boundary_.stops[i] - this->boundary_.starts[i]);
        });
        float tree_span_squared = tree_span * tree_span;

#pragma omp parallel for
        for (auto& point : points) {
            st::internal::unroll_for<Rank>([&](auto i) { point.acceleration[i] = 0; });
            const std::function<void(uint64_t, float)> update_recursively = [&](auto branch_index,
                                                                                auto span_squared) {
                const auto update_data = [&](const n_body_tree_data_external& data) {
                    float distance_squared =
                        st::internal::euclidean_distance_squared_arr<CoordinateType, Rank>(
                            point.position, data.position);

                    static constexpr float G = 1.1;
                    static constexpr float EpsilonSquared = 1;
                    distance_squared += EpsilonSquared;
                    float reciprocal = isqrt(distance_squared);
                    float reciprocal_squared = reciprocal * reciprocal;
                    float coeff = G * data.mass * reciprocal * reciprocal_squared;
                    st::internal::unroll_for<Rank>([&](auto i) {
                        point.acceleration[i] += (data.position[i] - point.position[i]) * coeff;
                    });
                };

                std::array<bool, spatial_tree_type::BranchingFactor> recurse;
                const auto& branch = this->branches_[branch_index];
                st::internal::unroll_for<spatial_tree_type::BranchingFactor>([&](auto quad) {
                    const auto& branch_ = this->branches_[branch.index_of_first_child + quad];
                    const auto& branch_data__ = branch_data_[branch.index_of_first_child + quad];
                    float       distance_squared =
                        st::internal::euclidean_distance_squared_arr<CoordinateType, Rank>(
                            point.position, branch_data__.position);
                    if (distance_squared >= span_squared) {
                        update_data(branch_data__);
                        recurse[quad] = false;
                    } else if (branch_.is_terminal()) {
                        const auto& leaf = this->leaves_[branch_.index()];
                        for (uint64_t i = 0; i < branch_.size; ++i) {
                            update_data(
                                n_body_tree_data_external{leaf.coordinates[i], leaf.items[i].data});
                        }
                        recurse[quad] = false;
                    } else {
                        recurse[quad] = true;
                    }
                });
                st::internal::unroll_for<spatial_tree_type::BranchingFactor>([&](auto quad) {
                    if (recurse[quad]) {
                        update_recursively(branch.index_of_first_child + quad, span_squared / 4);
                    }
                });
            };
            update_recursively(0, tree_span_squared / 4);
            st::internal::unroll_for<Rank>([&](auto i) {
                point.velocity[i] += point.acceleration[i] * dt;
                point.position[i] += point.velocity[i] * dt;
            });
        }
    }

    std::vector<n_body_tree_data_external> branch_data_;
};

static std::vector<n_body_data> generate_galaxy(uint64_t                         N,
                                                std::array<CoordinateType, Rank> window_size) {
    CoordinateType span = std::numeric_limits<CoordinateType>::max();
    st::internal::unroll_for<Rank>([&](auto i) { span = std::min(span, window_size[i]); });
    const double                           galaxy_radius = span / 2;
    const double                           arm_tightness = 10.0 / span;
    static constexpr int                   arms = 2;
    static constexpr double                spread = 0.4;
    static constexpr double                max_velocity = 10.0;
    static constexpr double                radius_scale = 5.0;
    std::default_random_engine             gen;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::normal_distribution<double>       normal(0.0, spread);
    std::vector<n_body_data>               points;

    for (uint64_t i = 0; i < N; ++i) {
        int    arm = i % arms;
        double r = galaxy_radius * std::sqrt(uniform(gen));
        double theta = r * arm_tightness + (2 * M_PI / arms) * arm;
        theta += normal(gen);
        r += normal(gen);
        double v_orbit = max_velocity * (1.0 - std::exp(-r / radius_scale));

        n_body_data data;
        data.mass = 1;

        st::internal::unroll_for<Rank>([&](auto i) {
            data.position[i] = window_size[i] / 2;
            if (i == 0) {
                data.position[i] += r * std::cos(theta);
                data.velocity[i] = -v_orbit * std::sin(theta);
            } else if (i == 1) {
                data.position[i] += r * std::sin(theta);
                data.velocity[i] = v_orbit * std::cos(theta);
            } else {
                data.position[i] = 0;
                data.velocity[i] = 0;
            }
        });

        points.push_back(data);
    }

    return points;
}

int main(int argc, char** argv) {
    assert(argc == 3);
    const uint32_t number_of_points = std::atoi(argv[1]);
    const bool     draw = std::atoi(argv[2]);

    static constexpr uint64_t        WindowWidth = 1620;
    static constexpr uint64_t        WindowHeight = 1080;
    std::array<CoordinateType, Rank> window_size = {WindowWidth, WindowHeight};
    if constexpr (Rank == 3) {
        window_size[2] = WindowHeight;
    }
    std::vector<n_body_data> galaxy = generate_galaxy(number_of_points, window_size);

    InitWindow(WindowWidth, WindowHeight, "2-galaxy problem");
    SetTargetFPS(60);

    BeginDrawing();
    ClearBackground(RAYWHITE);
    EndDrawing();

    barnes_hut solver;
    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);

        solver.clear();
        solver.build(galaxy);
        solver.propagate();
        solver.update(galaxy, 0.05);

        if (draw) {
            for (auto point : galaxy) {
                DrawPixel(point.position[0], point.position[1], BLUE);
            }
        }

        DrawFPS(10, 10);
        EndDrawing();
    }

    CloseWindow();

    return 0;
}
