#include <raylib.h>

#include <iostream>
#include <random>

#include "spatial-tree.h"

static constexpr int Rank = 2;
static constexpr int MaximumLeafSize = 128;
using CoordinateType = float;
static constexpr CoordinateType R = 1;
using MassType = float;

struct n_body_data {
    std::array<CoordinateType, Rank> position;
    std::array<CoordinateType, Rank> velocity;
    MassType                         mass;
};

struct n_body_tree_data {
    std::array<CoordinateType, Rank> velocity;
    MassType                         mass;
};

class barnes_hut : public st::internal::
                       spatial_tree<CoordinateType, n_body_tree_data, Rank, MaximumLeafSize, 32> {
    using spatial_tree_type =
        st::internal::spatial_tree<CoordinateType, n_body_tree_data, Rank, MaximumLeafSize, 32>;

public:
    barnes_hut() = default;
    ~barnes_hut() = default;

    void build(const auto& points) {
        st::bounding_box<CoordinateType, Rank> boundary;
        boundary.stops.fill(std::numeric_limits<CoordinateType>::lowest() / 2);
        boundary.starts.fill(std::numeric_limits<CoordinateType>::max() / 2);

        for (const auto& point : points) {
            assert(point.mass);
            st::internal::unroll_for<Rank>([&](auto i) {
                boundary.starts[i] = std::min(point.position[i], boundary.starts[i]);
                boundary.stops[i] = std::max(point.position[i], boundary.stops[i]);
            });
        }

        this->reset(boundary);
        for (const auto& point : points) {
            auto inserted =
                this->emplace(point.position, n_body_tree_data{point.velocity, point.mass});
            assert(inserted.second);
        }
    }

    void propagate() {
        branch_data_.resize(this->branches_.size());
        for (int64_t i = this->branches_.size() - 1; i >= 0; --i) {
            const auto& branch = this->branches_[i];
            auto&       branch_data_entry = branch_data_[i];
            branch_data_entry.mass = 0;
            branch_data_entry.position = {0};

            if (branch.is_terminal()) {
                auto& leaf = this->leaves_[branch.index()];
                std::for_each(leaf.coordinates.begin(),
                              leaf.coordinates.begin() + leaf.size,
                              [&](auto position) {
                                  st::internal::unroll_for<Rank>([&](auto i) {
                                      branch_data_entry.position[i] += position[i];
                                  });
                              });
                std::for_each(leaf.items.begin(),
                              leaf.items.begin() + leaf.size,
                              [&](const auto& data) { branch_data_entry.mass += data.data.mass; });

                st::internal::unroll_for<Rank>(
                    [&](auto i) { branch_data_entry.position[i] /= branch.size; });
                branch_data_entry.mass /= branch.size;
            } else {
                st::internal::unroll_for<spatial_tree_type::BranchingFactor>([&](auto i) {
                    const auto& child_branch_data_entry =
                        branch_data_[branch.index_of_first_child + i];
                    st::internal::unroll_for<Rank>([&](auto i) {
                        branch_data_entry.position[i] += child_branch_data_entry.position[i];
                    });
                    branch_data_entry.mass += child_branch_data_entry.mass;
                });
                st::internal::unroll_for<Rank>([&](auto i) {
                    branch_data_entry.position[i] /= spatial_tree_type::BranchingFactor;
                });
                branch_data_entry.mass /= spatial_tree_type::BranchingFactor;
            }
        }
    }

    void update(auto& points, float dt) {
        for (auto& point : points) {
            const std::function<void(uint64_t, CoordinateType)> update_recursively =
                [&](uint64_t branch_index, CoordinateType width) {
                    const auto& branch = this->branches_[branch_index];
                    const auto& branch_data = branch_data_[branch_index];

                    // static constexpr float G = 6.67e-11;
                    static constexpr float G = 1e5;
                    static constexpr float EPSILON = 1e-6;

                    const auto update_data = [&](const n_body_data& data) {
                        float distance_squared =
                            st::internal::euclidean_distance_squared_arr<CoordinateType, Rank>(
                                point.position, data.position);
                        float distance_3_2 = std::pow(distance_squared, 1.5);
                        float coeff = G * data.mass / (distance_3_2 + EPSILON);
                        st::internal::unroll_for<Rank>([&](auto i) {
                            float acceleration = (data.position[i] - point.position[i]) * coeff;
                            point.position[i] += point.velocity[i] * dt;
                            point.velocity[i] += acceleration * dt;
                        });
                    };

                    // Replace distance with fast inverse square root.
                    auto distance_squared =
                        st::internal::euclidean_distance_squared_arr<CoordinateType, Rank>(
                            point.position, branch_data.position);

                    if (distance_squared >= R * width) {
                        update_data(branch_data);
                        return;
                    }

                    if (branch.is_terminal()) {
                        const auto& leaf = this->leaves_[branch.index()];
                        for (uint64_t i = 0; i < branch.size; ++i) {
                            update_data(n_body_data{leaf.coordinates[i],
                                                    leaf.items[i].data.velocity,
                                                    leaf.items[i].data.mass});
                        }
                        return;
                    }

                    st::internal::unroll_for<spatial_tree_type::BranchingFactor>([&](auto i) {
                        update_recursively(branch.index_of_first_child + branch_index, width / 2);
                    });

                    return;
                };
            update_recursively(0, this->boundary_.stops[0] - this->boundary_.starts[0]);
        }
    }

    std::vector<n_body_data> branch_data_;
};

// Function to compute orbital velocity based on distance (flat rotation curve)
double compute_orbital_velocity(double r, double max_velocity = 5.0, double radius_scale = 5.0) {
    // Flat rotation curve: velocity becomes constant beyond a certain radius
    return max_velocity * (1.0 - std::exp(-r / radius_scale));
}

std::vector<n_body_data> generate_points(uint64_t size, double span) {
    const double            galaxy_radius = span / 2;
    const double            arm_tightness = 10.0 / span;
    static constexpr int    arms = 3;
    static constexpr double spread = 0.2;

    std::default_random_engine             gen;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::normal_distribution<double>       normal(0.0, spread);

    std::vector<n_body_data> points;
    for (uint64_t i = 0; i < size; ++i) {
        int    arm = i % arms;
        double r = galaxy_radius * std::sqrt(uniform(gen));
        double theta = r * arm_tightness + (2 * M_PI / arms) * arm;
        theta += normal(gen);
        r += normal(gen);
        double v_orbit = compute_orbital_velocity(r);

        n_body_data data{0};
        data.mass = 1;

        st::internal::unroll_for<Rank>([&](auto i) {
            if (i == 0) {
                data.position[i] = span / 2 + r * std::cos(theta);
                data.velocity[i] = -v_orbit * std::sin(theta);
            } else if (i == 1) {
                data.position[i] = span / 2 + r * std::sin(theta);
                data.velocity[i] = v_orbit * std::cos(theta);
            } else {
                data.position[i] = span / 2;
            }
        });

        points.push_back(data);
    }
    std::array<CoordinateType, Rank> black_hole;
    st::internal::unroll_for<Rank>([&](auto i) { black_hole[i] = CoordinateType(span / 2); });
    points.push_back({black_hole, {0, 0}, MassType(10 * size)});

    return points;
}

int main(int argc, char** argv) {
    assert(argc == 2);
    const uint32_t number_of_points = std::atoi(argv[1]);
    barnes_hut     solver;

    static constexpr uint64_t WindowWidth = 800;
    static constexpr uint64_t WindowHeight = 800;
    std::vector<n_body_data>  entities = generate_points(number_of_points, WindowWidth);

    InitWindow(WindowWidth, WindowHeight, "n-body");
    SetTargetFPS(60);

    BeginDrawing();
    ClearBackground(RAYWHITE);
    EndDrawing();

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);

        solver.clear();
        solver.build(entities);
        solver.propagate();
        solver.update(entities, GetFrameTime());

        int i = 0;
        for (auto point : entities) {
            DrawCircle(
                point.position[0], point.position[1], i++ == entities.size() - 1 ? 10 : 1, BLUE);
        }

        DrawFPS(10, 10);
        EndDrawing();
    }

    CloseWindow();

    return 0;
}
