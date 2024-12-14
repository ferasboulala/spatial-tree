#include <iostream>

#include "spatial-tree.h"

template <typename CoordinateType = float,
          typename MassType = float,
          uint64_t Rank = 2,
          uint16_t MaximumLeafSize = 128>
struct barnes_hut {
    barnes_hut() = default;
    ~barnes_hut() = default;
    static constexpr uint64_t BranchingFactor = st::internal::pow(2, Rank);

    void build(const auto& points) {
        branches_.clear();
        branches_.resize(1);
        leaves_.clear();
        leaves_.resize(1);
        freed_leaves_.clear();
        presence_.clear();

        st::bounding_box<CoordinateType, Rank> boundary_;
        boundary_ = decltype(boundary_)();
        for (const auto& [position, _] : points) {
            st::internal::unroll_for<Rank>([&](auto i) {
                boundary_.starts[i] = std::max(position[i], boundary_.starts[i]) - 1;
                boundary_.stops[i] = std::min(position[i], boundary_.stops[i]) + 1;
            });
        }

        std::cout << "Boundaries are: \n";
        st::internal::unroll_for<Rank>(
            [&](auto i) { std::cout << boundary_.starts[i] << " " << boundary_.stops[i] << "\n"; });

        for (const auto& [position, mass] : points) {
            auto inserted = emplace(position, mass);
            assert(inserted);
        }
    }

    void propagate();
    void update();

    struct tree_leaf {
        uint16_t                                                      size;
        std::array<std::array<CoordinateType, Rank>, MaximumLeafSize> coordinates;
        std::array<MassType, MaximumLeafSize>                         masses;

        inline tree_leaf() { reset(); }
        inline tree_leaf(const tree_leaf& other) : size(other.size) {
            st::internal::unroll_for<4>(uint16_t(0), size, [&](auto i) {
                coordinates[i] = other.coordinates[i];
                masses[i] = other.masses[i];
            });
        }
        inline tree_leaf(tree_leaf&& other) : size(other.size) {
            st::internal::unroll_for<4>(uint16_t(0), size, [&](auto i) {
                coordinates[i] = other.coordinates[i];
                masses[i] = other.masses[i];
            });
        }
        inline void reset() { size = 0; }
    };

    struct tree_branch {
        uint32_t size;
        int32_t  index_of_first_child;

        inline tree_branch() { reset(); }
        inline bool     is_terminal() const { return index_of_first_child < 0; }
        inline uint32_t index() const { return index_of_first_child * -1 - 1; }
        inline void     index(uint32_t idx) { index_of_first_child = int32_t(idx + 1) * -1; }
        inline void     reset() {
            size = 0;
            index_of_first_child = -1;
        }
    };

    inline bool emplace(std::array<CoordinateType, Rank> point, MassType mass) {
        assert(boundary_.contains(point));
        auto [_, inserted] = presence_.emplace(point);
        if (!inserted) [[unlikely]] {
            return false;
        }

        emplace_recursively(0, boundary_, point, mass);

        return true;
    }

    inline void emplace_recursively_helper(uint32_t index_of_first_child,
                                           const st::bounding_box<CoordinateType, Rank>& boundary,
                                           std::array<CoordinateType, Rank>              point,
                                           MassType                                      mass) {
        const auto [new_boundary, selected_quad] = boundary.recurse(point);
        return emplace_recursively(index_of_first_child + selected_quad, new_boundary, point, mass);
    }

    void emplace_recursively(uint32_t                                      branch_index,
                             const st::bounding_box<CoordinateType, Rank>& _boundary,
                             std::array<CoordinateType, Rank>              point,
                             MassType                                      mass) {
        st::bounding_box<CoordinateType, Rank> boundary = _boundary;
        while (true) {
            assert(_boundary.contains(point));
            assert(branch_index < branches_.size());

            tree_branch& branch = branches_[branch_index];
            ++branch.size;
            if (branch.is_terminal()) [[unlikely]] {
                break;
            }

            const auto [new_boundary, selected_quad] = boundary.recurse(point);
            branch_index = branch.index_of_first_child + selected_quad;
            boundary = new_boundary;
        }

        tree_branch& terminal_branch = branches_[branch_index];
        assert(terminal_branch.index() < leaves_.size());
        tree_leaf& leaf = leaves_[terminal_branch.index()];

        if (leaf.size < MaximumLeafSize) [[likely]] {
            st::internal::set<CoordinateType, Rank>(leaf.coordinates[leaf.size], point);
            leaf.masses[leaf.size] = mass;
            return;
        }

        uint64_t index_of_freed_leaf = terminal_branch.index();
        uint64_t new_index_of_first_child = branches_.size();
        branches_.resize(branches_.size() + BranchingFactor);
        assert((new_index_of_first_child - 1) % BranchingFactor == 0);
        tree_branch& branch = branches_[branch_index];
        branch.index_of_first_child = new_index_of_first_child;
        branch.size = MaximumLeafSize + 1;

        st::internal::unroll_for<BranchingFactor>([&](auto i) {
            tree_branch& terminal_branch = branches_[new_index_of_first_child + i];
            if (!freed_leaves_.empty()) {
                terminal_branch.index(freed_leaves_.back());
                assert(leaves_[terminal_branch.index()].size == 0);
                freed_leaves_.pop_back();
            } else {
                terminal_branch.index(leaves_.size());
                leaves_.resize(leaves_.size() + 1);
            }
        });

        st::internal::unroll_for<MaximumLeafSize>([&](auto i) {
            emplace_recursively_helper(new_index_of_first_child,
                                       boundary,
                                       leaves_[index_of_freed_leaf].coordinates[i],
                                       leaves_[index_of_freed_leaf].masses[i]);
        });

        leaves_[index_of_freed_leaf].reset();
        freed_leaves_.push_back(index_of_freed_leaf);

        return emplace_recursively_helper(new_index_of_first_child, boundary, point, mass);
    }

    st::bounding_box<CoordinateType, Rank> boundary_;
    std::vector<tree_branch>               branches_;
    std::vector<tree_leaf>                 leaves_;
    std::vector<uint32_t>                  freed_leaves_;

    robin_hood::unordered_set<std::array<CoordinateType, Rank>,
                              st::internal::point_hash<CoordinateType, Rank>,
                              st::internal::point_equal<CoordinateType, Rank>>
        presence_;
};

int main(int argc, char** argv) {
    assert(argc == 2);
    const uint32_t                                      number_of_points = std::atoi(argv[1]);
    barnes_hut                                          solver;
    std::vector<std::pair<std::array<float, 2>, float>> entities;
    entities.push_back({{0, 0}, 1.0});
    solver.build(entities);
    return 0;
}
