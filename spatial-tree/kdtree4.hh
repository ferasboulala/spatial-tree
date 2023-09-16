#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <tuple>
#include <vector>

#include "compilation.hh"
#include "euclide.hh"
#include "intervals.hh"

namespace st {

template <typename StorageType, typename CoordinateType = int>
class QuadTree {
public:
    QuadTree() = default;
    ~QuadTree() = default;

    __always_inline void     reserve(uint64_t size) { nodes_.reserve(size); }
    __always_inline uint64_t capacity() const { return nodes_.capacity(); }
    __always_inline uint64_t size() const { return nodes_.size(); }
    __always_inline bool     empty() const { return nodes_.empty(); }
    __always_inline void     clear() { nodes_.clear(); }
    __always_inline std::tuple<CoordinateType, CoordinateType, const StorageType &> operator[](
        uint64_t index) const {
        assert(index < nodes_.size());
        const Node &node = nodes_[index];

        return std::tuple<CoordinateType, CoordinateType, const StorageType &>(
            node.x, node.y, node.storage);
    }
    __always_inline std::tuple<CoordinateType, CoordinateType, StorageType &> operator[](
        uint64_t index) {
        assert(index < nodes_.size());
        Node &node = nodes_[index];

        return std::tuple<CoordinateType, CoordinateType, StorageType &>(
            node.x, node.y, node.storage);
    }

    __always_inline uint64_t index(const StorageType &value) const {
        const uint64_t offset = offsetof(Node, storage);
        const uint8_t *base_address = reinterpret_cast<uint8_t *>(&value);
        const uint8_t *node_address = base_address - offset;
        const auto     node = reinterpret_cast<const Node *>(node_address);

        return std::distance(nodes_.begin(), node);
    }

    template <typename... Args>
    StorageType &emplace(CoordinateType x, CoordinateType y, Args &&...args) {
        if (nodes_.empty()) {
            nodes_.emplace_back(x, y, std::forward<Args>(args)...);
        } else {
            emplace_recursively(0, x, y, std::forward<Args>(args)...);
        }

        return nodes_.back().storage;
    }

    __always_inline StorageType &insert(CoordinateType     x,
                                        CoordinateType     y,
                                        const StorageType &storage) {
        return emplace(x, y, storage);
    }

    std::tuple<CoordinateType, CoordinateType, StorageType &> nearest(CoordinateType x,
                                                                      CoordinateType y) {
        assert(!empty());
        Node          *nearest_node = &nodes_.front();
        CoordinateType nearest_distance_squared =
            euclidean_distance_squared(x, nearest_node->x, y, nearest_node->y);
        nearest_recursive(0, {x, y}, &nearest_node, nearest_distance_squared);

        return std::tuple<CoordinateType, CoordinateType, StorageType &>(
            nearest_node->x, nearest_node->y, nearest_node->storage);
    }

    std::vector<std::tuple<CoordinateType, CoordinateType, StorageType &>> find(
        const BoundingBox<CoordinateType> &bbox) {
        std::vector<std::tuple<CoordinateType, CoordinateType, StorageType &>> ret;
        if (!empty()) {
            find_recursively(0, bbox, ret);
        }

        return ret;
    }

private:
    static constexpr uint64_t EMPTY_INDEX = std::numeric_limits<uint64_t>::max();
    struct Node {
        CoordinateType          x;
        CoordinateType          y;
        std::array<uint64_t, 4> quadrants;
        StorageType             storage;

        template <typename... Args>
        Node(CoordinateType x_, CoordinateType y_, Args &&...args)
            : x(x_), y(y_), storage(std::forward<Args>(args)...) {
            std::fill_n(quadrants.begin(), 4, EMPTY_INDEX);
        }
        Node(CoordinateType x_, CoordinateType y_, const StorageType &storage_)
            : x(x_), y(y_), storage(storage_) {
            std::fill_n(quadrants.begin(), 4, EMPTY_INDEX);
        }
    };

    // NE: 0, NW: 1, SW: 2, SE: 3
    __always_inline uint64_t belongs_to_quadrant(const Node    &parent,
                                                 CoordinateType x,
                                                 CoordinateType y) {
        // Much faster than the branching variant. Very strange that the compiler does not do it.
        // Thank you Fedor Pikus.
        return ((y >= parent.y) * (x < parent.x)) +
               (y < parent.y) * ((x < parent.x) * 2 + (x >= parent.x) * 3);
    }

    template <typename... Args>
    void emplace_recursively(uint64_t       node_index,
                             CoordinateType x,
                             CoordinateType y,
                             Args &&...args) {
        assert(node_index != EMPTY_INDEX);
        assert(node_index < nodes_.size());

        Node &parent = nodes_[node_index];
        if (parent.x == x && parent.y == y) {
            parent.storage = std::move(StorageType(std::forward<Args>(args)...));
            return;
        }

        const uint64_t selected_quadrant = belongs_to_quadrant(parent, x, y);
        uint64_t      &child_index = parent.quadrants[selected_quadrant];
        if (child_index == EMPTY_INDEX) {
            child_index = nodes_.size();
            nodes_.emplace_back(x, y, std::forward<Args>(args)...);
        } else {
            emplace_recursively(child_index, x, y, std::forward<Args>(args)...);
        }
    }

    static __always_inline CoordinateType safe_delta(CoordinateType a, CoordinateType b) {
        return (a > b) * (a - b) + (a < b) * (b - a);
    }

    void nearest_recursive(uint64_t                             node_index,
                           const std::array<CoordinateType, 2> &xy,
                           Node                               **nearest_node,
                           CoordinateType                      &nearest_distance_squared) {
        if (node_index == EMPTY_INDEX) {
            return;
        }

        Node                &parent = nodes_[node_index];
        const CoordinateType distance_squared =
            euclidean_distance_squared(xy[0], parent.x, xy[1], parent.y);
        if (distance_squared < nearest_distance_squared) {
            nearest_distance_squared = distance_squared;
            *nearest_node = &parent;
        }

        const CoordinateType dx = safe_delta(xy[0], parent.x);
        const CoordinateType dxdx = dx * dx;
        const CoordinateType dy = safe_delta(xy[1], parent.y);
        const CoordinateType dydy = dy * dy;
        const CoordinateType dxdxdydy = dxdx + dydy;

        const uint64_t selected_quadrant = belongs_to_quadrant(parent, xy[0], xy[1]);
        const uint64_t previous_quadrant = (selected_quadrant + 4 - 1) % 4;
        const uint64_t next_quadrant = (selected_quadrant + 4 + 1) % 4;
        const uint64_t diagonal_quadrant = (selected_quadrant + 4 + 2) % 4;
        const uint64_t x_quadrant = previous_quadrant % 2 ? previous_quadrant : next_quadrant;
        const uint64_t y_quadrant = previous_quadrant % 2 ? next_quadrant : previous_quadrant;

        nearest_recursive(
            parent.quadrants[selected_quadrant], xy, nearest_node, nearest_distance_squared);
        if (dxdx < nearest_distance_squared) {
            nearest_recursive(
                parent.quadrants[x_quadrant], xy, nearest_node, nearest_distance_squared);
        }
        if (dydy < nearest_distance_squared) {
            nearest_recursive(
                parent.quadrants[y_quadrant], xy, nearest_node, nearest_distance_squared);
        }
        if (dxdxdydy < nearest_distance_squared) {
            nearest_recursive(
                parent.quadrants[diagonal_quadrant], xy, nearest_node, nearest_distance_squared);
        }
    }

    static __always_inline std::array<BoundingBox<CoordinateType>, 4> neighborhood(
        CoordinateType x, CoordinateType y) {
        return {BoundingBox<CoordinateType>(x,
                                            std::numeric_limits<CoordinateType>::max(),
                                            std::numeric_limits<CoordinateType>::max(),
                                            y),
                BoundingBox<CoordinateType>(std::numeric_limits<CoordinateType>::min(),
                                            std::numeric_limits<CoordinateType>::max(),
                                            x,
                                            y),
                BoundingBox<CoordinateType>(std::numeric_limits<CoordinateType>::min(),
                                            y,
                                            x,
                                            std::numeric_limits<CoordinateType>::min()),
                BoundingBox<CoordinateType>(x,
                                            y,
                                            std::numeric_limits<CoordinateType>::max(),
                                            std::numeric_limits<CoordinateType>::min())};
    }

    void find_recursively(
        uint64_t                                                                node_index,
        const BoundingBox<CoordinateType>                                      &target_bbox,
        std::vector<std::tuple<CoordinateType, CoordinateType, StorageType &>> &ret) {
        Node &node = nodes_[node_index];
        if (is_inside_bounding_box(node.x, node.y, target_bbox)) {
            ret.emplace_back(node.x, node.y, std::ref(node.storage));
        }

        const std::array<BoundingBox<CoordinateType>, 4> neighbors = neighborhood(node.x, node.y);

// Compilers do not always unroll the following loop.
#ifdef GNU_COMPILER
#pragma GCC unroll 4
#endif
#ifdef CLANG_COMPILER
#pragma clang loop unroll_count(4)
#endif

        for (uint64_t i = 0; i < 4; ++i) {
            const uint64_t child_index = node.quadrants[i];
            if (child_index != EMPTY_INDEX) {
                const BoundingBox<CoordinateType> &neighbor_bbox = neighbors[i];
                if (bounding_boxes_overlap(neighbor_bbox, target_bbox)) {
                    find_recursively(child_index, target_bbox, ret);
                }
            }
        }
    }

    std::vector<Node> nodes_;
};

}  // namespace st
