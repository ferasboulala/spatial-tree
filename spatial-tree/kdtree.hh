#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <limits>
#include <tuple>
#include <vector>

#include "euclide.hh"
#include "intervals.hh"

namespace st {
template <typename StorageType, typename CoordinateType = int>
class KDTree {
public:
    // For testing purposes.
    using __CoordinateType = CoordinateType;

    KDTree() = default;
    ~KDTree() = default;

    __always_inline void reserve(uint64_t size) {
        nodes_.reserve(size);
        storage_.reserve(size);
    }
    __always_inline uint64_t capacity() const { return nodes_.capacity(); }
    __always_inline uint64_t size() const { return nodes_.size(); }
    __always_inline bool     empty() const { return nodes_.empty(); }
    __always_inline void     clear() {
        nodes_.clear();
        storage_.clear();
    }

    struct Iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = Iterator;
        using pointer = value_type *;
        using reference = value_type &;

        Iterator(const KDTree<StorageType, CoordinateType> *tree, uint64_t index)
            : tree_(tree), index_(index), next_is_left_(true) {}
        Iterator(const KDTree<StorageType, CoordinateType> *tree) : tree_(tree) { end(); }
        ~Iterator() = default;

        __always_inline auto operator*() const { return tree_->operator[](index_); }

        Iterator &operator++() {
            assert(index_ != NO_INDEX);
            const auto visit = [&](int index, int side) {
                const Node &node = tree_->nodes_[index];
                if (node.children[side] != NO_INDEX) {
                    index_ = node.children[side];
                    next_is_left_ = true;
                    return true;
                }
                return false;
            };

            if (next_is_left_ && visit(index_, 0)) {
                return *this;
            }

            if (visit(index_, 1)) {
                return *this;
            }

            const Node *node = &tree_->nodes_[index_];
            while (node->parent != NO_PARENT) {
                const Node &parent_node = tree_->nodes_[node->parent];
                if (parent_node.children[0] == index_ && visit(node->parent, 1)) {
                    return *this;
                }
                index_ = node->parent;
                node = &parent_node;
            }

            end();

            return *this;
        }

        __always_inline bool operator==(const Iterator &other) const {
            assert(tree_ == other.tree_);
            return index_ == other.index_;
        }

        __always_inline bool operator!=(const Iterator &other) const { return !(*this == other); }

    private:
        __always_inline void end() { index_ = NO_INDEX; }

        const KDTree<StorageType, CoordinateType> *tree_;

        uint64_t index_;
        bool     next_is_left_;
    };

    void walk(const std::function<void(Iterator)> func) const {
        if (empty()) {
            return;
        }

        const std::function<void(uint64_t)> recurse = [&](uint64_t index) {
            if (index == NO_INDEX) {
                return;
            }

            func(Iterator(this, index));
            recurse(nodes_[index].children[0]);
            recurse(nodes_[index].children[1]);
        };

        recurse(0);
    }

    __always_inline Iterator end() const { return Iterator(this); }
    __always_inline Iterator begin() const {
        if (empty()) {
            return end();
        }

        return Iterator(this, 0);
    }

    template <typename... Args>
    std::pair<Iterator, bool> emplace(CoordinateType x, CoordinateType y, Args &&...args) {
        if (nodes_.empty()) {
            nodes_.emplace_back(x, y, 0);
            storage_.emplace_back(std::forward<Args>(args)...);
            return {Iterator(this, 0), true};
        }

        return emplace_recursively(0, {x, y}, std::forward<Args>(args)...);
    }

    __always_inline std::pair<Iterator, bool> insert(CoordinateType x,
                                                     CoordinateType y,
                                                     StorageType  &&storage) {
        return emplace(x, y, storage);
    }

    __always_inline std::pair<Iterator, bool> insert(CoordinateType     x,
                                                     CoordinateType     y,
                                                     const StorageType &storage) {
        return emplace(x, y, storage);
    }

    std::vector<Iterator> nearest(CoordinateType x, CoordinateType y) const {
        std::vector<Iterator> ret;

        if (empty()) {
            return ret;
        }

        CoordinateType nearest_distance_squared = std::numeric_limits<CoordinateType>::max();
        nearest_recursive(0, {x, y}, nearest_distance_squared, ret);

        return ret;
    }

    __always_inline void find(const BoundingBox<CoordinateType> &bbox,
                              std::function<void(Iterator)>      func) const {
        if (!empty()) {
            find_recursively(0, bbox, func);
        }
    }

    Iterator find(CoordinateType x, CoordinateType y) const {
        Iterator it = end();
        find({x, y, x, y}, [&](Iterator found) { it = found; });

        return it;
    }

private:
    static constexpr uint64_t NO_INDEX = std::numeric_limits<uint64_t>::max();
    static constexpr uint64_t NO_PARENT = std::numeric_limits<uint64_t>::max();

    struct Node {
        std::array<CoordinateType, 2> coordinates;
        // 0 is left, 1 is right
        std::array<uint64_t, 2> children;
        uint64_t                parent;

        // 0 is x, 1 is y
        uint8_t dimension;

        Node(CoordinateType x, CoordinateType y, uint8_t dim, uint64_t p = NO_PARENT)
            : parent(p), dimension(dim) {
            coordinates[0] = x;
            coordinates[1] = y;
            children[0] = NO_INDEX;
            children[1] = NO_INDEX;
        }
    };

    __always_inline std::tuple<CoordinateType, CoordinateType, const StorageType &> operator[](
        uint64_t index) const {
        assert(index < nodes_.size());
        const Node        &node = nodes_[index];
        const StorageType &storage = storage_[index];

        return std::tuple<CoordinateType, CoordinateType, const StorageType &>(
            node.coordinates[0], node.coordinates[1], storage);
    }

    template <typename... Args>
    std::pair<Iterator, bool> emplace_recursively(uint64_t                             node_index,
                                                  const std::array<CoordinateType, 2> &xy,
                                                  Args &&...args) {
        assert(node_index != NO_INDEX);
        assert(node_index < nodes_.size());

        Node &parent = nodes_[node_index];
        if (parent.coordinates[0] == xy[0] && parent.coordinates[1] == xy[1]) {
            return {Iterator(this, node_index), false};
        }

        const uint64_t left_or_right = xy[parent.dimension] > parent.coordinates[parent.dimension];
        uint64_t      &child_index = parent.children[left_or_right];
        if (child_index == NO_INDEX) {
            child_index = nodes_.size();
            nodes_.emplace_back(xy[0], xy[1], 1 - parent.dimension, node_index);
            storage_.emplace_back(std::forward<Args>(args)...);
            return {Iterator(this, nodes_.size() - 1), true};
        } else {
            return emplace_recursively(child_index, xy, std::forward<Args>(args)...);
        }
    }

    void nearest_recursive(uint64_t                             node_index,
                           const std::array<CoordinateType, 2> &xy,
                           CoordinateType                      &nearest_distance_squared,
                           std::vector<Iterator>               &nearest_points) const {
        if (node_index == NO_INDEX) {
            return;
        }

        const Node          &parent = nodes_[node_index];
        const CoordinateType distance_squared =
            euclidean_distance_squared(xy[0], parent.coordinates[0], xy[1], parent.coordinates[1]);
        if (distance_squared < nearest_distance_squared) {
            nearest_distance_squared = distance_squared;
            nearest_points.clear();
            nearest_points.push_back(Iterator(this, node_index));
        } else if (distance_squared == nearest_distance_squared) {
            nearest_points.push_back(Iterator(this, node_index));
        }

        const uint64_t left_or_right = xy[parent.dimension] > parent.coordinates[parent.dimension];
        const uint64_t child_index = parent.children[left_or_right];

        nearest_recursive(child_index, xy, nearest_distance_squared, nearest_points);

        const CoordinateType lower_bound_distance_squared =
            euclidean_distance_squared(xy[parent.dimension], parent.coordinates[parent.dimension]);
        if (lower_bound_distance_squared < nearest_distance_squared) {
            nearest_recursive(
                parent.children[1 - left_or_right], xy, nearest_distance_squared, nearest_points);
        }
    }

    void find_recursively(uint64_t                           node_index,
                          const BoundingBox<CoordinateType> &bbox,
                          std::function<void(Iterator)>      func) const {
        if (node_index == NO_INDEX) {
            return;
        }

        const Node &parent = nodes_[node_index];
        if (is_inside_bounding_box(parent.coordinates[0], parent.coordinates[1], bbox)) {
            func(Iterator(this, node_index));
        }

        const auto coords = parent.dimension
                                ? std::array<CoordinateType, 2>{bbox.bottom_y, bbox.top_y}
                                : std::array<CoordinateType, 2>{bbox.top_x, bbox.bottom_x};
        if (coords[1] < parent.coordinates[parent.dimension]) {
            find_recursively(parent.children[0], bbox, func);
        } else if (coords[0] > parent.coordinates[parent.dimension]) {
            find_recursively(parent.children[1], bbox, func);
        } else {
            find_recursively(parent.children[0], bbox, func);
            find_recursively(parent.children[1], bbox, func);
        }
    }

    std::vector<Node>        nodes_;
    std::vector<StorageType> storage_;
};
}  // namespace st