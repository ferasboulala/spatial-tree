#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>

#include "euclide.hh"
#include "helpers.hh"
#include "intervals.hh"

namespace st {
template <typename StorageType, typename CoordinateType = int, uint64_t RECURSION_CUTOFF = 4>
class QuadTree {
public:
    using __CoordinateType = CoordinateType;

    QuadTree(const BoundingBox<CoordinateType> &boundaries) : size_(0) {
        width_ = boundaries.bottom_x - boundaries.top_x;
        height_ = boundaries.top_y - boundaries.bottom_y;
        origin_x_ = (boundaries.top_x + boundaries.bottom_x) / 2;
        origin_y_ = (boundaries.top_y + boundaries.bottom_y) / 2;

        nodes_.resize(1);
    }

    ~QuadTree() = default;

    __always_inline uint64_t size() const { return size_; }
    __always_inline bool     empty() const { return size() == 0; }
    __always_inline void     clear() { nodes_.resize(1); }

    struct Iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = Iterator;
        using pointer = value_type *;
        using reference = value_type &;

        Iterator(const QuadTree<StorageType, CoordinateType> *tree,
                 uint64_t                                     node_index,
                 uint8_t                                      item_index)
            : tree_(tree), node_index_(node_index), item_index_(item_index) {}
        Iterator(const QuadTree<StorageType, CoordinateType> *tree) : tree_(tree) { end(); }
        ~Iterator() = default;

        __always_inline auto operator*() const { return tree_->operator()(node_index_, item_index_); }

        // TODO
        Iterator &operator++();

        __always_inline bool operator==(const Iterator &other) const {
            assert(tree_ == other.tree_);
            return node_index_ == other.index_ && item_index_ == other.item_index_;
        }

        __always_inline bool operator!=(const Iterator &other) const { return !(*this == other); }

    private:
        __always_inline void end() { node_index_ = NO_INDEX; }

        const QuadTree<StorageType, CoordinateType> *tree_;

        uint64_t node_index_;
        uint8_t  item_index_;
    };

    // TODO
    void walk(const std::function<void(Iterator)> func) const;

    __always_inline Iterator end() const { return Iterator(this); }
    __always_inline Iterator begin() const {
        if (empty()) {
            return end();
        }

        return Iterator(this, 0);
    }

    template <typename... Args>
    __always_inline std::pair<Iterator, bool> emplace(CoordinateType x,
                                                      CoordinateType y,
                                                      Args &&...args) {
        /// TODO: Check that coordinates fit the global bounding box.
        return emplace_recursively(
            0, height_, width_, origin_x_, origin_y_, x, y, std::forward<Args>(args)...);
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

    // TODO
    std::vector<Iterator> nearest(CoordinateType x, CoordinateType y) const;

    // TODO
    __always_inline void find(const BoundingBox<CoordinateType> &bbox,
                              std::function<void(Iterator)>      func) const;

    Iterator find(CoordinateType x, CoordinateType y) const {
        Iterator it = end();
        find({x, y, x, y}, [&](Iterator found) { it = found; });

        return it;
    }

private:
    static constexpr uint64_t NO_INDEX = std::numeric_limits<uint64_t>::max();
    static constexpr uint64_t NO_PARENT = std::numeric_limits<uint64_t>::max();
    // This is the equivalent of a namespace. Can't use a class enum because it does not let me
    // index into arrays.
    struct Quadrant {
        enum _Quadrant { NE = 0, NW = 1, SW = 2, SE = 3 };
    };

    struct NodeContent {
        // TODO: Check if having two arrays, one for coordinates and one for storage would improve
        // performance when manipulating coordinates.
        CoordinateType x;
        CoordinateType y;
        StorageType    storage;

        template <typename... Args>
        NodeContent(CoordinateType x_, CoordinateType y_, Args &&...args)
            : x(x_), y(y_), storage(std::forward<Args>(args)...) {}
    };

    struct NodeLeaves {
        std::array<NodeContent, RECURSION_CUTOFF> items;
        uint8_t                                   size;
    };

    struct Node {
        uint64_t parent;

        union {
            NodeLeaves              leaves;
            std::array<uint64_t, 4> children;
        };
        bool should_recurse;

        __always_inline Node() : should_recurse(false) { leaves.size = 0; }
    };

    __always_inline std::tuple<CoordinateType, CoordinateType, const StorageType &> operator()(
        uint64_t node_index, uint64_t item_index) const {
        assert(node_index < nodes_.size());
        const Node &node = nodes_[node_index];

        assert(!node.should_recurse);
        assert(item_index < node.leaves.size);

        const auto &[x, y, storage] = node.leaves.items[item_index];

        return std::tuple<CoordinateType, CoordinateType, const StorageType &>(x, y, storage);
    }

    /// TODO: Move to intervals.hh
    __always_inline auto belongs_to_quadrant(CoordinateType origin_x,
                                             CoordinateType origin_y,
                                             CoordinateType x,
                                             CoordinateType y) {
        return Quadrant::NE + ((y >= origin_y) * (x < origin_x)) * Quadrant::NW +
               (y < origin_y) * ((x < origin_x) * Quadrant::SW + (x >= origin_x) * Quadrant::SE);
    }

    __always_inline std::tuple<CoordinateType, CoordinateType, CoordinateType, CoordinateType>
                    child_location(uint64_t       quad,
                                   CoordinateType h,
                                   CoordinateType w,
                                   CoordinateType origin_x,
                                   CoordinateType origin_y) {
        static std::array<std::tuple<int, int>, 4> SIGNS;
        SIGNS[Quadrant::NE] = std::tuple<int, int>{1, 1};
        SIGNS[Quadrant::NW] = std::tuple<int, int>{-1, 1};
        SIGNS[Quadrant::SW] = std::tuple<int, int>{-1, -1};
        SIGNS[Quadrant::SE] = std::tuple<int, int>{1, -1};

        const auto [sign_x, sign_y] = SIGNS[quad];
        return {h / 2, w / 2, origin_x + sign_x * h / 2, origin_y + sign_y * w / 2};
    }

    template <typename... Args>
    __always_inline std::pair<Iterator, bool> emplace_recursively_helper(
        const std::array<uint64_t, 4> &children,
        CoordinateType                 h,
        CoordinateType                 w,
        CoordinateType                 origin_x,
        CoordinateType                 origin_y,
        CoordinateType                 x,
        CoordinateType                 y,
        Args &&...args) {
        assert(h || w);

        const auto selected_quadrant = belongs_to_quadrant(origin_x, origin_y, x, y);
        const auto [new_h, new_w, new_origin_x, new_origin_y] =
            child_location(selected_quadrant, h, w, origin_x, origin_y);

        return emplace_recursively(children[selected_quadrant],
                                   new_h,
                                   new_w,
                                   new_origin_x,
                                   new_origin_y,
                                   x,
                                   y,
                                   std::forward<Args>(args)...);
    }

    template <typename... Args>
    std::pair<Iterator, bool> emplace_recursively(uint64_t       node_index,
                                                  CoordinateType h,
                                                  CoordinateType w,
                                                  CoordinateType origin_x,
                                                  CoordinateType origin_y,
                                                  CoordinateType x,
                                                  CoordinateType y,
                                                  Args &&...args) {
        assert(node_index < nodes_.size());

        Node &node = nodes_[node_index];

        const bool should_split = node.leaves.size == RECURSION_CUTOFF;
        if (!node.should_recurse && !should_split) {
            uint64_t idx;
            bool     inserted = true;
            // TODO: Check if can optimize by removing the break and unrolling.
            for (idx = 0; idx < node.leaves.size; ++idx) {
                const auto &[x_, y_, storage_] = node.leaves.items[idx];
                if (x == x_ && y == y_) {
                    inserted = false;
                    break;
                }
            }

            if (inserted) {
                node.leaves.items[idx] = std::move(NodeContent(x, y, std::forward<Args>(args)...));
                ++node.leaves.size;
                ++size_;
            }

            return {Iterator(this, node_index, idx), inserted};
        }

        if (should_split) {
            // Subdivide the current node.
            std::array<uint64_t, 4> new_children;
            for (uint64_t i = 0; i < 4; ++i) {
                new_children[i] = nodes_.size();
                /// FIXME: Will not work if `StorageType` does not have a default constructor.
                nodes_.resize(nodes_.size() + 1);
            }

            /// Insert the nodes one after the other by selecting the appropriate quadrant every
            /// time.

            /// TODO: Optimization for the case where all children would be long to a deep
            /// square. Instead, find the smallest square that would contain them all and
            /// start inserting from there. Then, create all the parents that lead to it.
            for (auto &entry : node.leaves.items) {
                emplace_recursively_helper(
                    new_children, h, w, origin_x, origin_y, entry.x, entry.y, entry.storage);
            }

            node.should_recurse = true;
            node.children = new_children;
            size_ -= RECURSION_CUTOFF;
        }

        // Attempt a new insertion (could recurse again).
        /// FIXME: Use a while loop to avoid too much recursion.
        return emplace_recursively_helper(
            node.children, h, w, origin_x, origin_y, x, y, std::forward<Args>(args)...);
    }

    CoordinateType height_;
    CoordinateType width_;
    CoordinateType origin_x_;
    CoordinateType origin_y_;

    uint64_t size_;

    std::vector<Node> nodes_;
};
}  // namespace st