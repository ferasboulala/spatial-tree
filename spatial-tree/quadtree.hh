#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <tuple>
#include <vector>

#include "euclide.hh"
#include "helpers.hh"
#include "intervals.hh"

/// Far fetched optimization for Integer based coodinate types:
/// Can stop recursing when h * w <= MAXIMUM_NODE_SIZE

namespace st {
template <typename StorageType, typename CoordinateType = int, uint8_t MAXIMUM_NODE_SIZE = 4>
class QuadTree {
public:
    static_assert(MAXIMUM_NODE_SIZE < 255);

    using __CoordinateType = CoordinateType;

    QuadTree(const BoundingBox<CoordinateType> &boundaries) : boundaries_(boundaries) { clear(); }

    ~QuadTree() = default;

    __always_inline uint64_t size() const { return size_; }
    __always_inline bool     empty() const { return size() == 0; }
    __always_inline void     clear() {
        nodes_.resize(1);
        nodes_.front() = Node();
        size_ = 0;
    }

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

        __always_inline auto operator*() const {
            return tree_->operator()(node_index_, item_index_);
        }

        // TODO
        Iterator &operator++() {
            assert(node_index_ != NO_INDEX);

            const Node *node = &tree_->nodes_[node_index_];
            if (!node->is_a_branch && ++item_index_ < node->leaves.size) {
                return *this;
            }

            do {
                ++node_index_;
                node = &tree_->nodes_[node_index_];
                item_index_ = 0;
                /// TODO: Add a Node::empty function to use here.
            } while ((node->is_a_branch || !node->leaves.size) &&
                     node_index_ < tree_->nodes_.size());

            if (node_index_ == tree_->nodes_.size()) end();

            return *this;
        }

        __always_inline bool operator==(const Iterator &other) const {
            assert(tree_ == other.tree_);
            return node_index_ == other.node_index_ && item_index_ == other.item_index_;
        }

        __always_inline bool operator!=(const Iterator &other) const { return !(*this == other); }

    private:
        __always_inline void end() {
            node_index_ = NO_INDEX;
            item_index_ = 0;
        }

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

        Iterator ret(this, 0, 0);

        const Node &first_node = nodes_.front();

        if (!first_node.is_a_branch && first_node.leaves.size) return ret;

        return ++ret;
    }

    template <typename... Args>
    __always_inline std::pair<Iterator, bool> emplace(CoordinateType x,
                                                      CoordinateType y,
                                                      Args &&...args) {
        /// TODO: Check that coordinates fit the global bounding box.
        return emplace_recursively(0, boundaries_, x, y, std::forward<Args>(args)...);
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
        /// Can be faster by checking a point.
        find({x, y, x, y}, [&](Iterator found) { it = found; });

        return it;
    }

private:
    static constexpr uint64_t NO_INDEX = std::numeric_limits<uint64_t>::max();
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
        std::array<NodeContent, MAXIMUM_NODE_SIZE> items;
        uint8_t                                    size;
    };

    struct Node {
        union {
            NodeLeaves              leaves;
            std::array<uint64_t, 4> children;
        };
        bool is_a_branch;

        __always_inline Node() : is_a_branch(false) { leaves.size = 0; }
    };

    __always_inline std::tuple<CoordinateType, CoordinateType, const StorageType &> operator()(
        uint64_t node_index, uint64_t item_index) const {
        assert(node_index < nodes_.size());
        const Node &node = nodes_[node_index];

        assert(!node.is_a_branch);
        assert(item_index < node.leaves.size);

        const auto &[x, y, storage] = node.leaves.items[item_index];

        return std::tuple<CoordinateType, CoordinateType, const StorageType &>(x, y, storage);
    }

    /// TODO: Move to intervals.hh
    /// TODO: See if BoundingBox is sent through registers. Check what happens to perf if changed.
    static __always_inline auto belongs_to_quadrant(const BoundingBox<CoordinateType> &boundaries,
                                                    CoordinateType                     x,
                                                    CoordinateType                     y) {
        CoordinateType origin_x = (boundaries.top_x + boundaries.bottom_x) / 2;
        CoordinateType origin_y = (boundaries.top_y + boundaries.bottom_y) / 2;

        return Quadrant::NE + ((y >= origin_y) * (x < origin_x)) * Quadrant::NW +
               (y < origin_y) * ((x < origin_x) * Quadrant::SW +
                                 (y < origin_y) * (x >= origin_x) * Quadrant::SE);
    }

    static __always_inline BoundingBox<CoordinateType> compute_new_boundaries(
        uint64_t quad, const BoundingBox<CoordinateType> &boundaries) {
        static std::array<std::tuple<bool, bool>, 4> factors;
        factors[Quadrant::NE] = std::tuple<bool, bool>{1, 1};
        factors[Quadrant::NW] = std::tuple<bool, bool>{0, 1};
        factors[Quadrant::SW] = std::tuple<bool, bool>{0, 0};
        factors[Quadrant::SE] = std::tuple<bool, bool>{1, 0};

        const auto [fx, fy] = factors[quad];

        /// TODO: Move this logic to BoundingBox::
        CoordinateType width = (boundaries.bottom_x - boundaries.top_x);
        CoordinateType height = (boundaries.top_y - boundaries.bottom_y);
        CoordinateType new_width = width / 2;
        CoordinateType new_height = height / 2;
        CoordinateType width_remainder = width - 2 * new_width;
        CoordinateType height_remainder = height - 2 * new_height;

        /// TODO: Combine this logic with belongs_to_quadrant because it must match.
        CoordinateType top_x = boundaries.top_x + fx * new_width;
        CoordinateType top_y = boundaries.top_y - (1 - fy) * new_height;

        CoordinateType bottom_x = boundaries.bottom_x - (1 - fx) * (new_width + width_remainder);
        CoordinateType bottom_y = boundaries.bottom_y + fy * (new_height + height_remainder);

        return {top_x, top_y, bottom_x, bottom_y};
    }

    template <typename... Args>
    __always_inline std::pair<Iterator, bool> emplace_recursively_helper(
        const std::array<uint64_t, 4>     &children,
        const BoundingBox<CoordinateType> &boundaries,
        CoordinateType                     x,
        CoordinateType                     y,
        Args &&...args) {
        const auto selected_quadrant = belongs_to_quadrant(boundaries, x, y);
        const auto new_boundaries = compute_new_boundaries(selected_quadrant, boundaries);

        return emplace_recursively(
            children[selected_quadrant], new_boundaries, x, y, std::forward<Args>(args)...);
    }

    template <typename... Args>
    std::pair<Iterator, bool> emplace_recursively(uint64_t                           node_index,
                                                  const BoundingBox<CoordinateType> &boundaries,
                                                  CoordinateType                     x,
                                                  CoordinateType                     y,
                                                  Args &&...args) {
        assert(node_index < nodes_.size());

        Node &node = nodes_[node_index];

        if (node.is_a_branch) {
            return emplace_recursively_helper(
                node.children, boundaries, x, y, std::forward<Args>(args)...);
        }

        if (node.leaves.size) {
            auto beg = node.leaves.items.begin();
            auto end = node.leaves.items.begin() + node.leaves.size;
            auto it = std::find_if(beg, end, [&](auto &entry) {
                const auto &[x_, y_, storage_] = entry;
                return x == x_ && y == y_;
            });
            if (it != end) {
                return {Iterator(this, node_index, std::distance(beg, it)), false};
            }
        }

        const bool node_is_full = node.leaves.size == MAXIMUM_NODE_SIZE;
        if (!node_is_full) {
            node.leaves.items[node.leaves.size++] =
                std::move(NodeContent(x, y, std::forward<Args>(args)...));
            ++size_;
            return {Iterator(this, node_index, node.leaves.size - 1), true};
        }

        // Subdivide the current node.
        std::array<uint64_t, 4> new_children;
        std::iota(new_children.begin(), new_children.end(), nodes_.size());
        nodes_.resize(nodes_.size() + 4);

        /// Insert the nodes one after the other by selecting the appropriate quadrant every
        /// time.
        /// TODO: Optimization for the case where all children would be long to a deep
        /// square. Instead, find the smallest square that would contain them all and
        /// start inserting from there. Then, create all the parents that lead to it.
        Node &node_as_branch = nodes_[node_index];
        for (auto &entry : node_as_branch.leaves.items) {
            emplace_recursively_helper(new_children, boundaries, entry.x, entry.y, entry.storage);
        }

        node_as_branch.is_a_branch = true;
        node_as_branch.children = new_children;
        size_ -= MAXIMUM_NODE_SIZE;

        // Attempt a new insertion (could recurse again).
        return emplace_recursively_helper(
            new_children, boundaries, x, y, std::forward<Args>(args)...);
    }

    const BoundingBox<CoordinateType> boundaries_;
    uint64_t                          size_;
    std::vector<Node>                 nodes_;
};
}  // namespace st