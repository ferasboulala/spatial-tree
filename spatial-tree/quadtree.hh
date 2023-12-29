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

namespace st {
template <typename StorageType = void, typename CoordinateType = int, uint8_t MAXIMUM_NODE_SIZE = 32>
class QuadTree {
public:
    static_assert(MAXIMUM_NODE_SIZE > 0, "Maximum node size must be greater than 1");

    QuadTree()
        // Dividing by 2 to be able to do (upper bound +/- lower bound) without overflowing.
        : boundaries_(std::numeric_limits<CoordinateType>::lowest() / 2,
                      std::numeric_limits<CoordinateType>::max() / 2,
                      std::numeric_limits<CoordinateType>::max() / 2,
                      std::numeric_limits<CoordinateType>::lowest() / 2) {
        clear();
    }
    QuadTree(const BoundingBox<CoordinateType> &boundaries) : boundaries_(boundaries) {
        assert(boundaries_.top_x >= std::numeric_limits<CoordinateType>::lowest() / 2);
        assert(boundaries_.top_y <= std::numeric_limits<CoordinateType>::max() / 2);
        assert(boundaries_.bottom_x <= std::numeric_limits<CoordinateType>::max() / 2);
        assert(boundaries_.bottom_y >= std::numeric_limits<CoordinateType>::lowest() / 2);
        clear();
    }

    ~QuadTree() = default;

    __always_inline uint64_t size() const { return size_; }
    __always_inline bool     empty() const { return size() == 0; }
    __always_inline void     clear() {
        nodes_.resize(1);
        nodes_.front().reset();
        size_ = 0;
    }

    struct Iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = Iterator;
        using pointer = value_type *;
        using reference = value_type &;

        Iterator(const QuadTree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> *tree,
                 uint64_t                                                        node_index,
                 uint8_t                                                         item_index)
            : tree_(tree), node_index_(node_index), item_index_(item_index) {}
        Iterator(const QuadTree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> *tree)
            : tree_(tree) {
            end();
        }
        ~Iterator() = default;

        __always_inline auto operator*() const {
            return tree_->operator()(node_index_, item_index_);
        }

        Iterator &operator++() {
            assert(node_index_ != Node::NO_INDEX);
            assert(node_index_ < tree_->nodes_.size());

            const Node *node = &tree_->nodes_[node_index_];
            if (node->is_a_leaf() && ++item_index_ < node->leaves.size) {
                return *this;
            }

            do {
                ++node_index_;
                node = &tree_->nodes_[node_index_];
                item_index_ = 0;
            } while ((node->is_a_branch() || !node->leaves.size) &&
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
            node_index_ = Node::NO_INDEX;
            item_index_ = 0;
        }

        const QuadTree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> *tree_;

        uint64_t node_index_;
        uint8_t  item_index_;
    };

    __always_inline Iterator end() const { return Iterator(this); }
    __always_inline Iterator begin() const {
        if (empty()) {
            return end();
        }

        Iterator    ret(this, 0, 0);
        const Node &first_node = nodes_.front();
        if (first_node.is_a_leaf() && first_node.leaves.size) return ret;

        return ++ret;
    }

    template <typename... Args>
    __always_inline std::pair<Iterator, bool> emplace(CoordinateType x,
                                                      CoordinateType y,
                                                      Args &&...args) {
        assert(is_inside_bounding_box(x, y, boundaries_));

        return emplace_recursively(0, boundaries_, x, y, std::forward<Args>(args)...);
    }

    std::vector<Iterator> nearest(CoordinateType x, CoordinateType y) const {
        CoordinateType        nearest_distance_squared = std::numeric_limits<CoordinateType>::max();
        std::vector<Iterator> results;
        nearest_recursively(0, boundaries_, x, y, nearest_distance_squared, results);

        return results;
    }

    __always_inline void find(const BoundingBox<CoordinateType> &bbox,
                              std::function<void(Iterator)>      func) const {
        find_recursively(bbox, boundaries_, func, 0);
    }

    __always_inline Iterator find(CoordinateType x, CoordinateType y) const {
        return find_recursively(x, y, boundaries_, 0);
    }

private:
    struct Quadrant {
        enum _Quadrant { NE = 0, NW = 1, SW = 2, SE = 3 };
    };

    struct NodeContentNonVoid {
        CoordinateType x;
        CoordinateType y;
        StorageType    storage;

        template <typename... Args>
        __always_inline NodeContentNonVoid(CoordinateType x_, CoordinateType y_, Args &&...args)
            : x(x_), y(y_), storage(std::forward<Args>(args)...) {}
    };

    struct NodeContentVoid {
        CoordinateType x;
        CoordinateType y;

        __always_inline NodeContentVoid(CoordinateType x_, CoordinateType y_) : x(x_), y(y_) {}
    };

    struct NodeContent
        : std::conditional<std::is_void_v<StorageType>, NodeContentVoid, NodeContentNonVoid>::type {
    };

    struct NodeLeaves {
        std::array<NodeContent, MAXIMUM_NODE_SIZE> items;
        uint8_t                                    size;
    };

    struct Node {
        static constexpr uint64_t NO_INDEX = std::numeric_limits<uint64_t>::max();
        union {
            NodeLeaves              leaves;
            std::array<uint64_t, 4> children;
        };
        bool branch;

        __always_inline Node() { reset(); }
        __always_inline Node(Node &&other) : branch(other.branch) {
            if (is_a_branch()) {
                children = std::move(other.children);
            } else {
                leaves = std::move(other.leaves);
            }
        }
        __always_inline ~Node() {
            if constexpr (!std::is_void_v<StorageType>) {
                if (is_a_leaf()) {
                    UNROLL_4
                    for (uint64_t i = 0; i < MAXIMUM_NODE_SIZE; ++i) {
                        leaves.items[i].storage.~StorageType();
                    }
                }
            }
        }
        __always_inline bool is_a_branch() const { return branch; }
        __always_inline bool is_a_leaf() const { return !is_a_branch(); }
        __always_inline void reset() {
            leaves.size = 0;
            branch = false;
        }
    };

    __always_inline auto operator()(uint64_t node_index, uint64_t item_index) const {
        assert(node_index < nodes_.size());
        const Node &node = nodes_[node_index];

        assert(node.is_a_leaf());
        assert(item_index < node.leaves.size);

        if constexpr (std::is_void_v<StorageType>) {
            const auto [x, y] = node.leaves.items[item_index];
            return std::tuple<CoordinateType, CoordinateType>(x, y);
        } else {
            const auto &[x, y, storage] = node.leaves.items[item_index];
            return std::tuple<CoordinateType, CoordinateType, const StorageType &>(x, y, storage);
        }
    }

    static __always_inline auto belongs_to_quadrant(const BoundingBox<CoordinateType> &boundaries,
                                                    CoordinateType                     x,
                                                    CoordinateType                     y) {
        CoordinateType origin_x = (boundaries.bottom_x - boundaries.top_x) / 2 + boundaries.top_x;
        CoordinateType origin_y =
            (boundaries.top_y - boundaries.bottom_y) / 2 + boundaries.bottom_y;

        // This logic must match integer flooring (less-or-equal instead of greater-or-equal).
        return Quadrant::NE + ((y > origin_y) * (x <= origin_x)) * Quadrant::NW +
               (y <= origin_y) * ((x <= origin_x) * Quadrant::SW +
                                  (y <= origin_y) * (x > origin_x) * Quadrant::SE);
    }

    static __always_inline BoundingBox<CoordinateType> compute_new_boundaries(
        uint64_t quad, const BoundingBox<CoordinateType> &boundaries) {
        static std::array<std::tuple<bool, bool>, 4> factors;
        factors[Quadrant::NE] = std::tuple<bool, bool>{1, 1};
        factors[Quadrant::NW] = std::tuple<bool, bool>{0, 1};
        factors[Quadrant::SW] = std::tuple<bool, bool>{0, 0};
        factors[Quadrant::SE] = std::tuple<bool, bool>{1, 0};
        const auto [fx, fy] = factors[quad];

        CoordinateType width = (boundaries.bottom_x - boundaries.top_x);
        CoordinateType height = (boundaries.top_y - boundaries.bottom_y);
        CoordinateType new_width = width / 2;
        CoordinateType new_height = height / 2;
        CoordinateType width_remainder = width - 2 * new_width;
        CoordinateType height_remainder = height - 2 * new_height;
        CoordinateType top_x = boundaries.top_x + fx * new_width;
        CoordinateType top_y = boundaries.top_y - (1 - fy) * (new_height + height_remainder);
        CoordinateType bottom_x = boundaries.bottom_x - (1 - fx) * (new_width + width_remainder);
        CoordinateType bottom_y = boundaries.bottom_y + fy * new_height;

        return {top_x, top_y, bottom_x, bottom_y};
    }

    template <typename... Args>
    __always_inline std::pair<Iterator, bool> emplace_recursively_helper(
        const std::array<uint64_t, 4>     &children,
        const BoundingBox<CoordinateType> &boundaries,
        CoordinateType                     x,
        CoordinateType                     y,
        Args &&...args) {
        const auto selected_quad = belongs_to_quadrant(boundaries, x, y);
        const auto new_boundaries = compute_new_boundaries(selected_quad, boundaries);

        return emplace_recursively(
            children[selected_quad], new_boundaries, x, y, std::forward<Args>(args)...);
    }

    template <typename... Args>
    std::pair<Iterator, bool> emplace_recursively(uint64_t                           node_index,
                                                  const BoundingBox<CoordinateType> &boundaries,
                                                  CoordinateType                     x,
                                                  CoordinateType                     y,
                                                  Args &&...args) {
        assert(node_index < nodes_.size());
        assert(is_inside_bounding_box(x, y, boundaries));

        Node &node = nodes_[node_index];

        if (node.is_a_branch()) {
            return emplace_recursively_helper(
                node.children, boundaries, x, y, std::forward<Args>(args)...);
        }

        int64_t item_index = 0;
        for (uint64_t i = 0; i < node.leaves.size; ++i) {
            auto x_ = node.leaves.items[i].x;
            auto y_ = node.leaves.items[i].y;
            item_index += (i + 1) * (x == x_ && y == y_);
        }

        if (item_index) {
            return {Iterator(this, node_index, item_index - 1), false};
        }

        const bool node_is_full = node.leaves.size == MAXIMUM_NODE_SIZE;
        if (!node_is_full) {
            node.leaves.items[node.leaves.size].x = x;
            node.leaves.items[node.leaves.size].y = y;
            if constexpr (!std::is_void_v<StorageType>) {
#if defined(CLANG_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmaybe-uninitialized"
#elif defined(GNU_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
                node.leaves.items[node.leaves.size].storage =
                    std::move(StorageType(std::forward<Args>(args)...));
#if defined(CLANG_COMPILER)
#pragma clang diagnostic pop
#elif defined(GNU_COMPILER)
#pragma GCC diagnostic pop
#endif
            }
            ++node.leaves.size;
            ++size_;
            return {Iterator(this, node_index, node.leaves.size - 1), true};
        }

        std::array<uint64_t, 4> new_children;
        std::iota(new_children.begin(), new_children.end(), nodes_.size());
        nodes_.resize(nodes_.size() + 4);

        Node &node_as_branch = nodes_[node_index];

        UNROLL_4
        for (uint64_t i = 0; i < MAXIMUM_NODE_SIZE; ++i) {
            auto x_ = node_as_branch.leaves.items[i].x;
            auto y_ = node_as_branch.leaves.items[i].y;
            if constexpr (std::is_void_v<StorageType>) {
                emplace_recursively_helper(new_children, boundaries, x_, y_);
            } else {
                emplace_recursively_helper(
                    new_children, boundaries, x_, y_, node_as_branch.leaves.items[i].storage);
            }
        }

        node_as_branch.branch = true;
        node_as_branch.children = new_children;
        size_ -= MAXIMUM_NODE_SIZE;

        return emplace_recursively_helper(
            new_children, boundaries, x, y, std::forward<Args>(args)...);
    }

    Iterator find_recursively(CoordinateType                     x,
                              CoordinateType                     y,
                              const BoundingBox<CoordinateType> &boundaries,
                              uint64_t                           node_index) const {
        assert(node_index < nodes_.size());

        const Node &node = nodes_[node_index];
        if (node.is_a_leaf()) {
            for (uint64_t i = 0; i < node.leaves.size; ++i) {
                if (x == node.leaves.items[i].x && y == node.leaves.items[i].y) {
                    return Iterator(this, node_index, i);
                }
            }
            return end();
        }

        const auto selected_quad = belongs_to_quadrant(boundaries, x, y);
        return find_recursively(
            x, y, compute_new_boundaries(selected_quad, boundaries), node.children[selected_quad]);
    }

    void find_recursively(const BoundingBox<CoordinateType> &bbox,
                          const BoundingBox<CoordinateType> &boundaries,
                          std::function<void(Iterator)>      func,
                          uint64_t                           node_index) const {
        assert(node_index < nodes_.size());

        const Node &node = nodes_[node_index];
        if (node.is_a_leaf()) {
            for (uint64_t i = 0; i < node.leaves.size; ++i) {
                if (is_inside_bounding_box(node.leaves.items[i].x, node.leaves.items[i].y, bbox)) {
                    func(Iterator(this, node_index, i));
                }
            }
            return;
        }

        UNROLL_4
        for (uint64_t i = 0; i < 4; ++i) {
            auto new_boundaries = compute_new_boundaries(i, boundaries);
            if (bounding_boxes_overlap(bbox, new_boundaries)) {
                find_recursively(bbox, new_boundaries, func, node.children[i]);
            }
        }
    }

    static __always_inline CoordinateType smallest_distance_from_bounding_box(
        const BoundingBox<CoordinateType> &bbox, CoordinateType x, CoordinateType y) {
        CoordinateType up = (y > bbox.top_y) * (y - bbox.top_y);
        CoordinateType down = (y < bbox.bottom_y) * (bbox.bottom_y - y);
        CoordinateType left = (x < bbox.top_x) * (bbox.top_x - x);
        CoordinateType right = (x > bbox.bottom_x) * (x - bbox.bottom_x);

        return up * up + down * down + left * left + right * right;
    }

    void nearest_recursively(uint64_t                           node_index,
                             const BoundingBox<CoordinateType> &boundaries,
                             CoordinateType                     x,
                             CoordinateType                     y,
                             CoordinateType                    &nearest_distance_squared,
                             std::vector<Iterator>             &results) const {
        assert(node_index < nodes_.size());

        const Node &node = nodes_[node_index];

        if (node.is_a_leaf()) {
            for (uint64_t i = 0; i < node.leaves.size; ++i) {
                auto distance = euclidean_distance_squared(
                    x, node.leaves.items[i].x, y, node.leaves.items[i].y);
                if (distance < nearest_distance_squared) {
                    results.clear();
                    nearest_distance_squared = distance;
                    results.emplace_back(this, node_index, i);
                } else if (distance == nearest_distance_squared) {
                    results.emplace_back(this, node_index, i);
                }
            }
            return;
        }

        const auto selected_quad = belongs_to_quadrant(boundaries, x, y);
        nearest_recursively(node.children[selected_quad],
                            compute_new_boundaries(selected_quad, boundaries),
                            x,
                            y,
                            nearest_distance_squared,
                            results);

        UNROLL_3
        for (uint64_t i = 1; i <= 3; ++i) {
            uint64_t   quad = (i + selected_quad) % 4;
            const auto new_boundaries = compute_new_boundaries(quad, boundaries);
            if (smallest_distance_from_bounding_box(new_boundaries, x, y) <=
                nearest_distance_squared) {
                nearest_recursively(
                    node.children[quad], new_boundaries, x, y, nearest_distance_squared, results);
            }
        }
    }

    const BoundingBox<CoordinateType> boundaries_;
    uint64_t                          size_;
    std::vector<Node>                 nodes_;
};
}  // namespace st