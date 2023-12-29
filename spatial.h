#pragma once

#if ((defined(__GNUC__) || defined(__GNUG__)) && !defined(__clang__))
#define __gnu_compiler
#elif (defined(__clang__))
#define __clang_compiler
#endif

#if defined(__gnu_compiler)
#define __unroll_2 _Pragma("GCC unroll 2")
#define __unroll_3 _Pragma("GCC unroll 3")
#define __unroll_4 _Pragma("GCC unroll 4")
#define __unroll_8 _Pragma("GCC unroll 8")
#elif defined(__clang_compiler)
#define __unroll_2 _Pragma("clang loop unroll_count(2)")
#define __unroll_3 _Pragma("clang loop unroll_count(3)")
#define __unroll_4 _Pragma("clang loop unroll_count(4)")
#define __unroll_8 _Pragma("clang loop unroll_count(8)")
#endif

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <tuple>
#include <vector>

namespace st {

#ifndef __always_inline
#define __undef_always_inline
#define __always_inline __inline __attribute__((__always_inline__))
#endif

template <typename CoordinateType>
struct bounding_box {
    CoordinateType  top_x, top_y, bottom_x, bottom_y;
    __always_inline bounding_box()
        : top_x(std::numeric_limits<CoordinateType>::lowest() / 2),
          top_y(std::numeric_limits<CoordinateType>::max() / 2),
          bottom_x(std::numeric_limits<CoordinateType>::max() / 2),
          bottom_y(std::numeric_limits<CoordinateType>::lowest() / 2) {
        assert(top_x <= bottom_x);
        assert(top_y >= bottom_y);
    }
    __always_inline bounding_box(CoordinateType top_x_,
                                 CoordinateType top_y_,
                                 CoordinateType bottom_x_,
                                 CoordinateType bottom_y_)
        : top_x(top_x_), top_y(top_y_), bottom_x(bottom_x_), bottom_y(bottom_y_) {
        assert(top_x <= bottom_x);
        assert(top_y >= bottom_y);
    }
    __always_inline CoordinateType area() const { return (bottom_x - top_x) * (top_y - bottom_y); }
};

namespace internal {

template <typename AbsDiff, typename T, typename... Args>
__always_inline T euclidean_distance_squared_impl(AbsDiff) {
    return T(0);
}

template <typename AbsDiff, typename T, typename... Args>
__always_inline T euclidean_distance_squared_impl(AbsDiff absdiff, T x1, T x2, Args... xs) {
    const T sum = euclidean_distance_squared_impl<AbsDiff, T>(absdiff, xs...);
    const T dx = absdiff(x2, x1);

    return sum + dx * dx;
}

template <typename T>
struct SafeAbsDiff {
    static __always_inline T safe_absdiff(T x, T y) {
        const bool gt = x > y;
        return gt * (x - y) + (1 - gt) * (y - x);
    }
    __always_inline T operator()(T x, T y) { return safe_absdiff(x, y); }
};

template <typename T>
struct UnsafeAbsDiff {
    static __always_inline T unsafe_absdiff(T x, T y) { return x - y; }
    __always_inline T        operator()(T x, T y) { return unsafe_absdiff(x, y); }
};

template <typename T>
T absdiff(T x, T y) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::fabs(UnsafeAbsDiff<T>()(x, y));
    } else {
        return SafeAbsDiff<T>()(x, y);
    }
}

template <typename T, typename... Args>
__always_inline T euclidean_distance_squared(T x1, T x2, Args... xs) {
    if constexpr (std::is_floating_point_v<T> || std::is_signed_v<T>) {
        return euclidean_distance_squared_impl(UnsafeAbsDiff<T>(), x1, x2, xs...);
    } else {
        return euclidean_distance_squared_impl(SafeAbsDiff<T>(), x1, x2, xs...);
    }
}

template <typename CoordinateType>
__always_inline bool is_within_interval(CoordinateType x, CoordinateType beg, CoordinateType end) {
    assert(end >= beg);
    return x >= beg && x <= end;
}

template <typename CoordinateType>
static __always_inline bool is_inside_bounding_box(CoordinateType                      x,
                                                   CoordinateType                      y,
                                                   const bounding_box<CoordinateType> &bbox) {
    return is_within_interval(x, bbox.top_x, bbox.bottom_x) &&
           is_within_interval(y, bbox.bottom_y, bbox.top_y);
}

template <typename CoordinateType>
static __always_inline bool intervals_overlap(CoordinateType lhs_beg,
                                              CoordinateType lhs_end,
                                              CoordinateType rhs_beg,
                                              CoordinateType rhs_end) {
    assert(lhs_beg <= lhs_end);
    assert(rhs_beg <= rhs_end);

    return !((lhs_beg > rhs_end) || (lhs_end < rhs_beg));
}

template <typename CoordinateType>
static __always_inline bool bounding_boxes_overlap(const bounding_box<CoordinateType> &lhs,
                                                   const bounding_box<CoordinateType> &rhs) {
    return intervals_overlap(lhs.top_x, lhs.bottom_x, rhs.top_x, rhs.bottom_x) &&
           intervals_overlap(lhs.bottom_y, lhs.top_y, rhs.bottom_y, rhs.top_y);
}

template <typename StorageType = void,
          typename CoordinateType = int,
          uint8_t MAXIMUM_NODE_SIZE = 32>
class spatial_tree {
public:
    static_assert(MAXIMUM_NODE_SIZE > 0, "Maximum node size must be greater than 1");

    spatial_tree()
        // Dividing by 2 to be able to do (upper bound +/- lower bound) without overflowing.
        : boundaries_(std::numeric_limits<CoordinateType>::lowest() / 2,
                      std::numeric_limits<CoordinateType>::max() / 2,
                      std::numeric_limits<CoordinateType>::max() / 2,
                      std::numeric_limits<CoordinateType>::lowest() / 2) {
        clear();
    }
    spatial_tree(const bounding_box<CoordinateType> &boundaries) : boundaries_(boundaries) {
        assert(boundaries_.top_x >= std::numeric_limits<CoordinateType>::lowest() / 2);
        assert(boundaries_.top_y <= std::numeric_limits<CoordinateType>::max() / 2);
        assert(boundaries_.bottom_x <= std::numeric_limits<CoordinateType>::max() / 2);
        assert(boundaries_.bottom_y >= std::numeric_limits<CoordinateType>::lowest() / 2);
        clear();
    }

    ~spatial_tree() = default;

    __always_inline uint64_t size() const { return size_; }
    __always_inline bool     empty() const { return size() == 0; }
    __always_inline void     clear() {
        nodes_.resize(1);
        nodes_.front().reset();
        size_ = 0;
    }

    struct iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = iterator;
        using pointer = value_type *;
        using reference = value_type &;

        iterator(const spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> *tree,
                 uint64_t                                                            node_index,
                 uint8_t                                                             item_index)
            : tree_(tree), node_index_(node_index), item_index_(item_index) {}
        iterator(const spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> *tree)
            : tree_(tree) {
            end();
        }
        ~iterator() = default;

        __always_inline auto operator*() const {
            return tree_->operator()(node_index_, item_index_);
        }

        iterator &operator++() {
            assert(node_index_ != tree_node::NO_INDEX);
            assert(node_index_ < tree_->nodes_.size());

            const tree_node *node = &tree_->nodes_[node_index_];
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

        __always_inline bool operator==(const iterator &other) const {
            assert(tree_ == other.tree_);
            return node_index_ == other.node_index_ && item_index_ == other.item_index_;
        }

        __always_inline bool operator!=(const iterator &other) const { return !(*this == other); }

    private:
        __always_inline void end() {
            node_index_ = tree_node::NO_INDEX;
            item_index_ = 0;
        }

        const spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> *tree_;

        uint64_t node_index_;
        uint8_t  item_index_;
    };

    __always_inline iterator end() const { return iterator(this); }
    __always_inline iterator begin() const {
        if (empty()) {
            return end();
        }

        iterator         ret(this, 0, 0);
        const tree_node &first_node = nodes_.front();
        if (first_node.is_a_leaf() && first_node.leaves.size) return ret;

        return ++ret;
    }

    template <typename... Args>
    __always_inline std::pair<iterator, bool> emplace(CoordinateType x,
                                                      CoordinateType y,
                                                      Args &&...args) {
        assert(is_inside_bounding_box(x, y, boundaries_));

        return emplace_recursively(0, boundaries_, x, y, std::forward<Args>(args)...);
    }

    std::vector<iterator> nearest(CoordinateType x, CoordinateType y) const {
        CoordinateType        nearest_distance_squared = std::numeric_limits<CoordinateType>::max();
        std::vector<iterator> results;
        nearest_recursively(0, boundaries_, x, y, nearest_distance_squared, results);

        return results;
    }

    __always_inline void find(const bounding_box<CoordinateType> &bbox,
                              std::function<void(iterator)>       func) const {
        find_recursively(bbox, boundaries_, func, 0);
    }

    __always_inline iterator find(CoordinateType x, CoordinateType y) const {
        return find_recursively(x, y, boundaries_, 0);
    }

private:
    // Acts as a namespace and it is still possible to index without an explicit cast.
    struct cartesian_quadrant {
        enum _cartesian_quadrant { NE = 0, NW = 1, SW = 2, SE = 3 };
    };

    struct tree_node_with_storage {
        CoordinateType x;
        CoordinateType y;
        StorageType    storage;

        template <typename... Args>
        __always_inline tree_node_with_storage(CoordinateType x_, CoordinateType y_, Args &&...args)
            : x(x_), y(y_), storage(std::forward<Args>(args)...) {}
    };

    struct tree_node_no_storage {
        CoordinateType x;
        CoordinateType y;

        __always_inline tree_node_no_storage(CoordinateType x_, CoordinateType y_) : x(x_), y(y_) {}
    };

    struct tree_nodeContent : std::conditional<std::is_void_v<StorageType>,
                                               tree_node_no_storage,
                                               tree_node_with_storage>::type {};

    struct tree_nodeLeaves {
        std::array<tree_nodeContent, MAXIMUM_NODE_SIZE> items;
        uint8_t                                         size;
    };

    struct tree_node {
        static constexpr uint64_t NO_INDEX = std::numeric_limits<uint64_t>::max();
        union {
            tree_nodeLeaves         leaves;
            std::array<uint64_t, 4> children;
        };
        bool branch;

        __always_inline tree_node() { reset(); }
        __always_inline tree_node(tree_node &&other) : branch(other.branch) {
            if (is_a_branch()) {
                children = std::move(other.children);
            } else {
                leaves = std::move(other.leaves);
            }
        }
        __always_inline ~tree_node() {
            if constexpr (!std::is_void_v<StorageType>) {
                if (is_a_leaf()) {
                    __unroll_4 for (uint64_t i = 0; i < MAXIMUM_NODE_SIZE; ++i) {
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
        const tree_node &node = nodes_[node_index];

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

    static __always_inline auto belongs_to_quadrant(const bounding_box<CoordinateType> &boundaries,
                                                    CoordinateType                      x,
                                                    CoordinateType                      y) {
        CoordinateType origin_x = (boundaries.bottom_x - boundaries.top_x) / 2 + boundaries.top_x;
        CoordinateType origin_y =
            (boundaries.top_y - boundaries.bottom_y) / 2 + boundaries.bottom_y;

        // This logic must match integer flooring (less-or-equal instead of greater-or-equal).
        return cartesian_quadrant::NE +
               ((y > origin_y) * (x <= origin_x)) * cartesian_quadrant::NW +
               (y <= origin_y) * ((x <= origin_x) * cartesian_quadrant::SW +
                                  (y <= origin_y) * (x > origin_x) * cartesian_quadrant::SE);
    }

    static __always_inline bounding_box<CoordinateType> compute_new_boundaries(
        uint64_t quad, const bounding_box<CoordinateType> &boundaries) {
        static std::array<std::tuple<bool, bool>, 4> factors;
        factors[cartesian_quadrant::NE] = std::tuple<bool, bool>{1, 1};
        factors[cartesian_quadrant::NW] = std::tuple<bool, bool>{0, 1};
        factors[cartesian_quadrant::SW] = std::tuple<bool, bool>{0, 0};
        factors[cartesian_quadrant::SE] = std::tuple<bool, bool>{1, 0};
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
    __always_inline std::pair<iterator, bool> emplace_recursively_helper(
        const std::array<uint64_t, 4>      &children,
        const bounding_box<CoordinateType> &boundaries,
        CoordinateType                      x,
        CoordinateType                      y,
        Args &&...args) {
        const auto selected_quad = belongs_to_quadrant(boundaries, x, y);
        const auto new_boundaries = compute_new_boundaries(selected_quad, boundaries);

        return emplace_recursively(
            children[selected_quad], new_boundaries, x, y, std::forward<Args>(args)...);
    }

    template <typename... Args>
    std::pair<iterator, bool> emplace_recursively(uint64_t                            node_index,
                                                  const bounding_box<CoordinateType> &boundaries,
                                                  CoordinateType                      x,
                                                  CoordinateType                      y,
                                                  Args &&...args) {
        assert(node_index < nodes_.size());
        assert(is_inside_bounding_box(x, y, boundaries));

        tree_node &node = nodes_[node_index];

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
            return {iterator(this, node_index, item_index - 1), false};
        }

        const bool node_is_full = node.leaves.size == MAXIMUM_NODE_SIZE;
        if (!node_is_full) {
            node.leaves.items[node.leaves.size].x = x;
            node.leaves.items[node.leaves.size].y = y;
            if constexpr (!std::is_void_v<StorageType>) {
#if defined(__clang_compiler)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmaybe-uninitialized"
#elif defined(__gnu_compiler)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
                node.leaves.items[node.leaves.size].storage =
                    std::move(StorageType(std::forward<Args>(args)...));
#if defined(__clang_compiler)
#pragma clang diagnostic pop
#elif defined(__gnu_compiler)
#pragma GCC diagnostic pop
#endif
            }
            ++node.leaves.size;
            ++size_;
            return {iterator(this, node_index, node.leaves.size - 1), true};
        }

        std::array<uint64_t, 4> new_children;
        std::iota(new_children.begin(), new_children.end(), nodes_.size());
        nodes_.resize(nodes_.size() + 4);

        tree_node &node_as_branch = nodes_[node_index];

        __unroll_4 for (uint64_t i = 0; i < MAXIMUM_NODE_SIZE; ++i) {
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

    iterator find_recursively(CoordinateType                      x,
                              CoordinateType                      y,
                              const bounding_box<CoordinateType> &boundaries,
                              uint64_t                            node_index) const {
        assert(node_index < nodes_.size());

        const tree_node &node = nodes_[node_index];
        if (node.is_a_leaf()) {
            for (uint64_t i = 0; i < node.leaves.size; ++i) {
                if (x == node.leaves.items[i].x && y == node.leaves.items[i].y) {
                    return iterator(this, node_index, i);
                }
            }
            return end();
        }

        const auto selected_quad = belongs_to_quadrant(boundaries, x, y);
        return find_recursively(
            x, y, compute_new_boundaries(selected_quad, boundaries), node.children[selected_quad]);
    }

    void find_recursively(const bounding_box<CoordinateType> &bbox,
                          const bounding_box<CoordinateType> &boundaries,
                          std::function<void(iterator)>       func,
                          uint64_t                            node_index) const {
        assert(node_index < nodes_.size());

        const tree_node &node = nodes_[node_index];
        if (node.is_a_leaf()) {
            for (uint64_t i = 0; i < node.leaves.size; ++i) {
                if (is_inside_bounding_box(node.leaves.items[i].x, node.leaves.items[i].y, bbox)) {
                    func(iterator(this, node_index, i));
                }
            }
            return;
        }

        __unroll_4 for (uint64_t i = 0; i < 4; ++i) {
            auto new_boundaries = compute_new_boundaries(i, boundaries);
            if (bounding_boxes_overlap(bbox, new_boundaries)) {
                find_recursively(bbox, new_boundaries, func, node.children[i]);
            }
        }
    }

    static __always_inline CoordinateType smallest_distance_from_bounding_box(
        const bounding_box<CoordinateType> &bbox, CoordinateType x, CoordinateType y) {
        CoordinateType up = (y > bbox.top_y) * (y - bbox.top_y);
        CoordinateType down = (y < bbox.bottom_y) * (bbox.bottom_y - y);
        CoordinateType left = (x < bbox.top_x) * (bbox.top_x - x);
        CoordinateType right = (x > bbox.bottom_x) * (x - bbox.bottom_x);

        return up * up + down * down + left * left + right * right;
    }

    void nearest_recursively(uint64_t                            node_index,
                             const bounding_box<CoordinateType> &boundaries,
                             CoordinateType                      x,
                             CoordinateType                      y,
                             CoordinateType                     &nearest_distance_squared,
                             std::vector<iterator>              &results) const {
        assert(node_index < nodes_.size());

        const tree_node &node = nodes_[node_index];

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

        __unroll_3 for (uint64_t i = 1; i <= 3; ++i) {
            uint64_t   quad = (i + selected_quad) % 4;
            const auto new_boundaries = compute_new_boundaries(quad, boundaries);
            if (smallest_distance_from_bounding_box(new_boundaries, x, y) <=
                nearest_distance_squared) {
                nearest_recursively(
                    node.children[quad], new_boundaries, x, y, nearest_distance_squared, results);
            }
        }
    }

    const bounding_box<CoordinateType> boundaries_;
    uint64_t                           size_;
    std::vector<tree_node>             nodes_;
};
}  // namespace internal

template <typename CoordinateType, typename StorageType, uint8_t MAXIMUM_NODE_SIZE = 32>
class spatial_map : public internal::spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> {
    static_assert(!std::is_void_v<StorageType>, "For no storage type, use st::spatial_set");
};

template <typename CoordinateType = double, uint8_t MAXIMUM_NODE_SIZE = 32>
class spatial_set : public internal::spatial_tree<void, CoordinateType, MAXIMUM_NODE_SIZE> {};
}  // namespace st

#ifdef __gnu_compiler
#undef __gnu_compiler
#endif

#ifdef __clang_compiler
#undef __clang_compiler
#endif

#ifdef __unroll_2
#undef __unroll_2
#undef __unroll_3
#undef __unroll_4
#undef __unroll_8
#endif

#ifdef __undef_always_inline
#undef __always_inline
#endif