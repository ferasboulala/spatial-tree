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

template <typename CoordinateType>
struct bounding_box {
    CoordinateType top_x, top_y, bottom_x, bottom_y;
    inline bounding_box()
        : top_x(std::numeric_limits<CoordinateType>::lowest() / 2),
          top_y(std::numeric_limits<CoordinateType>::max() / 2),
          bottom_x(std::numeric_limits<CoordinateType>::max() / 2),
          bottom_y(std::numeric_limits<CoordinateType>::lowest() / 2) {
        assert(top_x <= bottom_x);
        assert(top_y >= bottom_y);
    }
    inline bounding_box(CoordinateType top_x_,
                        CoordinateType top_y_,
                        CoordinateType bottom_x_,
                        CoordinateType bottom_y_)
        : top_x(top_x_), top_y(top_y_), bottom_x(bottom_x_), bottom_y(bottom_y_) {
        assert(top_x <= bottom_x);
        assert(top_y >= bottom_y);
    }
    inline CoordinateType area() const { return (bottom_x - top_x) * (top_y - bottom_y); }
};

namespace internal {

template <typename AbsDiff, typename T, typename... Args>
static inline T euclidean_distance_squared_impl(AbsDiff) {
    return T(0);
}

template <typename AbsDiff, typename T, typename... Args>
static inline T euclidean_distance_squared_impl(AbsDiff absdiff, T x1, T x2, Args... xs) {
    const T sum = euclidean_distance_squared_impl<AbsDiff, T>(absdiff, xs...);
    const T dx = absdiff(x2, x1);

    return sum + dx * dx;
}

template <typename T>
struct SafeAbsDiff {
    static inline T safe_absdiff(T x, T y) {
        const bool gt = x > y;
        return gt * (x - y) + (1 - gt) * (y - x);
    }
    inline T operator()(T x, T y) { return safe_absdiff(x, y); }
};

template <typename T>
struct UnsafeAbsDiff {
    static inline T unsafe_absdiff(T x, T y) { return x - y; }
    inline T        operator()(T x, T y) { return unsafe_absdiff(x, y); }
};

template <typename T>
static inline T absdiff(T x, T y) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::fabs(UnsafeAbsDiff<T>()(x, y));
    } else {
        return SafeAbsDiff<T>()(x, y);
    }
}

template <typename T, typename... Args>
static inline T euclidean_distance_squared(T x1, T x2, Args... xs) {
    if constexpr (std::is_floating_point_v<T> || std::is_signed_v<T>) {
        return euclidean_distance_squared_impl(UnsafeAbsDiff<T>(), x1, x2, xs...);
    } else {
        return euclidean_distance_squared_impl(SafeAbsDiff<T>(), x1, x2, xs...);
    }
}

template <typename CoordinateType>
static inline bool is_within_interval(CoordinateType x, CoordinateType beg, CoordinateType end) {
    assert(end >= beg);
    return x >= beg && x <= end;
}

template <typename CoordinateType>
static inline bool is_inside_bounding_box(CoordinateType                      x,
                                          CoordinateType                      y,
                                          const bounding_box<CoordinateType> &bbox) {
    return is_within_interval(x, bbox.top_x, bbox.bottom_x) &&
           is_within_interval(y, bbox.bottom_y, bbox.top_y);
}

template <typename CoordinateType>
static inline bool intervals_overlap(CoordinateType lhs_beg,
                                     CoordinateType lhs_end,
                                     CoordinateType rhs_beg,
                                     CoordinateType rhs_end) {
    assert(lhs_beg <= lhs_end);
    assert(rhs_beg <= rhs_end);

    return !((lhs_beg > rhs_end) || (lhs_end < rhs_beg));
}

template <typename CoordinateType>
static inline bool bounding_boxes_overlap(const bounding_box<CoordinateType> &lhs,
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

    inline uint64_t size() const {
        assert(nodes_.size());
        const auto &root = nodes_.front();
        if (root.is_a_leaf()) {
            return root.leaf.size;
        }
        return root.branch.size;
    }
    inline bool empty() const { return size() == 0; }
    inline void clear() {
        nodes_.resize(1);
        nodes_.front().reset();
        freed_nodes_.clear();
    }

    template <typename Func>
    void walk(Func func) const {
        const std::function<void(uint64_t, const bounding_box<CoordinateType> &)> walk_recursively =
            [&](auto node_index, auto boundaries) {
                assert(node_index < nodes_.size());
                const tree_node &node = nodes_[node_index];
                func(boundaries, node.is_a_leaf());

                if (node.is_a_leaf()) {
                    return;
                }

                for (uint64_t quad = 0; quad < 4; ++quad) {
                    const uint64_t child_index = node.branch.index_of_first_child + quad;
                    const auto     new_boundaries = compute_new_boundaries(quad, boundaries);
                    walk_recursively(child_index, new_boundaries);
                }
            };

        walk_recursively(0, boundaries_);
    }

    struct iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = iterator;
        using pointer = value_type *;
        using reference = value_type &;

        inline iterator(const spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> *tree,
                        uint64_t node_index,
                        uint8_t  item_index)
            : tree_(tree), node_index_(node_index), item_index_(item_index) {}
        inline iterator(const spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> *tree)
            : tree_(tree) {
            end();
        }
        inline ~iterator() = default;

        inline auto operator*() const { return tree_->operator()(node_index_, item_index_); }
        inline auto operator*() {
            return const_cast<spatial_tree *>(tree_)->operator()(node_index_, item_index_);
        }
        inline auto operator->() const { return this->operator(); }
        inline auto operator->() { return this->operator(); }

        iterator &operator++() {
            assert(node_index_ != tree_node::NO_INDEX);
            assert(node_index_ < tree_->nodes_.size());

            const tree_node *node = &tree_->nodes_[node_index_];
            // TODO: Make the is_leaf member variable an enum for the freed case to skip 4 leaves at
            // a time.
            if (node->is_a_leaf() && ++item_index_ < node->leaf.size) {
                return *this;
            }

            do {
                ++node_index_;
                node = &tree_->nodes_[node_index_];
                item_index_ = 0;
            } while ((node->is_a_branch() || !node->leaf.size) &&
                     node_index_ < tree_->nodes_.size());

            if (node_index_ == tree_->nodes_.size()) end();

            return *this;
        }

        inline bool operator==(const iterator &other) const {
            assert(tree_ == other.tree_);
            return node_index_ == other.node_index_ && item_index_ == other.item_index_;
        }

        inline bool operator!=(const iterator &other) const { return !(*this == other); }

    private:
        inline void end() {
            node_index_ = tree_node::NO_INDEX;
            item_index_ = 0;
        }

        const spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> *tree_;

        uint64_t node_index_;
        uint8_t  item_index_;
    };

    inline iterator end() const { return iterator(this); }
    inline iterator begin() const {
        if (empty()) {
            return end();
        }

        iterator         ret(this, 0, 0);
        const tree_node &first_node = nodes_.front();
        if (first_node.is_a_leaf() && first_node.leaf.size) return ret;

        return ++ret;
    }

    template <typename... Args>
    inline std::pair<iterator, bool> emplace(CoordinateType x, CoordinateType y, Args &&...args) {
        assert(is_inside_bounding_box(x, y, boundaries_));

        return emplace_recursively(0, boundaries_, x, y, std::forward<Args>(args)...);
    }

    inline bool erase(CoordinateType x, CoordinateType y) {
        return erase_recursively(0, boundaries_, size(), x, y);
    }

    std::vector<iterator> nearest(CoordinateType x, CoordinateType y) const {
        CoordinateType        nearest_distance_squared = std::numeric_limits<CoordinateType>::max();
        std::vector<iterator> results;
        nearest_recursively(0, boundaries_, x, y, nearest_distance_squared, results);

        return results;
    }

    template <typename Func>
    inline void find(const bounding_box<CoordinateType> &bbox, Func func) const {
        find_recursively(bbox, boundaries_, func, 0);
    }

    inline iterator find(CoordinateType x, CoordinateType y) const {
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
        inline tree_node_with_storage(CoordinateType x_, CoordinateType y_, Args &&...args)
            : x(x_), y(y_), storage(std::forward<Args>(args)...) {}
    };

    struct tree_node_no_storage {
        CoordinateType x;
        CoordinateType y;

        inline tree_node_no_storage(CoordinateType x_, CoordinateType y_) : x(x_), y(y_) {}
    };

    struct tree_node_storage : std::conditional<std::is_void_v<StorageType>,
                                                tree_node_no_storage,
                                                tree_node_with_storage>::type {};

    struct tree_node_leaf {
        std::array<tree_node_storage, MAXIMUM_NODE_SIZE> items;
        uint8_t                                          size;
    };

    struct tree_node_branch {
        // Children are adjacent in memory.
        uint64_t index_of_first_child;
        uint64_t size;
    };

    struct tree_node {
        static constexpr uint64_t NO_INDEX = std::numeric_limits<uint64_t>::max();
        union {
            tree_node_leaf   leaf;
            tree_node_branch branch;
        };
        bool is_branch;

        inline tree_node() { reset(); }
        inline tree_node(tree_node &&other) : is_branch(other.is_branch) {
            if (is_a_branch()) {
                branch = std::move(other.branch);
            } else {
                leaf = std::move(other.leaf);
            }
        }
        inline ~tree_node() {
            if constexpr (!std::is_void_v<StorageType>) {
                if (is_a_leaf()) {
                    __unroll_4 for (uint64_t i = 0; i < MAXIMUM_NODE_SIZE; ++i) {
                        leaf.items[i].storage.~StorageType();
                    }
                }
            }
        }
        inline bool is_a_branch() const { return is_branch; }
        inline bool is_a_leaf() const { return !is_a_branch(); }
        inline void reset() {
            leaf.size = 0;
            is_branch = false;
        }
    };

    inline auto operator()(uint64_t node_index, uint64_t item_index) const {
        assert(node_index < nodes_.size());
        const tree_node &node = nodes_[node_index];

        assert(node.is_a_leaf());
        assert(item_index < node.leaf.size);

        if constexpr (std::is_void_v<StorageType>) {
            const auto [x, y] = node.leaf.items[item_index];
            return std::tuple<CoordinateType, CoordinateType>(x, y);
        } else {
            const auto &[x, y, storage] = node.leaf.items[item_index];
            return std::tuple<CoordinateType, CoordinateType, const StorageType &>(x, y, storage);
        }
    }

    inline auto operator()(uint64_t node_index, uint64_t item_index) {
        auto this_as_const = const_cast<const spatial_tree *>(this);
        if constexpr (std::is_void_v<StorageType>) {
            return this_as_const->operator()(node_index, item_index);
        } else {
            const auto &[x, y, storage] = this_as_const->operator()(node_index, item_index);
            return std::tuple<CoordinateType, CoordinateType, StorageType &>(
                x, y, const_cast<StorageType &>(storage));
        }
    }

    static inline auto belongs_to_quadrant(const bounding_box<CoordinateType> &boundaries,
                                           CoordinateType                      x,
                                           CoordinateType                      y) {
        CoordinateType origin_x = (boundaries.bottom_x - boundaries.top_x) / 2 + boundaries.top_x;
        CoordinateType origin_y =
            (boundaries.top_y - boundaries.bottom_y) / 2 + boundaries.bottom_y;

        // Faster than branchless.
        if (x <= origin_x) {
            if (y > origin_y) {
                return cartesian_quadrant::NW;
            }
            return cartesian_quadrant::SW;
        } else {
            if (y <= origin_y) {
                return cartesian_quadrant::SE;
            }
            return cartesian_quadrant::NE;
        }
    }

    static inline bounding_box<CoordinateType> compute_new_boundaries(
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

        CoordinateType top_x = boundaries.top_x;
        CoordinateType bottom_x = boundaries.bottom_x;
        CoordinateType top_y = boundaries.top_y;
        CoordinateType bottom_y = boundaries.bottom_y;

        // Faster than branchless.
        if (fx) {
            top_x += new_width;
        } else {
            bottom_x -= (new_width + width_remainder);
        }

        if (fy) {
            bottom_y += new_height;
        } else {
            top_y -= (new_height + height_remainder);
        }

        return {top_x, top_y, bottom_x, bottom_y};
    }

    template <typename... Args>
    inline std::pair<iterator, bool> emplace_recursively_helper(
        uint64_t                            index_of_first_child,
        const bounding_box<CoordinateType> &boundaries,
        CoordinateType                      x,
        CoordinateType                      y,
        Args &&...args) {
        const auto selected_quad = belongs_to_quadrant(boundaries, x, y);
        const auto new_boundaries = compute_new_boundaries(selected_quad, boundaries);

        return emplace_recursively(index_of_first_child + selected_quad,
                                   new_boundaries,
                                   x,
                                   y,
                                   std::forward<Args>(args)...);
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
            const auto ret = emplace_recursively_helper(
                node.branch.index_of_first_child, boundaries, x, y, std::forward<Args>(args)...);
            nodes_[node_index].branch.size += ret.second;

            return ret;
        }

        int64_t item_index = 0;
        for (uint64_t i = 0; i < node.leaf.size; ++i) {
            auto x_ = node.leaf.items[i].x;
            auto y_ = node.leaf.items[i].y;
            item_index += (i + 1) * (x == x_ && y == y_);
        }

        if (item_index) {
            return {iterator(this, node_index, item_index - 1), false};
        }

        const bool node_is_full = node.leaf.size == MAXIMUM_NODE_SIZE;
        if (!node_is_full) {
            node.leaf.items[node.leaf.size].x = x;
            node.leaf.items[node.leaf.size].y = y;
            if constexpr (!std::is_void_v<StorageType>) {
#if defined(__gnu_compiler)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
                node.leaf.items[node.leaf.size].storage =
                    std::move(StorageType(std::forward<Args>(args)...));
#if defined(__gnu_compiler)
#pragma GCC diagnostic pop
#endif
            }
            ++node.leaf.size;
            return {iterator(this, node_index, node.leaf.size - 1), true};
        }

        uint64_t new_index_of_first_child;
        if (freed_nodes_.size()) {
            new_index_of_first_child = freed_nodes_.back();
            freed_nodes_.pop_back();
        } else {
            new_index_of_first_child = nodes_.size();
            nodes_.resize(nodes_.size() + 4);
        }

        tree_node &node_as_branch = nodes_[node_index];
        __unroll_4 for (uint64_t i = 0; i < MAXIMUM_NODE_SIZE; ++i) {
            auto x_ = node_as_branch.leaf.items[i].x;
            auto y_ = node_as_branch.leaf.items[i].y;
            if constexpr (std::is_void_v<StorageType>) {
                emplace_recursively_helper(new_index_of_first_child, boundaries, x_, y_);
            } else {
                emplace_recursively_helper(new_index_of_first_child,
                                           boundaries,
                                           x_,
                                           y_,
                                           node_as_branch.leaf.items[i].storage);
            }
        }

        node_as_branch.is_branch = true;
        node_as_branch.branch.index_of_first_child = new_index_of_first_child;
        node_as_branch.branch.size = MAXIMUM_NODE_SIZE + 1;

        const auto ret = emplace_recursively_helper(
            new_index_of_first_child, boundaries, x, y, std::forward<Args>(args)...);
        assert(ret.second);

        return ret;
    }

    inline void recursively_gather_points(uint64_t node_index, tree_node_leaf &leaf) {
        assert(node_index < nodes_.size());

        tree_node &node = nodes_[node_index];
        if (node.is_a_leaf()) {
            for (uint8_t i = 0; i < node.leaf.size; ++i) {
                leaf.items[leaf.size++] = std::move(node.leaf.items[i]);
            }
            assert(leaf.size <= MAXIMUM_NODE_SIZE);
            node.reset();
            return;
        }

        // TODO: Specify that it is safe to access &node because there was no realloc.
        // TODO: Do it elsewhere too.
        for (uint64_t quad = 0; quad < 4; ++quad) {
            recursively_gather_points(node.branch.index_of_first_child + quad, leaf);
        }

        freed_nodes_.push_back(node.branch.index_of_first_child);
        node.reset();
    }

    inline bool erase_recursively(uint64_t                            node_index,
                                  const bounding_box<CoordinateType> &boundaries,
                                  uint64_t                            parent_size_after_removal,
                                  CoordinateType                      x,
                                  CoordinateType                      y) {
        assert(node_index < nodes_.size());

        tree_node &node = nodes_[node_index];
        if (node.is_a_leaf()) {
            for (uint8_t i = 0; i < node.leaf.size; ++i) {
                if (x == node.leaf.items[i].x && y == node.leaf.items[i].y) {
                    uint64_t last_index = node.leaf.size - 1;
                    if (last_index != i) {
                        node.leaf.items[i] = std::move(node.leaf.items[last_index]);
                    }
                    --node.leaf.size;

                    return true;
                }
            }

            return false;
        }

        uint64_t   size_after_removal = node.branch.size - 1;
        const auto selected_quad = belongs_to_quadrant(boundaries, x, y);
        const bool erased = erase_recursively(node.branch.index_of_first_child + selected_quad,
                                              compute_new_boundaries(selected_quad, boundaries),
                                              size_after_removal,
                                              x,
                                              y);

        if (!erased) {
            return false;
        }

        if (parent_size_after_removal <= MAXIMUM_NODE_SIZE) {
            return true;
        }

        if (--node.branch.size > MAXIMUM_NODE_SIZE) {
            return true;
        }

        // TODO: Specify that it is safe to access &node because there was no realloc.
        /// EXPLANATION: Not doing it at the root because of the union.
        uint64_t index_of_first_child = node.branch.index_of_first_child;
        node.reset();
        for (uint64_t quad = 0; quad < 4; ++quad) {
            recursively_gather_points(index_of_first_child + quad, node.leaf);
        }
        freed_nodes_.push_back(index_of_first_child);
        assert(node.leaf.size <= MAXIMUM_NODE_SIZE);

        return true;
    }

    inline iterator find_recursively(CoordinateType                      x,
                                     CoordinateType                      y,
                                     const bounding_box<CoordinateType> &boundaries,
                                     uint64_t                            node_index) const {
        assert(node_index < nodes_.size());

        const tree_node &node = nodes_[node_index];
        if (node.is_a_leaf()) {
            for (uint64_t i = 0; i < node.leaf.size; ++i) {
                if (x == node.leaf.items[i].x && y == node.leaf.items[i].y) {
                    return iterator(this, node_index, i);
                }
            }
            return end();
        }

        const auto selected_quad = belongs_to_quadrant(boundaries, x, y);
        return find_recursively(x,
                                y,
                                compute_new_boundaries(selected_quad, boundaries),
                                node.branch.index_of_first_child + selected_quad);
    }

    template <typename Func>
    inline void find_recursively(const bounding_box<CoordinateType> &bbox,
                                 const bounding_box<CoordinateType> &boundaries,
                                 Func                                func,
                                 uint64_t                            node_index) const {
        assert(node_index < nodes_.size());

        const tree_node &node = nodes_[node_index];
        if (node.is_a_leaf()) {
            for (uint64_t i = 0; i < node.leaf.size; ++i) {
                if (is_inside_bounding_box(node.leaf.items[i].x, node.leaf.items[i].y, bbox)) {
                    func(iterator(this, node_index, i));
                }
            }
            return;
        }

        __unroll_4 for (uint64_t i = 0; i < 4; ++i) {
            auto new_boundaries = compute_new_boundaries(i, boundaries);
            if (bounding_boxes_overlap(bbox, new_boundaries)) {
                find_recursively(bbox, new_boundaries, func, node.branch.index_of_first_child + i);
            }
        }
    }

    static inline CoordinateType smallest_distance_from_bounding_box(
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
            for (uint64_t i = 0; i < node.leaf.size; ++i) {
                auto distance =
                    euclidean_distance_squared(x, node.leaf.items[i].x, y, node.leaf.items[i].y);
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

        std::array<bounding_box<CoordinateType>, 4> new_boundaries;
        __unroll_4 for (uint64_t i = 0; i < 4; ++i) {
            new_boundaries[i] = compute_new_boundaries(i, boundaries);
        }

        const auto selected_quad = belongs_to_quadrant(boundaries, x, y);
        nearest_recursively(node.branch.index_of_first_child + selected_quad,
                            new_boundaries[selected_quad],
                            x,
                            y,
                            nearest_distance_squared,
                            results);

        __unroll_3 for (uint64_t i = 1; i <= 3; ++i) {
            uint64_t quad = (i + selected_quad) % 4;
            if (smallest_distance_from_bounding_box(new_boundaries[quad], x, y) <=
                nearest_distance_squared) {
                nearest_recursively(node.branch.index_of_first_child + quad,
                                    new_boundaries[quad],
                                    x,
                                    y,
                                    nearest_distance_squared,
                                    results);
            }
        }
    }

    const bounding_box<CoordinateType> boundaries_;
    std::vector<tree_node>             nodes_;
    std::vector<uint64_t>              freed_nodes_;
};
}  // namespace internal

template <typename CoordinateType, typename StorageType, uint8_t MAXIMUM_NODE_SIZE = 32>
class spatial_map : public internal::spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> {
public:
    static_assert(!std::is_void_v<StorageType>, "For no storage type, use st::spatial_set");
    spatial_map() : internal::spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE>() {}
    spatial_map(const bounding_box<CoordinateType> &boundaries)
        : internal::spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE>(boundaries) {}
};

template <typename CoordinateType = double, uint8_t MAXIMUM_NODE_SIZE = 32>
class spatial_set : public internal::spatial_tree<void, CoordinateType, MAXIMUM_NODE_SIZE> {
public:
    spatial_set() : internal::spatial_tree<void, CoordinateType, MAXIMUM_NODE_SIZE>() {}
    spatial_set(const bounding_box<CoordinateType> &boundaries)
        : internal::spatial_tree<void, CoordinateType, MAXIMUM_NODE_SIZE>(boundaries) {}
};

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