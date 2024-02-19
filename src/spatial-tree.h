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

namespace internal {

template <uint64_t factor, typename Func, typename InductionVarType = uint64_t>
inline void unroll_for(InductionVarType start, InductionVarType stop, const Func &func) {
    if constexpr (factor == 2) {
        __unroll_2 for (InductionVarType i = start; i < stop; ++i) { func(i); }
    } else if (factor == 3) {
        __unroll_3 for (InductionVarType i = start; i < stop; ++i) { func(i); }
    } else if (factor == 4) {
        __unroll_4 for (InductionVarType i = start; i < stop; ++i) { func(i); }
    } else if (factor == 8) {
        __unroll_8 for (InductionVarType i = start; i < stop; ++i) { func(i); }
    } else {
        __unroll_2 for (InductionVarType i = start; i < stop; ++i) { func(i); }
    }
}

template <typename CoordinateType, uint64_t rank>
inline bool equal(std::array<CoordinateType, rank> lhs, std::array<CoordinateType, rank> rhs) {
    bool same = true;
    unroll_for<rank>(uint64_t(0), rank, [&](auto i) { same &= (lhs[i] == rhs[i]); });

    return same;
}

template <typename CoordinateType, uint64_t rank>
inline void set(std::array<CoordinateType, rank> &dst, std::array<CoordinateType, rank> src) {
    unroll_for<rank>(uint64_t(0), rank, [&](auto i) { dst[i] = src[i]; });
}

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

template <typename CoordinateType, bool strict = false>
static inline bool is_within_interval(CoordinateType x, CoordinateType beg, CoordinateType end) {
    assert(end >= beg);

    if constexpr (strict) {
        return x > beg && x < end;
    } else {
        return x >= beg && x <= end;
    }
}

template <typename CoordinateType, bool strict = false>
static inline bool intervals_overlap(CoordinateType lhs_beg,
                                     CoordinateType lhs_end,
                                     CoordinateType rhs_beg,
                                     CoordinateType rhs_end) {
    assert(lhs_beg <= lhs_end);
    assert(rhs_beg <= rhs_end);

    if constexpr (strict) {
        return !((lhs_beg >= rhs_end) || (lhs_end <= rhs_beg));
    } else {
        return !((lhs_beg > rhs_end) || (lhs_end < rhs_beg));
    }
}

}  // namespace internal

// When replacing the other bounding_box, be sure to check operand order.
template <typename CoordinateType, uint64_t rank = 2>
struct __bounding_box {
    static_assert(rank > 0 && rank <= sizeof(uint64_t) * 8);

    std::array<CoordinateType, rank> starts;
    std::array<CoordinateType, rank> stops;

    inline __bounding_box() {
        starts.fill(std::numeric_limits<CoordinateType>::lowest() / 2);
        stops.fill(std::numeric_limits<CoordinateType>::max() / 2);
    }

    inline __bounding_box(std::array<CoordinateType, rank * 2> boundary) {
        internal::unroll_for<rank>(uint64_t(0), rank, [&](uint64_t i) {
            assert(boundary[i] >= std::numeric_limits<CoordinateType>::lowest() / 2);
            assert(boundary[i + rank] <= std::numeric_limits<CoordinateType>::max() / 2);

            starts[i] = boundary[i];
            stops[i] = boundary[i + rank];
        });
    }

    inline bool operator==(const __bounding_box<CoordinateType, rank> &other) const {
        return !std::memcmp(starts.begin(), other.starts.begin(), sizeof(starts)) &&
               !std::memcmp(stops.begin(), other.stops.begin(), sizeof(stops));
    }

    inline CoordinateType area() const {
        CoordinateType area = 1;
        internal::unroll_for<rank>(
            uint64_t(0), rank, [&](auto i) { area *= (stops[i] - starts[i]); });

        return area;
    }

    template <bool strict = false>
    inline bool contains(std::array<CoordinateType, rank> point) const {
        for (uint64_t i = 0; i < rank; ++i) {
            if (!internal::is_within_interval<CoordinateType, strict>(
                    point[i], starts[i], stops[i])) {
                return false;
            }
        }

        return true;
    }

    template <bool strict = false>
    inline bool overlaps(const __bounding_box<CoordinateType, rank> &other) const {
        for (uint64_t i = 0; i < rank; ++i) {
            if (!internal::intervals_overlap<CoordinateType, strict>(
                    starts[i], stops[i], other.starts[i], other.stops[i])) {
                return false;
            }
        }

        return true;
    }

    inline std::array<CoordinateType, rank> origin() const {
        std::array<CoordinateType, rank> origins;
        internal::unroll_for<rank>(uint64_t(0), rank, [&](auto i) {
            origins[i] = (stops[i] - starts[i]) / 2 + starts[i];
        });

        return origins;
    }

    inline uint64_t quadrant(std::array<CoordinateType, rank> point) const {
        assert(contains(point));

        // Contract: lower region is ceiled, upper region is floored. (<)
        uint64_t   quadrant = 0;
        const auto origins = origin();
        internal::unroll_for<rank>(
            uint64_t(0), rank, [&](auto i) { quadrant |= (point[i] > origins[i]) << i; });

        return quadrant;
    }

    inline std::tuple<__bounding_box<CoordinateType, rank>, uint64_t> recurse(
        std::array<CoordinateType, rank> point) const {
        assert(contains(point));

        const auto                           origins = origin();
        uint64_t                             quad = 0;
        std::array<CoordinateType, rank * 2> boundary;
        internal::unroll_for<rank>(uint64_t(0), rank, [&](auto i) {
            CoordinateType span = stops[i] - starts[i];
            CoordinateType new_span = span / 2;
            CoordinateType remainder = span - 2 * new_span;

            bool gt = point[i] > origins[i];
            if (gt) {
                boundary[i] = starts[i] + new_span;
                boundary[i + rank] = stops[i];
                quad |= (1 << i);
            } else {
                boundary[i] = starts[i];
                boundary[i + rank] = stops[i] - (new_span + remainder);
            }
        });

        return {__bounding_box<CoordinateType, rank>(boundary), quad};
    }

    inline __bounding_box<CoordinateType, rank> qrecurse(uint64_t quad) const {
        assert(quad < std::pow(2, rank));
        std::array<CoordinateType, rank * 2> boundary;
        internal::unroll_for<rank>(uint64_t(0), rank, [&](auto i) {
            CoordinateType span = stops[i] - starts[i];
            CoordinateType new_span = span / 2;
            CoordinateType remainder = span - 2 * new_span;

            if (quad & (1 << i)) {
                boundary[i] = starts[i] + new_span;
                boundary[i + rank] = stops[i];
            } else {
                boundary[i] = starts[i];
                boundary[i + rank] = stops[i] - (new_span + remainder);
            }
        });

        return __bounding_box<CoordinateType, rank>(boundary);
    }

    inline CoordinateType sdistance(std::array<CoordinateType, rank> point) const {
        assert(!contains<true>(point));

        CoordinateType distance_squared = 0;
        internal::unroll_for<rank>(uint64_t(0), rank, [&](auto i) {
            CoordinateType gt = point[i] - stops[i];
            CoordinateType lt = starts[i] - point[i];
            CoordinateType d = std::max(gt, lt);
            distance_squared += d * d;
        });

        return distance_squared;
    }
};

namespace internal {

template <typename StorageType = void,
          typename CoordinateType = int,
          uint64_t rank = 2,
          uint8_t  MAXIMUM_NODE_SIZE = 32>
class __spatial_tree {
public:
    static_assert(rank > 0 && rank <= sizeof(uint64_t) * 8,
                  "Rank must be greater than 0 and less than 64");
    static_assert(MAXIMUM_NODE_SIZE > 0, "Maximum node size must be greater than 1");

    __spatial_tree() { clear(); }
    __spatial_tree(const __bounding_box<CoordinateType, rank> &boundary) : boundary_(boundary) {
        clear();
    }

    ~__spatial_tree() = default;

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
        const std::function<void(uint64_t, const __bounding_box<CoordinateType> &)>
            walk_recursively = [&](auto node_index, auto boundary) {
                assert(node_index < nodes_.size());
                const tree_node &node = nodes_[node_index];
                func(boundary, node.is_a_leaf());

                if (node.is_a_leaf()) {
                    return;
                }

                for (uint64_t child = 0; child < std::pow(2, rank); ++child) {
                    const uint64_t child_index = node.branch.index_of_first_child + child;
                    const auto     new_boundary = compute_new_boundary(child, boundary);
                    walk_recursively(child_index, new_boundary);
                }
            };

        walk_recursively(0, boundary_);
    }

    struct iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = iterator;
        using pointer = value_type *;
        using reference = value_type &;

        inline iterator(
            const __spatial_tree<StorageType, CoordinateType, rank, MAXIMUM_NODE_SIZE> *tree,
            uint64_t                                                                    node_index,
            uint8_t                                                                     item_index)
            : tree_(tree), node_index_(node_index), item_index_(item_index) {}
        inline iterator(
            const __spatial_tree<StorageType, CoordinateType, rank, MAXIMUM_NODE_SIZE> *tree)
            : tree_(tree) {
            end();
        }
        inline ~iterator() = default;

        inline auto operator*() const { return tree_->operator()(node_index_, item_index_); }
        inline auto operator*() {
            return const_cast<__spatial_tree *>(tree_)->operator()(node_index_, item_index_);
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

        const __spatial_tree<StorageType, CoordinateType, rank, MAXIMUM_NODE_SIZE> *tree_;

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
    inline std::pair<iterator, bool> emplace(std::array<CoordinateType, rank> point,
                                             Args &&...args) {
        assert(boundary_.contains(point));

        return emplace_recursively(0, boundary_, point, std::forward<Args>(args)...);
    }

private:
    using cartesian_quadrant = int;

    struct tree_node_with_storage {
        std::array<CoordinateType, rank> coordinates;
        StorageType                      storage;

        template <typename... Args>
        inline tree_node_with_storage(std::array<CoordinateType, rank> coords, Args &&...args)
            : coordinates(coords), storage(std::forward<Args>(args)...) {}
    };

    struct tree_node_no_storage {
        std::array<CoordinateType, rank> coordinates;
        inline tree_node_no_storage(std::array<CoordinateType, rank> coords)
            : coordinates(coords) {}
    };

    struct tree_node_storage : std::conditional<std::is_void_v<StorageType>,
                                                tree_node_no_storage,
                                                tree_node_with_storage>::type {};

    struct tree_node_leaf {
        uint8_t                                          size;
        std::array<tree_node_storage, MAXIMUM_NODE_SIZE> items;

        tree_node_leaf(const tree_node_leaf &other) {
            for (uint64_t i = 0; i < size; ++i) {
                items[i] = other.items[i];
            }
            size = other.size;
        }

        tree_node_leaf(tree_node_leaf &&other) {
            for (uint64_t i = 0; i < size; ++i) {
                items[i] = std::move(other.items[i]);
            }
            size = other.size;
        }

        tree_node_leaf() = default;
        ~tree_node_leaf() {
            for (uint64_t i = 0; i < size; ++i) {
                items[i].~tree_node_storage();
            }
        }

        tree_node_leaf &operator=(tree_node_leaf &&) = default;
    };

    struct tree_node_branch {
        // Children are adjacent in memory.
        uint64_t index_of_first_child;
        uint64_t size;
    };

    struct tree_node {
        static constexpr uint64_t NO_INDEX = std::numeric_limits<uint64_t>::max();

        bool is_branch;
        union {
            tree_node_leaf   leaf;
            tree_node_branch branch;
        };

        inline tree_node() { reset(); }
        inline tree_node(tree_node &&other) : is_branch(other.is_branch) {
            if (is_a_branch()) {
                branch = std::move(other.branch);
            } else {
                leaf = std::move(other.leaf);
            }
        }
        inline ~tree_node() {
            if (is_a_leaf()) {
                leaf.~tree_node_leaf();
            }
        }
        inline bool is_a_branch() const { return is_branch; }
        inline bool is_a_leaf() const { return !is_a_branch(); }
        inline void reset() {
            leaf.size = 0;
            is_branch = false;
        }
    };

    template <typename... Args>
    inline std::pair<iterator, bool> emplace_recursively_helper(
        uint64_t                              index_of_first_child,
        const __bounding_box<CoordinateType> &boundary,
        std::array<CoordinateType, rank>      point,
        Args &&...args) {
        const auto [new_boundary, selected_quad] = boundary.recurse(point);

        return emplace_recursively(
            index_of_first_child + selected_quad, new_boundary, point, std::forward<Args>(args)...);
    }

    template <typename... Args>
    std::pair<iterator, bool> emplace_recursively(uint64_t                              node_index,
                                                  const __bounding_box<CoordinateType> &boundary,
                                                  std::array<CoordinateType, rank>      point,
                                                  Args &&...args) {
        assert(node_index < nodes_.size());
        assert(boundary.contains(point));

        tree_node &node = nodes_[node_index];
        if (node.is_a_branch()) {
            const auto ret = emplace_recursively_helper(
                node.branch.index_of_first_child, boundary, point, std::forward<Args>(args)...);
            nodes_[node_index].branch.size += ret.second;

            return ret;
        }

        // SIMD this shit.
        uint64_t item_index = 0;
        __unroll_4 for (uint64_t i = 0; i < node.leaf.size; ++i) {
            const bool same =
                internal::equal<CoordinateType, rank>(point, node.leaf.items[i].coordinates);
            item_index += (i + 1) * same;
        }

        if (item_index) {
            return {iterator(this, node_index, item_index - 1), false};
        }

        const bool node_is_full = node.leaf.size == MAXIMUM_NODE_SIZE;
        if (!node_is_full) {
            internal::set<CoordinateType, rank>(node.leaf.items[node.leaf.size].coordinates, point);
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
            nodes_.resize(nodes_.size() + std::pow(2, rank));
        }

        // nodes_.resize may have reallocated.
        tree_node &node_as_branch = nodes_[node_index];
        internal::unroll_for<MAXIMUM_NODE_SIZE>(uint8_t(0), MAXIMUM_NODE_SIZE, [&](auto i) {
            if constexpr (std::is_void_v<StorageType>) {
                emplace_recursively_helper(
                    new_index_of_first_child, boundary, node_as_branch.leaf.items[i].coordinates);
            } else {
                emplace_recursively_helper(new_index_of_first_child,
                                           boundary,
                                           node_as_branch.leaf.items[i].coordinates,
                                           node_as_branch.leaf.items[i].storage);
            }
        });

        node_as_branch.is_branch = true;
        node_as_branch.branch.index_of_first_child = new_index_of_first_child;
        node_as_branch.branch.size = MAXIMUM_NODE_SIZE + 1;

        const auto ret = emplace_recursively_helper(
            new_index_of_first_child, boundary, point, std::forward<Args>(args)...);
        assert(ret.second);

        return ret;
    }

    const __bounding_box<CoordinateType, rank> boundary_;
    std::vector<tree_node>                     nodes_;
    std::vector<uint64_t>                      freed_nodes_;
};

}  // namespace internal

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

template <typename CoordinateType>
static inline bool bounding_boxes_overlap(const bounding_box<CoordinateType> &lhs,
                                          const bounding_box<CoordinateType> &rhs) {
    return intervals_overlap(lhs.top_x, lhs.bottom_x, rhs.top_x, rhs.bottom_x) &&
           intervals_overlap(lhs.bottom_y, lhs.top_y, rhs.bottom_y, rhs.top_y);
}

template <typename CoordinateType>
static inline bool is_inside_bounding_box(CoordinateType                      x,
                                          CoordinateType                      y,
                                          const bounding_box<CoordinateType> &bbox) {
    return is_within_interval(x, bbox.top_x, bbox.bottom_x) &&
           is_within_interval(y, bbox.bottom_y, bbox.top_y);
}

template <typename StorageType = void,
          typename CoordinateType = int,
          uint8_t MAXIMUM_NODE_SIZE = 32>
class spatial_tree {
public:
    static_assert(MAXIMUM_NODE_SIZE > 0, "Maximum node size must be greater than 0");

    spatial_tree()
        // Dividing by 2 to be able to do (upper bound +/- lower bound) without overflowing.
        : boundary_(std::numeric_limits<CoordinateType>::lowest() / 2,
                    std::numeric_limits<CoordinateType>::max() / 2,
                    std::numeric_limits<CoordinateType>::max() / 2,
                    std::numeric_limits<CoordinateType>::lowest() / 2) {
        clear();
    }
    spatial_tree(const bounding_box<CoordinateType> &boundary) : boundary_(boundary) {
        assert(boundary_.top_x >= std::numeric_limits<CoordinateType>::lowest() / 2);
        assert(boundary_.top_y <= std::numeric_limits<CoordinateType>::max() / 2);
        assert(boundary_.bottom_x <= std::numeric_limits<CoordinateType>::max() / 2);
        assert(boundary_.bottom_y >= std::numeric_limits<CoordinateType>::lowest() / 2);
        clear();
    }

    ~spatial_tree() = default;

    inline void     reserve(uint64_t capacity) { nodes_.reserve(4 * capacity / MAXIMUM_NODE_SIZE); }
    inline uint64_t capacity() const { return nodes_.capacity() / 4 * MAXIMUM_NODE_SIZE; }
    inline uint64_t size() const {
        assert(nodes_.size());
        const auto &root = nodes_.front();
        if (root.is_a_leaf()) {
            return root.leaf.size;
        }
        return root.branch.size;
    }
    inline uint64_t bsize() const {
        return sizeof(boundary_) + sizeof(tree_node) * (nodes_.size() + freed_nodes_.size());
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
            [&](auto node_index, auto boundary) {
                assert(node_index < nodes_.size());
                const tree_node &node = nodes_[node_index];
                func(boundary, node.is_a_leaf());

                if (node.is_a_leaf()) {
                    return;
                }

                for (uint64_t quad = 0; quad < 4; ++quad) {
                    const uint64_t child_index = node.branch.index_of_first_child + quad;
                    const auto     new_boundary = compute_new_boundary(quad, boundary);
                    walk_recursively(child_index, new_boundary);
                }
            };

        walk_recursively(0, boundary_);
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
        assert(is_inside_bounding_box(x, y, boundary_));

        return emplace_recursively(0, boundary_, x, y, std::forward<Args>(args)...);
    }

    inline bool erase(CoordinateType x, CoordinateType y) {
        return erase_recursively(0, boundary_, size(), x, y);
    }

    std::vector<iterator> nearest(CoordinateType x, CoordinateType y) const {
        CoordinateType        nearest_distance_squared = std::numeric_limits<CoordinateType>::max();
        std::vector<iterator> results;
        nearest_recursively(0, boundary_, x, y, nearest_distance_squared, results);

        return results;
    }

    template <typename Func>
    inline void find(const bounding_box<CoordinateType> &bbox, Func func) const {
        find_recursively(bbox, boundary_, func, 0);
    }

    inline iterator find(CoordinateType x, CoordinateType y) const {
        return find_recursively(x, y, boundary_, 0);
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
        uint8_t                                          size;
        std::array<tree_node_storage, MAXIMUM_NODE_SIZE> items;

        tree_node_leaf(const tree_node_leaf &other) {
            for (uint64_t i = 0; i < size; ++i) {
                items[i] = other.items[i];
            }
            size = other.size;
        }

        tree_node_leaf(tree_node_leaf &&other) {
            for (uint64_t i = 0; i < size; ++i) {
                items[i] = std::move(other.items[i]);
            }
            size = other.size;
        }

        tree_node_leaf() = default;
        ~tree_node_leaf() {
            for (uint64_t i = 0; i < size; ++i) {
                items[i].~tree_node_storage();
            }
        }

        tree_node_leaf &operator=(tree_node_leaf &&) = default;
    };

    struct tree_node_branch {
        // Children are adjacent in memory.
        uint64_t index_of_first_child;
        uint64_t size;
    };

    struct tree_node {
        static constexpr uint64_t NO_INDEX = std::numeric_limits<uint64_t>::max();

        bool is_branch;
        union {
            tree_node_leaf   leaf;
            tree_node_branch branch;
        };

        inline tree_node() { reset(); }
        inline tree_node(tree_node &&other) : is_branch(other.is_branch) {
            if (is_a_branch()) {
                branch = std::move(other.branch);
            } else {
                leaf = std::move(other.leaf);
            }
        }
        inline ~tree_node() {
            if (is_a_leaf()) {
                leaf.~tree_node_leaf();
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

    static inline auto belongs_to_quadrant(const bounding_box<CoordinateType> &boundary,
                                           CoordinateType                      x,
                                           CoordinateType                      y) {
        CoordinateType origin_x = (boundary.bottom_x - boundary.top_x) / 2 + boundary.top_x;
        CoordinateType origin_y = (boundary.top_y - boundary.bottom_y) / 2 + boundary.bottom_y;

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

    static inline bounding_box<CoordinateType> compute_new_boundary(
        uint64_t quad, const bounding_box<CoordinateType> &boundary) {
        static std::array<std::tuple<bool, bool>, 4> factors;
        factors[cartesian_quadrant::NE] = std::tuple<bool, bool>{1, 1};
        factors[cartesian_quadrant::NW] = std::tuple<bool, bool>{0, 1};
        factors[cartesian_quadrant::SW] = std::tuple<bool, bool>{0, 0};
        factors[cartesian_quadrant::SE] = std::tuple<bool, bool>{1, 0};
        const auto [fx, fy] = factors[quad];

        CoordinateType width = (boundary.bottom_x - boundary.top_x);
        CoordinateType height = (boundary.top_y - boundary.bottom_y);
        CoordinateType new_width = width / 2;
        CoordinateType new_height = height / 2;
        CoordinateType width_remainder = width - 2 * new_width;
        CoordinateType height_remainder = height - 2 * new_height;

        CoordinateType top_x = boundary.top_x;
        CoordinateType bottom_x = boundary.bottom_x;
        CoordinateType top_y = boundary.top_y;
        CoordinateType bottom_y = boundary.bottom_y;

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
        const bounding_box<CoordinateType> &boundary,
        CoordinateType                      x,
        CoordinateType                      y,
        Args &&...args) {
        const auto selected_quad = belongs_to_quadrant(boundary, x, y);
        const auto new_boundary = compute_new_boundary(selected_quad, boundary);

        return emplace_recursively(
            index_of_first_child + selected_quad, new_boundary, x, y, std::forward<Args>(args)...);
    }

    template <typename... Args>
    std::pair<iterator, bool> emplace_recursively(uint64_t                            node_index,
                                                  const bounding_box<CoordinateType> &boundary,
                                                  CoordinateType                      x,
                                                  CoordinateType                      y,
                                                  Args &&...args) {
        assert(node_index < nodes_.size());
        assert(is_inside_bounding_box(x, y, boundary));

        tree_node &node = nodes_[node_index];

        if (node.is_a_branch()) {
            const auto ret = emplace_recursively_helper(
                node.branch.index_of_first_child, boundary, x, y, std::forward<Args>(args)...);
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
                emplace_recursively_helper(new_index_of_first_child, boundary, x_, y_);
            } else {
                emplace_recursively_helper(new_index_of_first_child,
                                           boundary,
                                           x_,
                                           y_,
                                           node_as_branch.leaf.items[i].storage);
            }
        }

        node_as_branch.is_branch = true;
        node_as_branch.branch.index_of_first_child = new_index_of_first_child;
        node_as_branch.branch.size = MAXIMUM_NODE_SIZE + 1;

        const auto ret = emplace_recursively_helper(
            new_index_of_first_child, boundary, x, y, std::forward<Args>(args)...);
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
                                  const bounding_box<CoordinateType> &boundary,
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
        const auto selected_quad = belongs_to_quadrant(boundary, x, y);
        const bool erased = erase_recursively(node.branch.index_of_first_child + selected_quad,
                                              compute_new_boundary(selected_quad, boundary),
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
                                     const bounding_box<CoordinateType> &boundary,
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

        const auto selected_quad = belongs_to_quadrant(boundary, x, y);
        return find_recursively(x,
                                y,
                                compute_new_boundary(selected_quad, boundary),
                                node.branch.index_of_first_child + selected_quad);
    }

    template <typename Func>
    inline void find_recursively(const bounding_box<CoordinateType> &bbox,
                                 const bounding_box<CoordinateType> &boundary,
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
            auto new_boundary = compute_new_boundary(i, boundary);
            if (bounding_boxes_overlap(bbox, new_boundary)) {
                find_recursively(bbox, new_boundary, func, node.branch.index_of_first_child + i);
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
                             const bounding_box<CoordinateType> &boundary,
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

        std::array<bounding_box<CoordinateType>, 4> new_boundary;
        __unroll_4 for (uint64_t i = 0; i < 4; ++i) {
            new_boundary[i] = compute_new_boundary(i, boundary);
        }

        const auto selected_quad = belongs_to_quadrant(boundary, x, y);
        nearest_recursively(node.branch.index_of_first_child + selected_quad,
                            new_boundary[selected_quad],
                            x,
                            y,
                            nearest_distance_squared,
                            results);

        __unroll_3 for (uint64_t i = 1; i <= 3; ++i) {
            uint64_t quad = (i + selected_quad) % 4;
            if (smallest_distance_from_bounding_box(new_boundary[quad], x, y) <=
                nearest_distance_squared) {
                nearest_recursively(node.branch.index_of_first_child + quad,
                                    new_boundary[quad],
                                    x,
                                    y,
                                    nearest_distance_squared,
                                    results);
            }
        }
    }

    const bounding_box<CoordinateType> boundary_;
    std::vector<tree_node>             nodes_;
    std::vector<uint64_t>              freed_nodes_;
};
}  // namespace internal

template <typename CoordinateType, typename StorageType, uint8_t MAXIMUM_NODE_SIZE = 32>
class spatial_map : public internal::spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> {
public:
    static_assert(!std::is_void_v<StorageType>, "For no storage type, use st::spatial_set");
    spatial_map() : internal::spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE>() {}
    spatial_map(const bounding_box<CoordinateType> &boundary)
        : internal::spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE>(boundary) {}
};

template <typename CoordinateType = double, uint8_t MAXIMUM_NODE_SIZE = 32>
class spatial_set : public internal::spatial_tree<void, CoordinateType, MAXIMUM_NODE_SIZE> {
public:
    spatial_set() : internal::spatial_tree<void, CoordinateType, MAXIMUM_NODE_SIZE>() {}
    spatial_set(const bounding_box<CoordinateType> &boundary)
        : internal::spatial_tree<void, CoordinateType, MAXIMUM_NODE_SIZE>(boundary) {}
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