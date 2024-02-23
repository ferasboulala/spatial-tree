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

#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <initializer_list>
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

template <typename T>
inline constexpr T pow(T base, uint64_t exponent) {
    for (uint64_t i = 1; i < exponent; ++i) {
        base *= exponent;
    }

    return base;
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
inline T euclidean_distance_squared_impl(AbsDiff) {
    return T(0);
}

template <typename AbsDiff, typename T, typename... Args>
inline T euclidean_distance_squared_impl(AbsDiff absdiff, T x1, T x2, Args... xs) {
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
inline T absdiff(T x, T y) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::fabs(UnsafeAbsDiff<T>()(x, y));
    } else {
        return SafeAbsDiff<T>()(x, y);
    }
}

template <typename T, typename... Args>
inline T euclidean_distance_squared(T x1, T x2, Args... xs) {
    if constexpr (std::is_floating_point_v<T> || std::is_signed_v<T>) {
        return euclidean_distance_squared_impl(UnsafeAbsDiff<T>(), x1, x2, xs...);
    } else {
        return euclidean_distance_squared_impl(SafeAbsDiff<T>(), x1, x2, xs...);
    }
}

template <typename AbsDiff, typename T, uint64_t rank>
inline T euclidean_distance_squared_arr_impl(AbsDiff             absdiff,
                                             std::array<T, rank> lhs,
                                             std::array<T, rank> rhs) {
    T sum = 0;
    unroll_for<rank>(uint64_t(0), rank, [&](auto i) {
        const T dx = absdiff(lhs[i], rhs[i]);
        sum += dx * dx;
    });

    return sum;
}

template <typename T, uint64_t rank>
inline T euclidean_distance_squared_arr(std::array<T, rank> lhs, std::array<T, rank> rhs) {
    if constexpr (std::is_floating_point_v<T> || std::is_signed_v<T>) {
        return euclidean_distance_squared_arr_impl<UnsafeAbsDiff<T>, T, rank>(
            UnsafeAbsDiff<T>(), lhs, rhs);
    } else {
        return euclidean_distance_squared_arr_impl<SafeAbsDiff<T>, T, rank>(
            SafeAbsDiff<T>(), lhs, rhs);
    }
}

template <typename CoordinateType, bool strict = false>
inline bool is_within_interval(CoordinateType x, CoordinateType beg, CoordinateType end) {
    assert(end >= beg);

    if constexpr (strict) {
        return x > beg && x < end;
    } else {
        return x >= beg && x <= end;
    }
}

template <typename CoordinateType, bool strict = false>
inline bool intervals_overlap(CoordinateType lhs_beg,
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
struct bounding_box {
    static_assert(rank > 0 && rank <= sizeof(uint64_t) * 8);

    static constexpr uint64_t BRANCHING_FACTOR = internal::pow(2, rank);

    std::array<CoordinateType, rank> starts;
    std::array<CoordinateType, rank> stops;

    inline bounding_box() {
        starts.fill(std::numeric_limits<CoordinateType>::lowest() / 2);
        stops.fill(std::numeric_limits<CoordinateType>::max() / 2);
    }

    inline bounding_box(std::array<CoordinateType, rank * 2> boundary) {
        internal::unroll_for<rank>(uint64_t(0), rank, [&](uint64_t i) {
            assert(boundary[i] >= std::numeric_limits<CoordinateType>::lowest() / 2);
            assert(boundary[i + rank] <= std::numeric_limits<CoordinateType>::max() / 2);

            starts[i] = boundary[i];
            stops[i] = boundary[i + rank];
        });
    }

    inline bounding_box(std::initializer_list<CoordinateType> boundary) {
        assert(boundary.size() == 2 * rank);
        internal::unroll_for<rank>(uint64_t(0), rank, [&](uint64_t i) {
            assert(*(boundary.begin() + i) >= std::numeric_limits<CoordinateType>::lowest() / 2);
            assert(*(boundary.begin() + i + rank) <=
                   std::numeric_limits<CoordinateType>::max() / 2);

            starts[i] = *(boundary.begin() + i);
            stops[i] = *(boundary.begin() + i + rank);
        });
    }

    inline bool operator==(const bounding_box<CoordinateType, rank> &other) const {
        return !std::memcmp(this, &other, sizeof(other));
    }

    inline CoordinateType area() const {
        CoordinateType area = 1;
        internal::unroll_for<rank>(
            uint64_t(0), rank, [&](auto i) { area *= stops[i] - starts[i]; });

        return area;
    }

    template <bool strict = false>
    inline bool contains(std::array<CoordinateType, rank> point) const {
        bool is_contained = true;
        internal::unroll_for<rank>(uint64_t(0), rank, [&](auto i) {
            is_contained &=
                internal::is_within_interval<CoordinateType, strict>(point[i], starts[i], stops[i]);
        });

        return is_contained;
    }

    template <bool strict = false>
    inline bool overlaps(const bounding_box<CoordinateType, rank> &other) const {
        bool is_overlapping = true;
        internal::unroll_for<rank>(uint64_t(0), rank, [&](auto i) {
            is_overlapping &= internal::intervals_overlap<CoordinateType, strict>(
                starts[i], stops[i], other.starts[i], other.stops[i]);
        });

        return is_overlapping;
    }

    inline std::array<CoordinateType, rank> origin() const {
        std::array<CoordinateType, rank> origins;
        internal::unroll_for<rank>(uint64_t(0), rank, [&](auto i) {
            origins[i] = (stops[i] - starts[i]) / 2 + starts[i];
        });

        return origins;
    }

    inline uint64_t quadrant(std::array<CoordinateType, rank> point) const {
        // Contract: lower region is ceiled, upper region is floored. (>)
        // Assumes the bbox is infinite.
        uint64_t   quadrant = 0;
        const auto origins = origin();
        internal::unroll_for<rank>(
            uint64_t(0), rank, [&](auto i) { quadrant |= (point[i] > origins[i]) << i; });

        return quadrant;
    }

    inline std::tuple<bounding_box<CoordinateType, rank>, uint64_t> recurse(
        std::array<CoordinateType, rank> point) const {
        assert(contains(point));

        const auto                           origins = origin();
        uint64_t                             quad = 0;
        std::array<CoordinateType, rank * 2> boundary;
        internal::unroll_for<rank>(uint64_t(0), rank, [&](auto i) {
            CoordinateType span = stops[i] - starts[i];
            CoordinateType new_span = span / 2;
            // Contract: lower region is ceiled, upper region is floored. (>)
            bool gt = point[i] > origins[i];

            // Faster than branchless.
            if (gt) {
                boundary[i] = starts[i] + new_span;
                boundary[i + rank] = stops[i];
                quad |= (1 << i);
            } else {
                CoordinateType remainder = span - 2 * new_span;
                boundary[i] = starts[i];
                boundary[i + rank] = stops[i] - (new_span + remainder);
            }
        });

        return {bounding_box<CoordinateType, rank>(boundary), quad};
    }

    inline bounding_box<CoordinateType, rank> qrecurse(uint64_t quad) const {
        assert(quad < BRANCHING_FACTOR);

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

        return bounding_box<CoordinateType, rank>(boundary);
    }

    inline CoordinateType sdistance(std::array<CoordinateType, rank> point) const {
        CoordinateType distance_squared = 0;
        internal::unroll_for<rank>(uint64_t(0), rank, [&](auto i) {
            CoordinateType gt = (point[i] > stops[i]) * (point[i] - stops[i]);
            CoordinateType lt = (point[i] < starts[i]) * (starts[i] - point[i]);
            distance_squared += gt * gt + lt * lt;
        });

        return distance_squared;
    }
};

namespace internal {

template <typename StorageType = void,
          typename CoordinateType = int,
          uint64_t rank = 2,
          uint8_t  MAXIMUM_NODE_SIZE = 32>
class spatial_tree {
public:
    static_assert(rank > 0 && rank <= sizeof(uint64_t) * 8,
                  "Rank must be greater than 0 and less than 64");
    static_assert(MAXIMUM_NODE_SIZE > 0, "Maximum node size must be greater than 1");

    static constexpr uint64_t BRANCHING_FACTOR = internal::pow(2, rank);

    spatial_tree() { clear(); }
    spatial_tree(const bounding_box<CoordinateType, rank> &boundary) : boundary_(boundary) {
        clear();
    }
    spatial_tree(std::initializer_list<CoordinateType> boundary) : boundary_(boundary) { clear(); }

    ~spatial_tree() = default;

    inline void reserve(uint64_t capacity) {
        nodes_.reserve(BRANCHING_FACTOR * capacity / MAXIMUM_NODE_SIZE);
    }
    inline uint64_t capacity() const {
        return nodes_.capacity() / BRANCHING_FACTOR * MAXIMUM_NODE_SIZE;
    }
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
            [&](auto node_index, auto boundary) {
                assert(node_index < nodes_.size());
                const tree_node &node = nodes_[node_index];
                func(boundary, node.is_a_leaf());

                if (node.is_a_leaf()) {
                    return;
                }

                for (uint64_t child = 0; child < BRANCHING_FACTOR; ++child) {
                    const uint64_t child_index = node.branch.index_of_first_child + child;
                    const auto     new_boundary = boundary.qrecurse(child);
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
            const spatial_tree<StorageType, CoordinateType, rank, MAXIMUM_NODE_SIZE> *tree,
            uint64_t                                                                  node_index,
            uint8_t                                                                   item_index)
            : tree_(tree), node_index_(node_index), item_index_(item_index) {}
        inline iterator(
            const spatial_tree<StorageType, CoordinateType, rank, MAXIMUM_NODE_SIZE> *tree)
            : tree_(tree) {
            end();
        }
        inline ~iterator() = default;

        inline auto operator*() const { return tree_->operator()(node_index_, item_index_); }
        inline auto operator*() {
            return const_cast<spatial_tree<StorageType, CoordinateType, rank, MAXIMUM_NODE_SIZE> *>(
                       tree_)
                ->
                operator()(node_index_, item_index_);
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

        const spatial_tree<StorageType, CoordinateType, rank, MAXIMUM_NODE_SIZE> *tree_;

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

    inline bool erase(std::array<CoordinateType, rank> point) {
        return erase_recursively(0, boundary_, size(), point);
    }

    template <typename Func>
    inline void find(const bounding_box<CoordinateType> &bbox, Func func) const {
        find_recursively(bbox, boundary_, func, 0);
    }

    inline iterator find(std::array<CoordinateType, rank> point) const {
        return find_recursively(point, boundary_, 0);
    }

    std::vector<iterator> nearest(std::array<CoordinateType, rank> point) const {
        CoordinateType        nearest_distance_squared = std::numeric_limits<CoordinateType>::max();
        std::vector<iterator> results;
        nearest_recursively(0, boundary_, point, nearest_distance_squared, results);

        return results;
    }

private:
    struct tree_node_with_storage {
        std::array<CoordinateType, rank> coordinates;
        StorageType                      storage;

        inline tree_node_with_storage() = default;

        template <typename... Args>
        inline tree_node_with_storage(std::array<CoordinateType, rank> coords, Args &&...args)
            : coordinates(coords), storage(std::forward<Args>(args)...) {}
    };

    struct tree_node_no_storage {
        std::array<CoordinateType, rank> coordinates;

        inline tree_node_no_storage() = default;

        inline tree_node_no_storage(std::array<CoordinateType, rank> coords)
            : coordinates(coords) {}
    };

    using tree_node_storage_impl = std::conditional<std::is_void_v<StorageType>,
                                                tree_node_no_storage,
                                                tree_node_with_storage>;

    struct tree_node_storage : tree_node_storage_impl::type {
        inline tree_node_storage() {
            if constexpr (!std::is_void_v<StorageType>) {
                this->coordinates = {};
            }
        }
    };

    struct tree_node_leaf {
        std::array<tree_node_storage, MAXIMUM_NODE_SIZE> items;
        uint8_t                                          size;

        inline tree_node_leaf(const tree_node_leaf &other) {
            for (uint64_t i = 0; i < size; ++i) {
                items[i] = other.items[i];
            }
            size = other.size;
        }

        inline tree_node_leaf(tree_node_leaf &&other) {
            for (uint64_t i = 0; i < size; ++i) {
                items[i] = std::move(other.items[i]);
            }
            size = other.size;
        }

        inline tree_node_leaf() : items{} { size = 0; }

        inline ~tree_node_leaf() {
            for (uint64_t i = 0; i < size; ++i) {
                items[i].~tree_node_storage();
            }
        }

        inline tree_node_leaf &operator=(tree_node_leaf &&) = default;
    };

    struct tree_node_branch {
        // Children are adjacent in memory.
        uint64_t size;
        uint64_t index_of_first_child;
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
            if (is_a_leaf()) {
                leaf.~tree_node_leaf();
            }
        }
        inline bool is_a_branch() const { return is_branch; }
        inline bool is_a_leaf() const { return !is_a_branch(); }
        inline void reset() {
            if (is_a_leaf()) {
                leaf.~tree_node_leaf();
            }
            leaf = tree_node_leaf();
            is_branch = false;
        }
    };

    inline auto operator()(uint64_t node_index, uint64_t item_index) const {
        assert(node_index < nodes_.size());
        const tree_node &node = nodes_[node_index];

        assert(node.is_a_leaf());
        assert(item_index < node.leaf.size);

        if constexpr (std::is_void_v<StorageType>) {
            return node.leaf.items[item_index].coordinates;
        } else {
            auto               coordinates = node.leaf.items[item_index].coordinates;
            const StorageType &storage = node.leaf.items[item_index].storage;
            return std::tuple<decltype(coordinates), const StorageType &>(coordinates, storage);
        }
    }

    inline auto operator()(uint64_t node_index, uint64_t item_index) {
        auto this_as_const = const_cast<const spatial_tree *>(this);
        if constexpr (std::is_void_v<StorageType>) {
            return this_as_const->operator()(node_index, item_index);
        } else {
            const auto [coordinates, storage] = this_as_const->operator()(node_index, item_index);
            return std::tuple<decltype(coordinates), StorageType &>(
                coordinates, const_cast<StorageType &>(storage));
        }
    }

    template <typename... Args>
    inline std::pair<iterator, bool> emplace_recursively_helper(
        uint64_t                            index_of_first_child,
        const bounding_box<CoordinateType> &boundary,
        std::array<CoordinateType, rank>    point,
        Args &&...args) {
        const auto [new_boundary, selected_quad] = boundary.recurse(point);

        return emplace_recursively(
            index_of_first_child + selected_quad, new_boundary, point, std::forward<Args>(args)...);
    }

    template <typename... Args>
    std::pair<iterator, bool> emplace_recursively(uint64_t                            node_index,
                                                  const bounding_box<CoordinateType> &boundary,
                                                  std::array<CoordinateType, rank>    point,
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

        uint64_t item_index = 0;
        for (uint64_t i = 0; i < node.leaf.size; ++i) {
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
            nodes_.resize(nodes_.size() + BRANCHING_FACTOR);
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
                                  std::array<CoordinateType, rank>    point) {
        assert(node_index < nodes_.size());

        tree_node &node = nodes_[node_index];
        if (node.is_a_leaf()) {
            for (uint8_t i = 0; i < node.leaf.size; ++i) {
                if (internal::equal<CoordinateType, rank>(point, node.leaf.items[i].coordinates)) {
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

        uint64_t size_after_removal = node.branch.size - 1;
        const auto [new_boundary, selected_quad] = boundary.recurse(point);
        const bool erased = erase_recursively(node.branch.index_of_first_child + selected_quad,
                                              new_boundary,
                                              size_after_removal,
                                              point);

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
        static constexpr uint64_t N_CHILDREN = BRANCHING_FACTOR;
        internal::unroll_for<N_CHILDREN>(uint64_t(0), N_CHILDREN, [&](auto quad) {
            recursively_gather_points(index_of_first_child + quad, node.leaf);
        });
        freed_nodes_.push_back(index_of_first_child);
        assert(node.leaf.size <= MAXIMUM_NODE_SIZE);

        return true;
    }

    inline iterator find_recursively(std::array<CoordinateType, rank>    point,
                                     const bounding_box<CoordinateType> &boundary,
                                     uint64_t                            node_index) const {
        assert(node_index < nodes_.size());

        const tree_node &node = nodes_[node_index];
        if (node.is_a_leaf()) {
            for (uint64_t i = 0; i < node.leaf.size; ++i) {
                if (internal::equal<CoordinateType, rank>(point, node.leaf.items[i].coordinates)) {
                    return iterator(this, node_index, i);
                }
            }
            return end();
        }

        const auto [new_boundary, selected_quad] = boundary.recurse(point);
        return find_recursively(
            point, new_boundary, node.branch.index_of_first_child + selected_quad);
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
                if (bbox.contains(node.leaf.items[i].coordinates)) {
                    func(iterator(this, node_index, i));
                }
            }
            return;
        }

        static constexpr uint64_t N_CHILDREN = BRANCHING_FACTOR;
        internal::unroll_for<N_CHILDREN>(uint64_t(0), N_CHILDREN, [&](auto quad) {
            auto new_boundary = boundary.qrecurse(quad);
            if (bbox.overlaps(new_boundary)) {
                find_recursively(bbox, new_boundary, func, node.branch.index_of_first_child + quad);
            }
        });
    }

    void nearest_recursively(uint64_t                            node_index,
                             const bounding_box<CoordinateType> &boundary,
                             std::array<CoordinateType, rank>    point,
                             CoordinateType                     &nearest_distance_squared,
                             std::vector<iterator>              &results) const {
        assert(node_index < nodes_.size());

        const tree_node &node = nodes_[node_index];

        if (node.is_a_leaf()) {
            for (uint64_t i = 0; i < node.leaf.size; ++i) {
                auto distance = euclidean_distance_squared_arr<CoordinateType, rank>(
                    point, node.leaf.items[i].coordinates);
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

        std::array<bounding_box<CoordinateType>, BRANCHING_FACTOR> new_boundaries;
        internal::unroll_for<BRANCHING_FACTOR>(uint64_t(0), BRANCHING_FACTOR, [&](auto quad) {
            new_boundaries[quad] = boundary.qrecurse(quad);
        });

        const auto selected_quad = boundary.quadrant(point);
        nearest_recursively(node.branch.index_of_first_child + selected_quad,
                            new_boundaries[selected_quad],
                            point,
                            nearest_distance_squared,
                            results);

        internal::unroll_for<BRANCHING_FACTOR - 1>(uint64_t(1), BRANCHING_FACTOR, [&](auto i) {
            uint64_t quad = (i + selected_quad) % BRANCHING_FACTOR;
            if (new_boundaries[quad].sdistance(point) <= nearest_distance_squared) {
                nearest_recursively(node.branch.index_of_first_child + quad,
                                    new_boundaries[quad],
                                    point,
                                    nearest_distance_squared,
                                    results);
            }
        });
    }

    const bounding_box<CoordinateType, rank> boundary_;
    std::vector<tree_node>                   nodes_;
    std::vector<uint64_t>                    freed_nodes_;
};

}  // namespace internal

template <typename CoordinateType,
          typename StorageType,
          uint64_t rank = 2,
          uint8_t  MAXIMUM_NODE_SIZE = 32>
class spatial_map : public internal::spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> {
public:
    static_assert(!std::is_void_v<StorageType>, "For no storage type, use st::spatial_set");
    spatial_map()
        : internal::spatial_tree<StorageType, CoordinateType, rank, MAXIMUM_NODE_SIZE>() {}
    spatial_map(const bounding_box<CoordinateType, rank> &boundary)
        : internal::spatial_tree<StorageType, CoordinateType, rank, MAXIMUM_NODE_SIZE>(boundary) {}
};

template <typename CoordinateType = double, uint64_t rank = 2, uint8_t MAXIMUM_NODE_SIZE = 32>
class spatial_set : public internal::spatial_tree<void, CoordinateType, rank, MAXIMUM_NODE_SIZE> {
public:
    spatial_set() : internal::spatial_tree<void, CoordinateType, rank, MAXIMUM_NODE_SIZE>() {}
    spatial_set(const bounding_box<CoordinateType, rank> &boundary)
        : internal::spatial_tree<void, CoordinateType, rank, MAXIMUM_NODE_SIZE>(boundary) {}
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