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

#define NDEBUG
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "robin_hood.h"

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
inline void hash_combine(std::size_t &seed, const T &v) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <typename CoordinateType, uint64_t rank>
struct point_hash {
    inline size_t operator()(const std::array<CoordinateType, rank> &key) const {
        size_t hash = 0x778abe;
        internal::unroll_for<rank>(uint64_t(0), rank, [&](auto i) { hash_combine(hash, key[i]); });
        return hash;
    }
};

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
struct point_equal {
    inline size_t operator()(std::array<CoordinateType, rank> lhs,
                             std::array<CoordinateType, rank> rhs) const {
        return equal<CoordinateType, rank>(lhs, rhs);
    }
};

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
          uint16_t MAXIMUM_NODE_SIZE = 64>
class spatial_tree {
private:
    static_assert(rank > 0 && rank <= sizeof(uint64_t) * 8,
                  "Rank must be greater than 0 and less than 64");
    static_assert(MAXIMUM_NODE_SIZE > 0, "Maximum node size must be greater than 1");

    static constexpr uint64_t BRANCHING_FACTOR = internal::pow(2, rank);

    struct storage_pool_empty {};
    struct storage_pool_full {
        std::vector<uint64_t> vec;
    };
    struct storage_pool : std::conditional<std::is_void_v<StorageType>,
                                           storage_pool_empty,
                                           storage_pool_full>::type {};

    struct storage_container_empty {};
    struct storage_container_full {
        std::vector<StorageType> vec;
    };
    struct storage_container : std::conditional<std::is_void_v<StorageType>,
                                                storage_container_empty,
                                                storage_container_full>::type {};
    using hash_table_type = typename std::conditional<
        std::is_void_v<StorageType>,
        robin_hood::unordered_set<std::array<CoordinateType, rank>,
                                  internal::point_hash<CoordinateType, rank>,
                                  internal::point_equal<CoordinateType, rank>>,
        std::unordered_map<std::array<CoordinateType, rank>,
                           StorageType,
                           internal::point_hash<CoordinateType, rank>,
                           internal::point_equal<CoordinateType, rank>>>::type;

    using iterator = hash_table_type::iterator;
    using const_iterator = hash_table_type::const_iterator;

public:
    spatial_tree() { clear(); }
    spatial_tree(const bounding_box<CoordinateType, rank> &boundary) : boundary_(boundary) {
        clear();
    }
    spatial_tree(std::initializer_list<CoordinateType> boundary) : boundary_(boundary) { clear(); }

    ~spatial_tree() = default;

    inline void reserve(uint64_t capacity) {
        nodes_.reserve(BRANCHING_FACTOR * capacity / MAXIMUM_NODE_SIZE);
        if constexpr (!std::is_void_v<StorageType>) {
            storage_.vec.reserve(capacity);
        }
        presence_.reserve(capacity);
    }
    inline uint64_t capacity() const {
        return nodes_.capacity() / BRANCHING_FACTOR * MAXIMUM_NODE_SIZE;
    }
    inline uint64_t size() const { return presence_.size(); }
    inline bool     empty() const { return size() == 0; }
    inline void     clear() {
        nodes_.clear();
        nodes_.resize(1);
        freed_nodes_.clear();
        if constexpr (!std::is_void_v<StorageType>) {
            storage_.vec.clear();
            pool_.vec.clear();
        }
        presence_.clear();
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

    inline iterator begin() const { return presence_.begin(); }
    inline iterator end() const { return presence_.end(); }

    template <typename... Args>
    inline std::pair<iterator, bool> emplace(std::array<CoordinateType, rank> point,
                                             Args &&...args) {
        assert(boundary_.contains(point));

        if constexpr (std::is_void_v<StorageType>) {
            auto [it, inserted] = presence_.emplace(point);
            if (!inserted) {
                return {it, inserted};
            }

            emplace_recursively(0, boundary_, point);

            return {it, inserted};
        } else {
            auto [it, inserted] = presence_.emplace(point, storage_.vec.size());
            auto [_, idx] = *it;

            if (!inserted) {
                return {it, false};
            }

            if (pool_.vec.empty()) {
                storage_.vec.emplace_back(std::forward<Args>(args)...);
            } else {
                // FIXME: This will probably call the destructor.
                storage_.vec[pool_.back()] = std::move(StorageType(std::forward<Args>(args)...));
                pool_.pop_back();
            }

            emplace_recursively(0, boundary_, point, idx);

            return {it, true};
        }
    }

    inline bool erase(std::array<CoordinateType, rank> point) {
        assert(boundary_.contains(point));

        auto it = presence_.find(point);
        if (it == presence_.end()) {
            return false;
        }

        if constexpr (std::is_void_v<StorageType>) {
            presence_.erase(it);
            return true;
        } else {
            auto [_, idx] = *it;
            presence_.erase(it);

            storage_.vec[idx].~StorageType();
            pool_.push_back(idx);

            erase_recursively(0, boundary_, size(), point);

            return true;
        }
    }

    template <typename Func>
    inline void find(const bounding_box<CoordinateType> &bbox, Func func) const {
        find_recursively(bbox, boundary_, func, 0);
    }

    inline const_iterator find(std::array<CoordinateType, rank> point) const {
        return presence_.find(point);
    }

    std::vector<const_iterator> nearest(std::array<CoordinateType, rank> point) const {
        CoordinateType nearest_distance_squared = std::numeric_limits<CoordinateType>::max();
        std::vector<const_iterator> results;
        nearest_recursively(0, boundary_, point, nearest_distance_squared, results);

        return results;
    }

private:
    struct tree_node_leaf {
        std::array<std::array<CoordinateType, rank>, MAXIMUM_NODE_SIZE> items;
        uint16_t                                                        size;
        tree_node_leaf() : size(0) {}
    };

    struct tree_node_branch {
        uint64_t size;
        // Children are adjacent in memory.
        uint64_t index_of_first_child;
        tree_node_branch() : size(0) {}
    };

    struct tree_node {
        static constexpr uint64_t NO_INDEX = std::numeric_limits<uint64_t>::max();

        union {
            tree_node_leaf   leaf;
            tree_node_branch branch;
        };
        bool is_branch;

        inline tree_node() noexcept : leaf(), is_branch(false) {}
        inline bool is_a_branch() const { return is_branch; }
        inline bool is_a_leaf() const { return !is_a_branch(); }
    };

    inline void emplace_recursively_helper(uint64_t                            index_of_first_child,
                                           const bounding_box<CoordinateType> &boundary,
                                           std::array<CoordinateType, rank>    point) {
        const auto [new_boundary, selected_quad] = boundary.recurse(point);
        emplace_recursively(index_of_first_child + selected_quad, new_boundary, point);
    }

    void emplace_recursively(uint64_t                            node_index,
                             const bounding_box<CoordinateType> &boundary,
                             std::array<CoordinateType, rank>    point) {
        assert(node_index < nodes_.size());
        assert(boundary.contains(point));

        tree_node &node = nodes_[node_index];
        if (node.is_a_branch()) {
            ++node.branch.size;
            emplace_recursively_helper(node.branch.index_of_first_child, boundary, point);
            return;
        }

        if (node.leaf.size < MAXIMUM_NODE_SIZE) {
            internal::set<CoordinateType, rank>(node.leaf.items[node.leaf.size++], point);
            return;
        }

        uint64_t new_index_of_first_child;
        if (freed_nodes_.size()) {
            new_index_of_first_child = freed_nodes_.back();
            assert(nodes_[new_index_of_first_child].is_a_leaf());
            freed_nodes_.pop_back();
        } else {
            new_index_of_first_child = nodes_.size();
            nodes_.resize(nodes_.size() + BRANCHING_FACTOR);
        }
        assert((new_index_of_first_child - 1) % BRANCHING_FACTOR == 0);

        // nodes_.resize may have reallocated.
        tree_node &node_as_branch = nodes_[node_index];
        internal::unroll_for<MAXIMUM_NODE_SIZE>(uint16_t(0), MAXIMUM_NODE_SIZE, [&](auto i) {
            emplace_recursively_helper(
                new_index_of_first_child, boundary, node_as_branch.leaf.items[i]);
        });
        node_as_branch.branch.index_of_first_child = new_index_of_first_child;
        node_as_branch.branch.size = MAXIMUM_NODE_SIZE + 1;
        node_as_branch.is_branch = true;
        emplace_recursively_helper(new_index_of_first_child, boundary, point);
    }

    inline void recursively_gather_points(uint64_t node_index, tree_node_leaf &leaf) {
        assert(node_index < nodes_.size());

        tree_node &node = nodes_[node_index];
        if (node.is_a_leaf()) {
            for (uint16_t i = 0; i < node.leaf.size; ++i) {
                leaf.items[leaf.size++] = std::move(node.leaf.items[i]);
            }
            assert(leaf.size <= MAXIMUM_NODE_SIZE);
            node.leaf.size = 0;
            return;
        }

        // TODO: Specify that it is safe to access &node because there was no realloc.
        // TODO: Do it elsewhere too.
        for (uint64_t quad = 0; quad < 4; ++quad) {
            recursively_gather_points(node.branch.index_of_first_child + quad, leaf);
        }

        freed_nodes_.push_back(node.branch.index_of_first_child);
        node = tree_node();
    }

    inline void erase_recursively(uint64_t                            node_index,
                                  const bounding_box<CoordinateType> &boundary,
                                  uint64_t                            parent_size_after_removal,
                                  std::array<CoordinateType, rank>    point) {
        assert(node_index < nodes_.size());

        tree_node &node = nodes_[node_index];
        if (node.is_a_leaf()) {
            for (uint16_t i = 0; i < node.leaf.size; ++i) {
                if (internal::equal<CoordinateType, rank>(point, node.leaf.items[i])) {
                    uint64_t last_index = node.leaf.size - 1;
                    if (last_index != i) {
                        node.leaf.items[i] = std::move(node.leaf.items[last_index]);
                    }
                    --node.leaf.size;

                    return;
                }
            }
            assert(false && "Unreachable");
        }

        uint64_t size_after_removal = node.branch.size - 1;
        const auto [new_boundary, selected_quad] = boundary.recurse(point);
        erase_recursively(node.branch.index_of_first_child + selected_quad,
                          new_boundary,
                          size_after_removal,
                          point);

        if (parent_size_after_removal <= MAXIMUM_NODE_SIZE) {
            return;
        }

        if (--node.branch.size > MAXIMUM_NODE_SIZE) {
            return;
        }

        // TODO: Specify that it is safe to access &node because there was no realloc.
        /// EXPLANATION: Not doing it at the root because of the union.
        uint64_t index_of_first_child = node.branch.index_of_first_child;
        node = tree_node();
        internal::unroll_for<BRANCHING_FACTOR>(uint64_t(0), BRANCHING_FACTOR, [&](auto quad) {
            recursively_gather_points(index_of_first_child + quad, node.leaf);
        });
        freed_nodes_.push_back(index_of_first_child);
        assert(node.leaf.size <= MAXIMUM_NODE_SIZE);
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
                if (bbox.contains(node.leaf.items[i])) {
                    func(node.leaf.items[i]);
                }
            }
            return;
        }

        internal::unroll_for<BRANCHING_FACTOR>(uint64_t(0), BRANCHING_FACTOR, [&](auto quad) {
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
                             std::vector<const_iterator>        &results) const {
        assert(node_index < nodes_.size());

        const tree_node &node = nodes_[node_index];

        if (node.is_a_leaf()) {
            for (uint64_t i = 0; i < node.leaf.size; ++i) {
                auto distance =
                    euclidean_distance_squared_arr<CoordinateType, rank>(point, node.leaf.items[i]);
                if (distance < nearest_distance_squared) {
                    results.clear();
                    nearest_distance_squared = distance;
                    results.push_back(presence_.find(node.leaf.items[i]));
                } else if (distance == nearest_distance_squared) {
                    results.push_back(presence_.find(node.leaf.items[i]));
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
    storage_container                        storage_;
    storage_pool                             pool_;
    hash_table_type                          presence_;
};

}  // namespace internal

template <typename CoordinateType,
          typename StorageType,
          uint64_t rank = 2,
          uint16_t MAXIMUM_NODE_SIZE = 64>
class spatial_map : public internal::spatial_tree<StorageType, CoordinateType, MAXIMUM_NODE_SIZE> {
public:
    static_assert(!std::is_void_v<StorageType>, "For no storage type, use st::spatial_set");
    spatial_map()
        : internal::spatial_tree<StorageType, CoordinateType, rank, MAXIMUM_NODE_SIZE>() {}
    spatial_map(const bounding_box<CoordinateType, rank> &boundary)
        : internal::spatial_tree<StorageType, CoordinateType, rank, MAXIMUM_NODE_SIZE>(boundary) {}
};

template <typename CoordinateType = double, uint64_t rank = 2, uint16_t MAXIMUM_NODE_SIZE = 64>
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