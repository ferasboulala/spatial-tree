#pragma once

#include <cmath>

namespace {
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
    __always_inline T safe_absdiff(T x, T y) {
        const bool gt = x > y;
        return gt * (x - y) + (1 - gt) * (y - x);
    }
    __always_inline T operator()(T x, T y) { return safe_absdiff(x, y); }
};

template <typename T>
struct UnsafeAbsDiff {
    __always_inline T unsafe_absdiff(T x, T y) { return x - y; }
    __always_inline T operator()(T x, T y) { return unsafe_absdiff(x, y); }
};

}  // namespace

namespace st {

template <typename T, typename... Args>
__always_inline T euclidean_distance_squared(T x, Args... xs) {
    if constexpr (std::is_floating_point_v<T> || std::is_signed_v<T>) {
        return euclidean_distance_squared_impl(UnsafeAbsDiff<T>(), x, xs...);
    } else {
        return euclidean_distance_squared_impl(SafeAbsDiff<T>(), x, xs...);
    }
}

template <typename T>
__always_inline double euclidean_distance(T x1, T y1, T x2, T y2) {
    return std::sqrt(euclidean_distance_squared(x1, y1, x2, y2));
}

}  // namespace st