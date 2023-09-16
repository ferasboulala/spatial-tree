#pragma once

#include <cassert>
#include <limits>

namespace st {
template <typename CoordinateType>
struct BoundingBox {
    CoordinateType  top_x, top_y, bottom_x, bottom_y;
    __always_inline BoundingBox()
        : top_x(std::numeric_limits<CoordinateType>::min()),
          top_y(std::numeric_limits<CoordinateType>::max()),
          bottom_x(std::numeric_limits<CoordinateType>::max()),
          bottom_y(std::numeric_limits<CoordinateType>::min()) {
        assert(top_x <= bottom_x);
        assert(top_y >= bottom_y);
    }
    __always_inline BoundingBox(CoordinateType top_x_,
                                CoordinateType top_y_,
                                CoordinateType bottom_x_,
                                CoordinateType bottom_y_)
        : top_x(top_x_), top_y(top_y_), bottom_x(bottom_x_), bottom_y(bottom_y_) {
        assert(top_x <= bottom_x);
        assert(top_y >= bottom_y);
    }
};

template <typename CoordinateType>
__always_inline bool is_within_interval(CoordinateType x, CoordinateType beg, CoordinateType end) {
    assert(end >= beg);
    return x >= beg && x <= end;
}

template <typename CoordinateType>
static __always_inline bool is_inside_bounding_box(CoordinateType                     x,
                                                   CoordinateType                     y,
                                                   const BoundingBox<CoordinateType> &bbox) {
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
static __always_inline bool bounding_boxes_overlap(const BoundingBox<CoordinateType> &lhs,
                                                   const BoundingBox<CoordinateType> &rhs) {
    return intervals_overlap(lhs.top_x, lhs.bottom_x, rhs.top_x, rhs.bottom_x) &&
           intervals_overlap(lhs.bottom_y, lhs.top_y, rhs.bottom_y, rhs.top_y);
}

}  // namespace st