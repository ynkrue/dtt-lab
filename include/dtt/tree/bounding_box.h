#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

namespace dtt::tree {

struct BoundingBox {
    double min_x;
    double min_y;
    double max_x;
    double max_y;

    void expand(double x, double y);

    void pad(double eps);

    double width() const;
    double height() const;
    double max_extent() const;

    bool valid() const;

    void clamp(double &x, double &y) const;
    void make_square();
    double distance_lower_bound_sq(const BoundingBox &other) const;
};

} // namespace dtt::tree