#include "dtt/tree/bounding_box.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace dtt::tree {

void BoundingBox::expand(double x, double y) {
    min_x = std::min(min_x, x);
    min_y = std::min(min_y, y);
    max_x = std::max(max_x, x);
    max_y = std::max(max_y, y);
}

void BoundingBox::pad(double eps) {
    min_x -= eps;
    min_y -= eps;
    max_x += eps;
    max_y += eps;
}

double BoundingBox::width() const {
    return max_x - min_x;
}

double BoundingBox::height() const {
    return max_y - min_y;
}

double BoundingBox::max_extent() const {
    return std::max(width(), height());
}

bool BoundingBox::valid() const {
    return min_x <= max_x && min_y <= max_y;
}

void BoundingBox::clamp(double &x, double &y) const {
    x = std::max(min_x, std::min(max_x, x));
    y = std::max(min_y, std::min(max_y, y));
}

void BoundingBox::make_square() {
    double w = width();
    double h = height();
    double max_side = std::max(w, h);
    double cx = 0.5 * (min_x + max_x);
    double cy = 0.5 * (min_y + max_y);
    min_x = cx - 0.5 * max_side;
    max_x = cx + 0.5 * max_side;
    min_y = cy - 0.5 * max_side;
    max_y = cy + 0.5 * max_side;
}

double BoundingBox::distance_lower_bound_sq(const BoundingBox &other) const {
    double dx = 0.0;
    if (max_x < other.min_x) {
        dx = other.min_x - max_x;
    } else if (other.max_x < min_x) {
        dx = min_x - other.max_x;
    }

    double dy = 0.0;
    if (max_y < other.min_y) {
        dy = other.min_y - max_y;
    } else if (other.max_y < min_y) {
        dy = min_y - other.max_y;
    }

    return dx * dx + dy * dy;
}

} // namespace dtt::tree