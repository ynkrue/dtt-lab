#include "dtt/tree/bounding-box.h"

namespace dtt::tree {

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

bool BoundingBox::valid() const {
    return (max_x >= min_x) && (max_y >= min_y);
}

void BoundingBox::clamp(double &x, double &y) const {
    if (x < min_x)
        x = min_x;
    else if (x > max_x)
        x = max_x;
    if (y < min_y)
        y = min_y;
    else if (y > max_y)
        y = max_y;
}

void BoundingBox::make_square() {
    double w = width();
    double h = height();
    if (w > h) {
        double diff = w - h;
        min_y -= diff / 2.0;
        max_y += diff / 2.0;
    } else {
        double diff = h - w;
        min_x -= diff / 2.0;
        max_x += diff / 2.0;
    }
}

double BoundingBox::distance_sq_lower_bound(const BoundingBox &other) const {
    double dx = 0.0;
    if (max_x < other.min_x)
        dx = other.min_x - max_x;
    else if (other.max_x < min_x)
        dx = min_x - other.max_x;

    double dy = 0.0;
    if (max_y < other.min_y)
        dy = other.min_y - max_y;
    else if (other.max_y < min_y)
        dy = min_y - other.max_y;

    return dx * dx + dy * dy;
}

} // namespace dtt::tree