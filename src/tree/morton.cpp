#include "dtt/tree/morton.h"

namespace dtt::tree {

void normalize_point(const BoundingBox &root, double x, double y, double &nx, double &ny) {
    nx = (x - root.min_x) / (root.max_x - root.min_x);
    ny = (y - root.min_y) / (root.max_y - root.min_y);
}

uint64_t morton_from_normalized(double nx, double ny, uint32_t levels) {
    uint32_t max_coord = 1u << levels;
    uint32_t ix = static_cast<uint32_t>(
        std::min(std::max(nx * max_coord, 0.0), static_cast<double>(max_coord - 1)));
    uint32_t iy = static_cast<uint32_t>(
        std::min(std::max(ny * max_coord, 0.0), static_cast<double>(max_coord - 1)));
    return morton_encode_2d(ix, iy);
}

} // namespace dtt::tree