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

void radix_sort_morton(std::vector<std::pair<uint64_t, std::size_t>> &keys) {
    const std::size_t n = keys.size();
    std::vector<std::pair<uint64_t, std::size_t>> temp(n);
    const int bits_per_pass = 16;
    const int num_passes = (64 + bits_per_pass - 1) / bits_per_pass;
    const int bucket_size = 1 << bits_per_pass;

    for (int pass = 0; pass < num_passes; ++pass) {
        std::vector<int> count(bucket_size, 0);
        int shift = pass * bits_per_pass;

        // counting occurrences
        for (std::size_t i = 0; i < n; ++i) {
            uint64_t key = (keys[i].first >> shift) & (bucket_size - 1);
            count[key]++;
        }

        // compute prefix sums
        for (int i = 1; i < bucket_size; ++i) {
            count[i] += count[i - 1];
        }

        // build output array
        for (std::size_t i = n; i-- > 0;) {
            uint64_t key = (keys[i].first >> shift) & (bucket_size - 1);
            temp[--count[key]] = keys[i];
        }

        // copy back to original array
        keys.swap(temp);
    }
}

} // namespace dtt::tree