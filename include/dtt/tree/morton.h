#pragma once

#include "dtt/tree/bounding_box.h"

#include <cstdint>

namespace dtt::tree {

inline uint64_t morton_encode_2d(uint32_t x, uint32_t y) {
    uint64_t answer = 0;
    for (uint64_t i = 0; i < (sizeof(uint32_t) * 8); ++i) {
        answer |= ((x & (1ULL << i)) << i) | ((y & (1ULL << i)) << (i + 1));
    }
    return answer;
}

// for debugging purposes
inline void morton_decode_2d(uint64_t code, uint32_t &x, uint32_t &y) {
    x = 0;
    y = 0;
    for (uint64_t i = 0; i < (sizeof(uint32_t) * 8); ++i) {
        x |= (code & (1ULL << (2 * i))) >> i;
        y |= (code & (1ULL << (2 * i + 1))) >> (i + 1);
    }
}

void normalize_point(const BoundingBox &root, double x, double y, double &nx, double &ny);

uint64_t morton_from_normalized(double nx, double ny, uint32_t levels);

} // namespace dtt::tree