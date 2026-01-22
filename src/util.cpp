#include "dtt/util.h"
#include "dtt/core/memory.h"

#include <cstddef>
#include <cstring>
#include <limits>
#include <tuple>
#include <vector>

namespace dtt::util {

std::tuple<double, double, double, double> minmax_box(const double *xs, const double *ys,
                                                      std::size_t count) {
    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();

#pragma omp parallel for reduction(min : min_x, min_y) reduction(max : max_x, max_y)
    for (std::size_t i = 0; i < count; ++i) {
        min_x = std::min(min_x, xs[i]);
        max_x = std::max(max_x, xs[i]);
        min_y = std::min(min_y, ys[i]);
        max_y = std::max(max_y, ys[i]);
    }
    return {min_x, max_x, min_y, max_y};
}

void morton_encode(uint64_t *keys, std::size_t *values, const double *xs, const double *ys,
                   std::size_t count, double min_x, double min_y, double max_x, double max_y) {
    const double scale_x = 1.0 / (max_x - min_x);
    const double scale_y = 1.0 / (max_y - min_y);

#pragma omp parallel for if (count >= 1000)
    for (std::size_t i = 0; i < count; ++i) {
        double nx = (xs[i] - min_x) * scale_x;
        double ny = (ys[i] - min_y) * scale_y;
        nx = std::min(std::max(nx, 0.0), 1.0);
        ny = std::min(std::max(ny, 0.0), 1.0);
        uint32_t ix = static_cast<uint32_t>(nx * ((1u << 21) - 1));
        uint32_t iy = static_cast<uint32_t>(ny * ((1u << 21) - 1));
        keys[i] = morton_encode_2d(ix, iy);
        values[i] = i;
    }
}

void radix_sort(uint64_t *keys, std::size_t *values, std::size_t count) {
    dtt::core::Memory<uint64_t> temp_keys_mem = dtt::core::Memory<uint64_t>::allocate(count);
    dtt::core::Memory<std::size_t> temp_values_mem =
        dtt::core::Memory<std::size_t>::allocate(count);
    uint64_t *keys_cur = keys;
    uint64_t *keys_tmp = temp_keys_mem.data();
    std::size_t *vals_cur = values;
    std::size_t *vals_tmp = temp_values_mem.data();

    const int bits_per_pass = 16;
    const int num_passes = (64 + bits_per_pass - 1) / bits_per_pass;
    const int bucket_size = 1 << bits_per_pass;

    for (int pass = 0; pass < num_passes; ++pass) {
        std::vector<int> bucket_count(bucket_size, 0);
        int shift = pass * bits_per_pass;

// counting occurrences
#pragma omp parallel for if (count >= 1000)
        for (std::size_t i = 0; i < count; ++i) {
            uint64_t key = (keys[i] >> shift) & (bucket_size - 1);
            bucket_count[key]++;
        }

        // compute prefix sums
        for (int i = 1; i < bucket_size; ++i) {
            bucket_count[i] += bucket_count[i - 1];
        }

// build output array
#pragma omp parallel for if (count >= 1000)
        for (std::size_t i = count; i-- > 0;) {
            uint64_t key = (keys_cur[i] >> shift) & (bucket_size - 1);
            const int pos = --bucket_count[key];
            keys_tmp[pos] = keys_cur[i];
            vals_tmp[pos] = vals_cur[i];
        }

        // swap buffers for next pass
        std::swap(keys_cur, keys_tmp);
        std::swap(vals_cur, vals_tmp);
    }

    // copy back to the original pointers (after an even number of swaps)
    if (keys_cur != keys) {
        std::memcpy(keys, keys_cur, count * sizeof(uint64_t));
        std::memcpy(values, vals_cur, count * sizeof(std::size_t));
    }
}

} // namespace dtt::util
