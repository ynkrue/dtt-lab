#pragma once

#include "dtt/core/memory.h"

#include <cstddef>

namespace dtt::core {

struct ForcesView;
struct ConstForcesView;

/**
 * @file forces.h
 * @brief Definitions related to forces and forcefield data structures.
 * Defines owning structures for managing forces acting on particles.
 * @author Yannik Rüfenacht
 */
struct ForcesBuffer {
    dtt::core::Memory<double> fx;
    dtt::core::Memory<double> fy;

    std::size_t count{};

    static ForcesBuffer make(std::size_t count,
                             dtt::core::MemoryType type = dtt::core::MemoryType::HOST,
                             std::size_t alignment = 64);

    ForcesesView view();
    ConstForcesView const_view() const;

    // default ctors, disable copy, allow move
    ForcesBuffer() = default;
    ForcesBuffer(ForcesBuffer &&other) noexcept = default;
    ForcesBuffer &operator=(ForcesBuffer &&other) noexcept = default;
    ForcesBuffer(const ForcesBuffer &other) = delete;
    ForcesBuffer &operator=(const ForcesBuffer &other) = delete;
    ~ForcesBuffer() = default;
};

/**
 * @brief A view structure for read-only access to force data.
 */
struct ConstForcesView {
    const double *fx{};
    const double *fy{};

    std::size_t count{0};
    bool valid() const;
};

/**
 * @brief A view structure for read-write access to force data.
 */
struct ForcesView {
    double *fx{};
    double *fy{};

    std::size_t count{0};
    bool valid() const;
};

} // namespace dtt::core