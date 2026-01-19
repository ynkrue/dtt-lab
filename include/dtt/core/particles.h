#pragma once

#include "dtt/core/memory.h"

#include <cstddef>

namespace dtt::core {
struct ParticlesView;
struct ConstParticlesView;

/**
 * @file particles.h
 * @brief Definitions related to particle data structures.
 * Defines owning structures for managing particle attributes such as position, velocity, and mass.
 * @author Yannik Rüfenacht
 */
struct ParticlesBuffer {
    // position and velocity
    dtt::core::Memory<double> x;
    dtt::core::Memory<double> y;
    dtt::core::Memory<double> vx;
    dtt::core::Memory<double> vy;
    dtt::core::Memory<double> mass;

    // number of particles
    std::size_t count{};

    static ParticlesBuffer make(std::size_t count,
                                dtt::core::MemoryType type = dtt::core::MemoryType::HOST,
                                std::size_t alignment = 64);

    ParticlesView view();
    ConstParticlesView const_view() const;

    // default ctors, disable copy, allow move
    ParticlesBuffer() = default;
    ParticlesBuffer(ParticlesBuffer &&other) noexcept = default;
    ParticlesBuffer &operator=(ParticlesBuffer &&other) noexcept = default;
    ParticlesBuffer(const ParticlesBuffer &other) = delete;
    ParticlesBuffer &operator=(const ParticlesBuffer &other) = delete;
    ~ParticlesBuffer() = default;
};

/**
 * @brief A view structure for read-only access to particle data.
 */
struct ConstParticlesView {
    const double *x{};
    const double *y{};
    const double *vx{};
    const double *vy{};
    const double *mass{};

    std::size_t count{0};
    bool valid() const;
};

/**
 * @brief A view structure for read-write access to particle data.
 */
struct ParticlesView {
    double *x{};
    double *y{};
    double *vx{};
    double *vy{};
    double *mass{};

    std::size_t count{0};
    bool valid() const;
};

} // namespace dtt::core
