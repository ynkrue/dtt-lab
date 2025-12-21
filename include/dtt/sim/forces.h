#pragma once

#include "dtt/sim/particles.h"

#include <array>
#include <optional>
#include <vector>

namespace dtt::sim {

using Force2D = std::array<double, 2>;
using ForceField = std::vector<Force2D>;

struct ForceParams {
    double softening;
    std::optional<double> cutoff;
};

// naive O(N^2) all-pairs force computation
ForceField compute_forces_naive(const Particles &particles, ForceParams params);

// simple euler step: x += vx * dt, vx += ax*dt, etc.
void euler_step(Particles &particles, const ForceField &forces, double dt);

} // namespace dtt::sim