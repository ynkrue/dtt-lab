#pragma once

#include "dtt/sim/particles.h"
#include "dtt/tree/quadtree.h"

#include <array>
#include <optional>
#include <vector>

namespace dtt::sim {

using Force2D = std::array<double, 2>;
using ForceField = std::vector<Force2D>;

struct ForceParams {
    double softening;
    std::optional<double> cutoff;
    double gravity = 5.0; // scaling factor for interaction strength
};

// naive O(N^2) all-pairs force computation
ForceField compute_forces_naive(const Particles &particles, ForceParams params);

// tree-based O(N log N) force computation
ForceField compute_forces_tree(const Particles &particles, const tree::Tree &tree,
                               ForceParams params, tree::MAC mac);

// Optional reflecting boundary for integration.
struct Boundary {
    double xmin, xmax;
    double ymin, ymax;
    double restitution = 0.8; // 1.0 = perfect reflection, <1 dampens velocity on bounce
};

// simple euler step: x += vx * dt, vx += ax*dt, etc. Bounds reflect if provided.
void euler_step(Particles &particles, const ForceField &forces, double dt,
                const Boundary *bounds = nullptr);

} // namespace dtt::sim
