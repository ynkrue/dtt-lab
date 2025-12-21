#include "dtt/sim/forces.h"

#include <cmath>
#include <limits>

namespace dtt::sim {

ForceField compute_forces_naive(const Particles &p, ForceParams params) {
    const std::size_t n = p.x.size();
    ForceField forces(n, {0.0, 0.0});
    const double eps2 = params.softening * params.softening;
    const double cutoff2 = params.cutoff ? (*params.cutoff) * (*params.cutoff)
                                         : std::numeric_limits<double>::infinity();

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const double dx = p.x[j] - p.x[i];
            const double dy = p.y[j] - p.y[i];
            const double r2 = dx * dx + dy * dy + eps2;
            if (r2 > cutoff2)
                continue;

            const double inv_r = 1.0 / std::sqrt(r2);
            const double inv_r3 = inv_r * inv_r * inv_r;
            const double scale = p.mass[i] * p.mass[j] * inv_r3;
            const double fx = scale * dx;
            const double fy = scale * dy;

            forces[i][0] += fx;
            forces[i][1] += fy;
            forces[j][0] -= fx;
            forces[j][1] -= fy;
        }
    }
    return forces;
}

void euler_step(Particles &particles, const ForceField &forces, double dt) {
    const std::size_t n = particles.x.size();
    for (std::size_t i = 0; i < n; ++i) {
        const double ax = forces[i][0] / particles.mass[i];
        const double ay = forces[i][1] / particles.mass[i];

        particles.vx[i] += ax * dt;
        particles.vy[i] += ay * dt;

        particles.x[i] += particles.vx[i] * dt;
        particles.y[i] += particles.vy[i] * dt;
    }
}

} // namespace dtt::sim