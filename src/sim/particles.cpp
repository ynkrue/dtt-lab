#include "dtt/sim/particles.h"

#include <random>

namespace dtt::sim {

bool validate_particles(const Particles &particles) {
    const std::size_t n = particles.x.size();
    return particles.y.size() == n && particles.vx.size() == n && particles.vy.size() == n &&
           particles.mass.size() == n;
}

Particles create_rnd_particles(std::size_t n, unsigned seed, double min_x, double max_x,
                               double min_y, double max_y, double mass_min, double mass_max) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist_x(min_x, max_x);
    std::uniform_real_distribution<double> dist_y(min_y, max_y);
    std::uniform_real_distribution<double> dist_mass(mass_min, mass_max);

    Particles particles;
    particles.x.resize(n);
    particles.y.resize(n);
    particles.vx.assign(n, 0.0);
    particles.vy.assign(n, 0.0);
    particles.mass.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        particles.x[i] = dist_x(rng);
        particles.y[i] = dist_y(rng);
        particles.mass[i] = dist_mass(rng);
    }

    return particles;
}

} // namespace dtt::sim
