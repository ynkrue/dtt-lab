#include "dtt/sim/particles.h"

#include <random>

namespace dtt::sim {

bool validate_particles(const Particles &particles) {
    const std::size_t n = particles.x.size();
    return particles.y.size() == n && particles.vx.size() == n && particles.vy.size() == n &&
           particles.mass.size() == n;
}

Particles create_rnd_particles(std::size_t n, unsigned seed, double min_x, double max_x,
                               double min_y, double max_y, double mass_min, double mass_max,
                               ParticleDistribution distribution, double velocity_std) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist_x(min_x, max_x);
    std::uniform_real_distribution<double> dist_y(min_y, max_y);
    std::uniform_real_distribution<double> dist_mass(mass_min, mass_max);
    std::normal_distribution<double> dist_pos_x(0.0, 0.25 * (max_x - min_x));
    std::normal_distribution<double> dist_pos_y(0.0, 0.25 * (max_y - min_y));
    std::normal_distribution<double> dist_vel(0.0, velocity_std);

    Particles particles;
    particles.x.resize(n);
    particles.y.resize(n);
    particles.vx.assign(n, 0.0);
    particles.vy.assign(n, 0.0);
    particles.mass.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        switch (distribution) {
        case ParticleDistribution::kNormal:
            particles.x[i] = dist_pos_x(rng);
            particles.y[i] = dist_pos_y(rng);
            break;
        case ParticleDistribution::kUniform:
        default:
            particles.x[i] = dist_x(rng);
            particles.y[i] = dist_y(rng);
            break;
        }
        particles.mass[i] = dist_mass(rng);
        if (velocity_std > 0.0) {
            particles.vx[i] = dist_vel(rng);
            particles.vy[i] = dist_vel(rng);
        }
    }

    return particles;
}

} // namespace dtt::sim
