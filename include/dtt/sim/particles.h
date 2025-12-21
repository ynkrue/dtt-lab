#pragma once

#include <vector>

namespace dtt::sim {

struct Particles {
    std::vector<double> x, y;
    std::vector<double> vx, vy;
    std::vector<double> mass;
};

enum class ParticleDistribution {
    kUniform,
    kNormal
};

// ensure particle data has same size
bool validate_particles(const Particles &particles);

// resize and random init helper
// Positions follow distribution; velocity_std > 0 seeds velocities from N(0, std).
Particles create_rnd_particles(std::size_t n, unsigned seed, double min_x, double max_x,
                               double min_y, double max_y, double mass_min = 1.0,
                               double mass_max = 1.0,
                               ParticleDistribution distribution = ParticleDistribution::kUniform,
                               double velocity_std = 0.0);
} // namespace dtt::sim
