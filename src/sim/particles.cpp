#include "dtt/sim/particles.h"

#include <algorithm>
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
    std::normal_distribution<double> dist_vel(0.0, velocity_std);
    const double extent = std::max(max_x - min_x, max_y - min_y);
    const double cluster_sigma = 0.1 * extent; // controls cluster tightness

    Particles particles;
    particles.x.resize(n);
    particles.y.resize(n);
    particles.vx.assign(n, 0.0);
    particles.vy.assign(n, 0.0);
    particles.mass.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        switch (distribution) {
        case ParticleDistribution::kCluster: {
            // Cluster centered in the middle of the domain with Gaussian falloff.
            const double cx = 0.5 * (min_x + max_x);
            const double cy = 0.5 * (min_y + max_y);
            std::normal_distribution<double> cluster_x(cx, cluster_sigma);
            std::normal_distribution<double> cluster_y(cy, cluster_sigma);
            particles.x[i] = std::clamp(cluster_x(rng), min_x, max_x);
            particles.y[i] = std::clamp(cluster_y(rng), min_y, max_y);
            break;
        }
        case ParticleDistribution::kCurl: {
            // Same spatial distribution as cluster, but add tangential velocity to swirl.
            const double cx = 0.5 * (min_x + max_x);
            const double cy = 0.5 * (min_y + max_y);
            std::normal_distribution<double> cluster_x(cx, cluster_sigma);
            std::normal_distribution<double> cluster_y(cy, cluster_sigma);
            particles.x[i] = std::clamp(cluster_x(rng), min_x, max_x);
            particles.y[i] = std::clamp(cluster_y(rng), min_y, max_y);
            const double dx = particles.x[i] - cx;
            const double dy = particles.y[i] - cy;
            const double r = std::hypot(dx, dy);
            if (r > 1e-12) {
                const double tx = -dy / r;
                const double ty = dx / r;
                const double omega = 1000.0; // stronger spin
                const double v = omega * r;
                particles.vx[i] = v * tx;
                particles.vy[i] = v * ty;
            }
            break;
        }
        case ParticleDistribution::kUniform:
        default:
            particles.x[i] = dist_x(rng);
            particles.y[i] = dist_y(rng);
            break;
        }
        particles.mass[i] = dist_mass(rng);
        if (velocity_std > 0.0 && distribution != ParticleDistribution::kCurl) {
            particles.vx[i] = dist_vel(rng);
            particles.vy[i] = dist_vel(rng);
        }
    }

    return particles;
}

} // namespace dtt::sim
