#include "dtt/sim/forces.h"
#include "dtt/tree/quadtree.h"

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
            const double scale = params.gravity * p.mass[i] * p.mass[j] * inv_r3;
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

ForceField compute_forces_tree(const Particles &particles, const tree::Tree &tree,
                               ForceParams params, tree::MAC mac) {
    const std::size_t n = particles.x.size();
    ForceField forces(n, {0.0, 0.0});

    // cutoff nodes with large distance
    auto prune = [&](const tree::Node &a, const tree::Node &b) {
        if (!params.cutoff) return false;
        const double cutoff_sq= (*params.cutoff) * (*params.cutoff);
        return a.bbox.distance_lower_bound_sq(b.bbox) > cutoff_sq;
    };

    // force helper function
    auto apply_force = [&](std::size_t i, std::size_t j, double dx, double dy) {
        const double r2 = dx * dx + dy * dy + params.softening * params.softening;
        const double inv_r = 1.0 / std::sqrt(r2);
        const double inv_r3 = inv_r * inv_r * inv_r;
        const double scale = params.gravity * particles.mass[i] * particles.mass[j] * inv_r3;
        const double fx = scale * dx;
        const double fy = scale * dy;
        forces[i][0] += fx; forces[i][1] += fy;
        forces[j][0] -= fx; forces[j][1] -= fy;
    };


    auto on_leaf_pair = [&](const tree::Node& a, const tree::Node& b) {
        const int a_first = a.first, a_count = a.count;
        const int b_first = b.first, b_count = b.count;
        if (&a == &b) {
            // self-interaction, loop over unique pairs only
            for (int i = 0; i < a_count; ++i) {
                std::size_t idx_i = tree.indices[a_first + i];
                for (int j = i + 1; j < a_count; ++j) {
                    std::size_t idx_j = tree.indices[a_first + j];
                    apply_force(idx_i, idx_j, particles.x[idx_j] - particles.x[idx_i],
                                particles.y[idx_j] - particles.y[idx_i]);
                }
            }
        } else {
            // full cross interaction
            for (int i = 0; i < a_count; ++i) {
                std::size_t idx_i = tree.indices[a_first + i];
                for (int j = 0; j < b_count; ++j) {
                    std::size_t idx_j = tree.indices[b_first + j];
                    apply_force(idx_i, idx_j, particles.x[idx_j] - particles.x[idx_i],
                                particles.y[idx_j] - particles.y[idx_i]);
                }
            }
        }
    };

    // approximate cluster-cluster interaction (point at COM)
    auto on_accepted_pair = [&](const tree::Node& a, const tree::Node& b) {
        const double dx = b.com_x - a.com_x;
        const double dy = b.com_y - a.com_y;
        const double r2 = dx * dx + dy * dy + params.softening * params.softening;
        const double inv_r = 1.0 / std::sqrt(r2);
        const double inv_r3 = inv_r * inv_r * inv_r;
        const double scale = params.gravity * a.mass * b.mass * inv_r3;
        const double fx = scale * dx;
        const double fy = scale * dy;

        // Distribute force to all particles in the nodes proportionally to their mass
        for (int k = 0; k < a.count; ++k) {
            std::size_t idx = tree.indices[a.first + k];
            const double m_frac = particles.mass[idx] / a.mass;
            forces[idx][0] += m_frac * fx;
            forces[idx][1] += m_frac * fy;
        }
        for (int k = 0; k < b.count; ++k) {
            std::size_t idx = tree.indices[b.first + k];
            const double m_frac = particles.mass[idx] / b.mass;
            forces[idx][0] -= m_frac * fx;
            forces[idx][1] -= m_frac * fy;
        }
    };

    tree::TraversalStats stats;
    tree::dual_tree_traversal(tree, mac, on_leaf_pair, on_accepted_pair, prune, &stats);
    return forces;
}

void euler_step(Particles &particles, const ForceField &forces, double dt, const Boundary *bounds) {
    const std::size_t n = particles.x.size();
    for (std::size_t i = 0; i < n; ++i) {
        const double ax = forces[i][0] / particles.mass[i];
        const double ay = forces[i][1] / particles.mass[i];

        particles.vx[i] += ax * dt;
        particles.vy[i] += ay * dt;

        particles.x[i] += particles.vx[i] * dt;
        particles.y[i] += particles.vy[i] * dt;

        if (bounds) {
            if (particles.x[i] < bounds->xmin) {
                particles.x[i] = bounds->xmin;
                particles.vx[i] = -particles.vx[i] * bounds->restitution;
            } else if (particles.x[i] > bounds->xmax) {
                particles.x[i] = bounds->xmax;
                particles.vx[i] = -particles.vx[i] * bounds->restitution;
            }
            if (particles.y[i] < bounds->ymin) {
                particles.y[i] = bounds->ymin;
                particles.vy[i] = -particles.vy[i] * bounds->restitution;
            } else if (particles.y[i] > bounds->ymax) {
                particles.y[i] = bounds->ymax;
                particles.vy[i] = -particles.vy[i] * bounds->restitution;
            }
        }
    }
}

} // namespace dtt::sim
