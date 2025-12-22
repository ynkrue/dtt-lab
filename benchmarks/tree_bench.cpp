#include "dtt/sim/forces.h"
#include "dtt/sim/particles.h"
#include "dtt/tree/quadtree.h"

#include <benchmark/benchmark.h>

namespace {

using dtt::sim::ForceField;
using dtt::sim::ForceParams;
using dtt::sim::Particles;
using dtt::sim::compute_forces_tree;
using dtt::sim::create_rnd_particles;
using dtt::tree::BuildParams;
using dtt::tree::MAC;
using dtt::tree::Tree;
using dtt::tree::build_quadtree;

Particles make_particles(std::size_t n) {
    const double half_extent = 0.75 * std::sqrt(static_cast<double>(n));
    return create_rnd_particles(n, /*seed=*/1234, -half_extent, half_extent, -half_extent,
                               half_extent, 5.0, 50.0,
                               dtt::sim::ParticleDistribution::kUniform, 0.02);
}

}  // namespace

static void BM_DualTree(benchmark::State& state) {
    const std::size_t n = static_cast<std::size_t>(state.range(0));
    Particles particles = make_particles(n);
    ForceParams params{.softening = 1e-3, .cutoff = std::nullopt, .gravity = 1.0};
    BuildParams build_params{.max_leaf_size = 16, .max_depth = 20, .padding = 1e-6};
    MAC mac{.theta = 0.7};

    for (auto _ : state) {
        Tree tree = build_quadtree(particles, build_params);
        ForceField forces = compute_forces_tree(particles, tree, params, mac);
        benchmark::DoNotOptimize(forces);
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * static_cast<long>(n));
}
BENCHMARK(BM_DualTree)->RangeMultiplier(2)->Range(2048, 2 << 17);
