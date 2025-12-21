#include "dtt/sim/forces.h"
#include "dtt/sim/particles.h"

#include <benchmark/benchmark.h>

namespace {

using dtt::sim::ForceParams;
using dtt::sim::Particles;
using dtt::sim::compute_forces_naive;
using dtt::sim::create_rnd_particles;

Particles make_particles(std::size_t n) {
    // Deterministic cloud in a unit square, unit masses.
    return create_rnd_particles(n, /*seed=*/1234, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0);
}

}  // namespace

static void BM_NaiveAllPairs(benchmark::State& state) {
    const std::size_t n = static_cast<std::size_t>(state.range(0));
    const Particles particles = make_particles(n);
    const ForceParams params{.softening = 1e-3, .cutoff = std::nullopt};

    for (auto _ : state) {
        auto forces = compute_forces_naive(particles, params);
        benchmark::DoNotOptimize(forces);
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * static_cast<long>(n) * static_cast<long>(n));
}
BENCHMARK(BM_NaiveAllPairs)->RangeMultiplier(2)->Range(256, 2048);

static void BM_NaiveWithCutoff(benchmark::State& state) {
    const std::size_t n = static_cast<std::size_t>(state.range(0));
    const Particles particles = make_particles(n);
    const ForceParams params{.softening = 1e-3, .cutoff = 0.1};

    for (auto _ : state) {
        auto forces = compute_forces_naive(particles, params);
        benchmark::DoNotOptimize(forces);
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * static_cast<long>(n) * static_cast<long>(n));
}
BENCHMARK(BM_NaiveWithCutoff)->RangeMultiplier(2)->Range(256, 2048);

BENCHMARK_MAIN();
