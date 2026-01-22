#include "dtt/core/particles.h"
#include "dtt/tree/build.h"

#include <benchmark/benchmark.h>

#include <random>
#include <vector>

using dtt::core::ParticlesBuffer;
using dtt::tree::BuildParams;
using dtt::tree::Tree;

namespace {

ParticlesBuffer make_random_particles(std::size_t n, unsigned seed = 1234) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    ParticlesBuffer buf = ParticlesBuffer::make(n);
    for (std::size_t i = 0; i < n; ++i) {
        buf.x[i] = dist(rng);
        buf.y[i] = dist(rng);
        buf.mass[i] = 1.0;
        buf.vx[i] = buf.vy[i] = 0.0;
    }
    return buf;
}

} // namespace

static void BM_BuildQuadtree(benchmark::State &state) {
    const std::size_t n = static_cast<std::size_t>(state.range(0));
    ParticlesBuffer buf = make_random_particles(n);
    BuildParams params;
    for (auto _ : state) {
        Tree t = dtt::tree::build_quadtree(buf.const_view(), params);
        benchmark::DoNotOptimize(t);
    }
    state.SetItemsProcessed(state.iterations() * static_cast<long>(n));
}
BENCHMARK(BM_BuildQuadtree)->RangeMultiplier(4)->Range(1 << 11, 1 << 23);

BENCHMARK_MAIN();
