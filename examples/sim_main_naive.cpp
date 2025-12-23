#include "dtt/io/vtk.h"
#include "dtt/sim/forces.h"
#include "dtt/sim/particles.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

int main(int argc, char** argv) {
    using dtt::sim::ForceParams;
    using dtt::sim::ForceField;
    using dtt::sim::Particles;
    using dtt::sim::compute_forces_naive;
    using dtt::sim::create_rnd_particles;
    using dtt::sim::euler_step;

    const std::size_t n = (argc > 1) ? static_cast<std::size_t>(std::stoul(argv[1])) : 512;
    const int steps = (argc > 2) ? std::stoi(argv[2]) : 20;
    const double dt = (argc > 3) ? std::stod(argv[3]) : 0.01;
    const std::string dist = (argc > 4) ? std::string(argv[4]) : "uniform";
    const double cutoff = (argc > 5) ? std::stod(argv[5]) : std::numeric_limits<double>::quiet_NaN();
    const double gravity = (argc > 6) ? std::stod(argv[6]) : 9.0;
    const std::string out_dir = (argc > 7) ? argv[7] : "output_naive";

    std::cout << "Running naive simulation with n=" << n << ", steps=" << steps
              << ", dt=" << dt << ", distribution=" << dist << ", gravity=" << gravity;
    if (!std::isnan(cutoff)) {
        std::cout << ", cutoff=" << cutoff;
    }
    std::cout << "\n";

    std::filesystem::create_directories(out_dir);

    // Scale initial domain with sqrt(n) to keep density roughly constant as n grows.
    const double half_extent = 0.75 * std::sqrt(static_cast<double>(n));
    auto distribution = dtt::sim::ParticleDistribution::kUniform;
    if (dist == "cluster") {
        distribution = dtt::sim::ParticleDistribution::kCluster;
    } else if (dist == "curl") {
        distribution = dtt::sim::ParticleDistribution::kCurl;
    }
    Particles particles =
        create_rnd_particles(n, /*seed=*/1234, -half_extent, half_extent, -half_extent,
                             half_extent, 1.0, 5.0, distribution, 0.02);
    std::optional<double> cutoff_opt = std::isnan(cutoff) ? std::nullopt : std::make_optional(cutoff);
    const ForceParams params{.softening = 1e-3, .cutoff = cutoff_opt, .gravity = gravity};
    const dtt::sim::Boundary bounds{.xmin = -half_extent,
                                    .xmax = half_extent,
                                    .ymin = -half_extent,
                                    .ymax = half_extent,
                                    .restitution = 0.9};

    std::vector<std::pair<double, std::string>> frames;
    frames.reserve(static_cast<std::size_t>(steps));

    const auto start = std::chrono::steady_clock::now();
    const int progress_stride = std::max(1, steps / 10);  // ~10% checkpoints

    // main simulation loop
    for (int step = 0; step < steps; ++step) {
        ForceField forces = compute_forces_naive(particles, params);
        euler_step(particles, forces, dt, &bounds);

        std::ostringstream oss;
        oss << out_dir << "/frame_" << std::setw(4) << std::setfill('0') << step << ".vtp";
        const std::string vtp_path = oss.str();
        if (!dtt::io::write_vtp_polydata(particles, vtp_path, &forces)) {
            std::cerr << "Failed to write " << vtp_path << "\n";
            return 1;
        }
        frames.emplace_back(step * dt, std::filesystem::path(vtp_path).filename().string());

        if ((step + 1) % progress_stride == 0 || step + 1 == steps) {
            const auto now = std::chrono::steady_clock::now();
            const double elapsed_s =
                std::chrono::duration_cast<std::chrono::duration<double>>(now - start).count();
            const double percent = 100.0 * (step + 1) / static_cast<double>(steps);
            const double avg_per_step = elapsed_s / static_cast<double>(step + 1);
            const double eta = avg_per_step * (steps - (step + 1));

            const int bar_width = 20;
            const int filled = static_cast<int>(percent / 100.0 * bar_width + 0.5);
            std::string bar(filled, '#');
            bar.resize(bar_width, '.');

            std::cout << "[" << bar << "] step " << (step + 1) << "/" << steps << " ("
                      << "elapsed " << std::setprecision(3) << elapsed_s << "s, "
                      << "avg " << std::setprecision(4) << avg_per_step << "s/step, "
                      << "ETA " << std::setprecision(3) << eta << "s\n";
        }
    }

    // write to file
    const std::string pvd_path = out_dir + "/frames.pvd";
    if (!dtt::io::write_pvd_collection(pvd_path, frames)) {
        std::cerr << "Failed to write " << pvd_path << "\n";
        return 1;
    }

    std::cout << "Wrote " << steps << " frames and collection: " << pvd_path << "\n";
    std::cout << "Load frames.pvd in ParaView and press Play.\n";
    return 0;
}
