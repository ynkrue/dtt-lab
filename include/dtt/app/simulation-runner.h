#pragma once

#include "dtt/core/forces.h"
#include "dtt/core/particles.h"
#include "dtt/tree/force-engine.h"

namespace dtt::app {

/**
 * @file simulation-runner.h
 * @brief Definitions related to simulation runner.
 * Defines functions for running particle simulations using tree-based force calculations.
 * @author Yannik Rüfenacht
 */

void run_simulation(const dtt::core::ConstParticlesView &initial_particles,
                    dtt::core::ParticlesView &particles,
                    const dtt::tree::ForceEngineParams &force_params, const double time_step,
                    const int num_steps);

} // namespace dtt::app