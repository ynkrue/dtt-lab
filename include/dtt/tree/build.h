#pragma once

#include "dtt/core/particles.h"
#include "dtt/tree/bounding-box.h"
#include "dtt/tree/tree.h"

namespace dtt::tree {
/**
 * @file build.h
 * @brief Definitions related to tree building functions.
 * Defines functions for constructing spatial partitioning trees.
 * @author Yannik Rüfenacht
 */

/**
 * Parameters for building the tree.
 * @note For large n override max depth: example ceil(log4(n/max_leaf_size))
 */
struct BuildParams {
    std::size_t max_leaf_size = 16; // maximal particles per node for accuracy
    std::size_t max_depth = 20;     // maximum tree depth for performance
    double padding = 1e-6;
};

/// @brief Simple basline implmenetation without morton-code ordering
dtt::tree::Tree build_quadtree_base(const dtt::core::ConstParticlesView &particles,
                                    const BuildParams &params);

/// @brief Quadtree build with mortoncode ordering and optional OpenMP parallelization and SIMD
/// vectorization
dtt::tree::Tree build_quadtree(const dtt::core::ConstParticlesView &particles,
                               const BuildParams &params);

/// @brief OpenMP parallelized quadtree builder using producer-consumer pattern
dtt::tree::Tree build_quadtree_omp(const dtt::core::ConstParticlesView &particles,
                                   const BuildParams &params);

#ifdef DTT_ENABLE_CUDA
namespace cu {
/// @brief CUDA parallelized quadtree builder using local producer-consumer pattern
dtt::tree::Tree build_quadtree_cuda(const dtt::core::ConstParticlesView &particles,
                                    const BuildParams &params);
} // namespace cu
#endif
} // namespace dtt::tree