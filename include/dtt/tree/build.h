#pragma once

#include "dtt/tree/bounding-box.h"
#include "dtt/tree/tree.h"
#include "dtt/core/particles.h"

namespace dtt::tree
{
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
        int max_leaf_size = 16; // maximal particles per node for accuracy
        int max_depth = 20; // maximum tree depth for performance
        double padding = 1e-6;
    };

    dtt::tree::Tree build_quadtree(const dtt::core::ConstParticlesView &particles, const BuildParams &params);
}