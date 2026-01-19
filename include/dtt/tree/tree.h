#pragma once

#include "dtt/core/memory.h"
#include "dtt/tree/bounding-box.h"

#include <cstddef>

namespace dtt::tree {
/**
 * @file tree.h
 * @brief Definitions related to tree data structures.
 * Defines structures and functions for managing spatial partitioning trees.
 * @author Yannik Rüfenacht
 */

struct Node {
    // indices of the 4 children, -1 if no child
    std::ptrdiff_t children[4] = {-1, -1, -1, -1};
    BoundingBox box{0, 0, 0, 0};
    std::size_t depth{0};

    std::size_t first{0}; // start index of particles index of this node
    std::size_t count{0}; // #particles in this node

    double mass{0.0};
    double com_x{0.0}, com_y{0.0}; // center of mass

    bool is_leaf() const {
        return children[0] == -1 && children[1] == -1 && children[2] == -1 && children[3] == -1;
    }
};

struct Tree {
    dtt::core::Memory<Node> nodes;
    dtt::core::Memory<std::size_t> indices; // particles permutation ordered by leafs
    std::size_t node_count{0};
    std::size_t index_count{0};

    std::size_t root = -1;
    BoundingBox root_box;
};
} // namespace dtt::tree