#pragma once

#include "dtt/sim/particles.h"
#include "dtt/tree/bounding_box.h"
#include "dtt/tree/morton.h"

#include <cstddef>
#include <vector>

namespace dtt::tree {

struct Node {
    BoundingBox bbox;
    int children[4]; // indices of the 4 children, -1 if no child
    int first;
    int count;
    double mass;
    double com_x, com_y;
    int level;
    bool is_leaf() const {
        return children[0] == -1 && children[1] == -1 && children[2] == -1 && children[3] == -1;
    }
};

struct Tree {
    std::vector<Node> nodes;
    std::vector<std::size_t> indices; // particles permutation ordered by leafs

    int root = -1;
    BoundingBox root_bbox;
};

struct BuildParams {
    int max_leaf_size = 16;
    int max_depth = 20;
    double padding = 1e-6;
};

Tree build_quadtree(const sim::Particles &particles, const BuildParams &params);

BoundingBox compute_root_bbox(const sim::Particles &particles, double padding);

struct TraversalStats {
    uint64_t popped = 0;
    uint64_t refined = 0;
    uint64_t accepted = 0;
    uint64_t pruned = 0;
    uint64_t leaf_leaf = 0;
};

struct MAC {
    double theta;
};

bool accept(const Node &a, const Node &b, const MAC &mac);

template <typename LeafFunc, typename AcceptFunc, typename PruneFunc>
void dual_tree_traversal(const Tree &tree, const MAC &mac, LeafFunc on_leaf_pair,
                         AcceptFunc on_accepted_pair, PruneFunc leaf,
                         TraversalStats *stats = nullptr) {
    (void)stats;
    
}

} // namespace dtt::tree