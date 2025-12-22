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

BoundingBox compute_root_bbox(const sim::Particles &particles, double padding);

Tree build_quadtree(const sim::Particles &particles, const BuildParams &params);


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
                         AcceptFunc on_accepted_pair, PruneFunc prune,
                         TraversalStats *stats = nullptr) {
    if (tree.root == -1)
        return;
    std::vector<std::pair<int, int>> stack;
    stack.reserve(64);
    stack.emplace_back(tree.root, tree.root);

    auto push_if_valid = [&](int ia, int ib) {
        if (ia != -1 && ib != -1) {
            stack.emplace_back(ia, ib);
        }
    };

    while (!stack.empty()) {
        const auto [ia, ib] = stack.back();
        stack.pop_back();
        if (stats)
            stats->popped++;

        const Node &a = tree.nodes[ia];
        const Node &b = tree.nodes[ib];

        // prune with cutoff
        if (prune(a, b)) {
            if (stats)
                stats->pruned++;
            continue;
        }

        const bool a_leaf = a.is_leaf();
        const bool b_leaf = b.is_leaf();

        // cluster-cluster interaction
        if (accept(a, b, mac)) {
            on_accepted_pair(a, b);
            if (stats) stats->accepted++;
            continue;
        }

        // leaf-leaf interaction
        if (a_leaf && b_leaf) {
            on_leaf_pair(a, b);
            if (stats) stats->leaf_leaf++;
            continue;
        }

        // refine
        const bool refine_a = !a_leaf && (b_leaf || a.bbox.max_extent() >= b.bbox.max_extent() ||
                                          a.count >= b.count);
        if (ia == ib) {
            // self-interaction, generate unique pairs only
            for (int i = 0; i < 4; ++i) {
                const int ca = a.children[i];
                if (ca == -1) continue;
                for (int j = i; j < 4; ++j) {
                    const int cb = a.children[j];
                    if (cb == -1) continue;
                    push_if_valid(ca, cb);
                }
            }
        } else if (refine_a) {
            // refine A
            for (int i = 0; i < 4; ++i) {
                const int ca = a.children[i];
                push_if_valid(ca, ib);
            }
        } else {
            // refine B
            for (int i = 0; i < 4; ++i) {
                const int cb = b.children[i];
                push_if_valid(ia, cb);
            }
        }
        if (stats) stats->refined++;
    }
}

} // namespace dtt::tree
