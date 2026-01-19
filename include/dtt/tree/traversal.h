#pragma once

#include "dtt/tree/bounding-box.h"
#include "dtt/tree/tree.h"

#include <algorithm>
#include <cstddef>

namespace dtt::tree {
/**
 * @file traversal.h
 * @brief Definitions related to tree traversal functions.
 * Defines functions for traversing spatial partitioning trees.
 * @author Yannik Rüfenacht
 */

/**
 * @brief Determines whether two nodes should be accepted for interaction based on the opening angle
 * theta.
 * @param a First node.
 * @param b Second node.
 * @param theta Opening angle parameter for the multipole acceptance criterion.
 * @return True if the nodes satisfy the acceptance criterion, false otherwise.
 */
__inline__ bool accept(const Node &a, const Node &b, const double theta) {
    if (&a == &b)
        return false;
    const double dx = a.com_x - b.com_x;
    const double dy = a.com_y - b.com_y;
    const double dist_sq = dx * dx + dy * dy + 1e-18;
    const double size = std::max({a.box.width(), a.box.height(), b.box.width(), b.box.height()});
    if (theta <= 0.0)
        return false;
    return (size * size) < (theta * theta * dist_sq);
}

/**
 * @brief Determines whether two nodes should be pruned based on a cutoff distance.
 * @param a First node.
 * @param b Second node.
 * @param cutoff Cutoff distance for pruning.
 * @return True if the nodes are beyond the cutoff distance, false otherwise.
 */
__inline__ bool prune(const Node &a, const Node &b, const double cutoff) {
    if (cutoff <= 0.0)
        return false;
    return a.box.distance_sq_lower_bound(b.box) > cutoff * cutoff;
}

template <typename LeafFunc, typename AcceptFunc>
void dual_tree_traversal(const Tree &tree, const double theta, const double cutoff,
                         LeafFunc on_leaf_pair, AcceptFunc on_accepted_pair) {
    // TODO: implement traversal
    (void)tree;
    (void)theta;
    (void)cutoff;
    (void)on_leaf_pair;
    (void)on_accepted_pair;
}

template <typename LeafFunc, typename AcceptFunc>
void dual_tree_traversal_omp(const Tree &tree, const double theta, const double cutoff,
                             LeafFunc on_leaf_pair, AcceptFunc on_accepted_pair) {
    // TODO: implement traversal
    (void)tree;
    (void)theta;
    (void)cutoff;
    (void)on_leaf_pair;
    (void)on_accepted_pair;
}

template <typename Func, typename AcceptFunc>
void dual_tree_traversal_mpi(const Tree &tree, const double theta, const double cutoff,
                             Func on_leaf_pair, AcceptFunc on_accepted_pair, int mpi_rank,
                             int mpi_size) {
    // TODO: implement traversal
    (void)tree;
    (void)theta;
    (void)cutoff;
    (void)on_leaf_pair;
    (void)on_accepted_pair;
    (void)mpi_rank;
    (void)mpi_size;
}

#ifdef DTT_ENABLE_CUDA
template <typename LeafFunc, typename AcceptFunc>
void cu::dual_tree_traversal_cuda(const Tree tree, const double theta, const double cutoff,
                                  LeafFunc on_leaf_pair, AcceptFunc on_accepted_pair);
#endif
} // namespace dtt::tree