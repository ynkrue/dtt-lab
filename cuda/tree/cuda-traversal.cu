#include "dtt/tree/cu/cuda-traversal.h"
#include "dtt/tree/cu/tree.h"

namespace dtt::tree::cu {

template <typename LeafFunc, typename AcceptFunc>
__global__ void
dual_tree_traversal_cuda_kernel(Node *nodes, std::size_t *node_count, std::size_t *indices,
                                std::size_t *index_count, const double *theta, const double *cutoff,
                                LeafFunc on_leaf_pair, AcceptFunc on_accepted_pair) {
    // CUDA kernel implementation would go here
    (void)nodes;
    (void)node_count;
    (void)indices;
    (void)index_count;
    (void)theta;
    (void)cutoff;
    (void)on_leaf_pair;
    (void)on_accepted_pair;
}

template <typename LeafFunc, typename AcceptFunc>
void dual_tree_traversal_cuda(const Tree *tree, const double theta, const double cutoff,
                              LeafFunc on_leaf_pair, AcceptFunc on_accepted_pair) {
    (void)tree;
    (void)theta;
    (void)cutoff;
    (void)on_leaf_pair;
    (void)on_accepted_pair;
}

} // namespace dtt::tree::cu