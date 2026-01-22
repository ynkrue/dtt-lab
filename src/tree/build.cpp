#include "dtt/tree/build.h"
#include "dtt/tree/bounding-box.h"
#include "dtt/tree/tree.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <vector>

namespace dtt::tree {

dtt::tree::Tree build_quadtree(const dtt::core::ConstParticlesView &particles,
                               const BuildParams &params) {
    Tree tree;
    if (particles.count == 0)
        return tree;

    tree.root = 0;

    // allocate nodes with a generous upper bound to avoid reallocs
    std::size_t capacity = 4 * (particles.count / params.max_leaf_size + 1);
    tree.nodes = dtt::core::Memory<Node>::allocate(capacity);

    // compute root bounding box
    auto [min_x, max_x] = std::minmax_element(particles.x, particles.x + particles.count);
    auto [min_y, max_y] = std::minmax_element(particles.y, particles.y + particles.count);
    BoundingBox box{.min_x = *min_x, .min_y = *min_y, .max_x = *max_x, .max_y = *max_y};
    box.make_square();
    box.pad(params.padding);
    tree.root_box = box;

    // init particles indices
    tree.index_count = particles.count;
    tree.indices = dtt::core::Memory<std::size_t>::allocate(particles.count);

    #pragma omp parallel for if (particles.count >= 1000)
    for (std::size_t i = 0; i < particles.count; ++i)
        tree.indices[i] = i;

    // dtt::core::Memory<uint64_t> morton_codes =
    //     dtt::core::Memory<uint64_t>::allocate(particles.count);

    // util::morton_encode(morton_codes.data(), tree.indices.data(), particles.x, particles.y,
    //                     particles.count, *min_x, *min_y, *max_x, *max_y);
    // util::radix_sort(morton_codes.data(), tree.indices.data(), particles.count);

    // create root node
    Node root{.box = box, .depth = 0, .first = 0, .count = particles.count};
    tree.nodes[0] = root;
    tree.node_count = 1;

    // workspace for permutations and manual stack
    dtt::core::Memory<std::size_t> scratch =
        dtt::core::Memory<std::size_t>::allocate(particles.count);
    std::size_t stack_capacity = capacity;
    auto stack = dtt::core::Memory<std::size_t>::allocate(stack_capacity);
    std::size_t sp = 0;
    stack[sp++] = tree.root;

    while (sp > 0) {
        std::size_t node_idx = stack[--sp];
        while (tree.node_count + 4 > capacity) {
            tree.nodes.extend(capacity << 1);
            capacity <<= 1;
            // keep stack capacity in sync if we need more space
            stack.extend(stack_capacity << 1);
            stack_capacity <<= 1;
        }
        Node &node = tree.nodes[node_idx];
        if (node.count <= params.max_leaf_size || node.depth >= params.max_depth) {
            double mass = 0.0, com_x = 0.0, com_y = 0.0;
            #pragma omp parallel for reduction(+:mass,com_x,com_y) if (node.count >= 2048)
            for (std::size_t i = 0; i < node.count; ++i) {
                std::size_t p_idx = tree.indices[node.first + i];
                const double m = particles.mass[p_idx];
                mass += m;
                com_x += m * particles.x[p_idx];
                com_y += m * particles.y[p_idx];
            }
            node.mass = mass;
            if (mass > 0.0) {
                node.com_x = com_x / mass;
                node.com_y = com_y / mass;
            }
            continue; // leaf node
        }

        const double mid_x = 0.5 * (node.box.min_x + node.box.max_x);
        const double mid_y = 0.5 * (node.box.min_y + node.box.max_y);

        std::size_t counts[4] = {0, 0, 0, 0};
        #pragma omp parallel for reduction(+:counts[:4]) if (node.count >= 2048)
        for (std::size_t i = 0; i < node.count; ++i) {
            std::size_t p_idx = tree.indices[node.first + i];
            int quad = (particles.x[p_idx] >= mid_x) | ((particles.y[p_idx] >= mid_y) << 1);
            counts[quad]++;
        }

        std::size_t offs[4];
        offs[0] = 0;
        for (int q = 1; q < 4; ++q)
            offs[q] = offs[q - 1] + counts[q - 1];

        double child_mass[4] = {0, 0, 0, 0};
        double child_cx[4] = {0, 0, 0, 0};
        double child_cy[4] = {0, 0, 0, 0};
        std::size_t cursor[4] = {offs[0], offs[1], offs[2], offs[3]};
        for (std::size_t i = 0; i < node.count; ++i) {
            std::size_t p_idx = tree.indices[node.first + i];
            int quad = (particles.x[p_idx] >= mid_x) | ((particles.y[p_idx] >= mid_y) << 1);
            scratch[node.first + cursor[quad]++] = p_idx;
            double m = particles.mass[p_idx];
            child_mass[quad] += m;
            child_cx[quad] += m * particles.x[p_idx];
            child_cy[quad] += m * particles.y[p_idx];
        }

        std::memcpy(tree.indices.data() + node.first, scratch.data() + node.first,
                    node.count * sizeof(std::size_t));

        std::array<BoundingBox, 4> child_boxes;
        child_boxes[0] = BoundingBox{node.box.min_x, node.box.min_y, mid_x, mid_y}; // SW
        child_boxes[1] = BoundingBox{mid_x, node.box.min_y, node.box.max_x, mid_y}; // SE
        child_boxes[2] = BoundingBox{node.box.min_x, mid_y, mid_x, node.box.max_y}; // NW
        child_boxes[3] = BoundingBox{mid_x, mid_y, node.box.max_x, node.box.max_y}; // NE

        std::size_t start = node.first;
        double parent_mass = 0.0, parent_com_x = 0.0, parent_com_y = 0.0;
        #pragma unroll
        for (int q = 0; q < 4; ++q) {
            if (counts[q] == 0) {
                node.children[q] = -1;
                continue;
            }
            std::size_t child_tree_idx = tree.node_count++;
            Node &child = tree.nodes[child_tree_idx];
            child.box = child_boxes[q];
            child.depth = node.depth + 1;
            child.first = start;
            child.count = counts[q];
            child.mass = child_mass[q];
            if (child.mass > 0.0) {
                child.com_x = child_cx[q] / child.mass;
                child.com_y = child_cy[q] / child.mass;
            }
            node.children[q] = child_tree_idx;
            if (sp >= stack_capacity) {
                stack.extend(stack_capacity << 1);
                stack_capacity <<= 1;
            }
            stack[sp++] = child_tree_idx;
            parent_mass += child.mass;
            parent_com_x += child.mass * child.com_x;
            parent_com_y += child.mass * child.com_y;
            start += counts[q];
        }
        node.mass = parent_mass;
        if (parent_mass > 0.0) {
            node.com_x = parent_com_x / parent_mass;
            node.com_y = parent_com_y / parent_mass;
        }
    }
    return tree;
}

dtt::tree::Tree build_quadtree_omp(const dtt::core::ConstParticlesView &particles,
                                   const BuildParams &params) {
    // Placeholder for OpenMP parallelized quadtree builder implementation
    // Currently, it calls the serial version
    return build_quadtree(particles, params);
}

} // namespace dtt::tree
