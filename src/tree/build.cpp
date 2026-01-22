#include "dtt/tree/build.h"
#include "dtt/tree/bounding-box.h"
#include "dtt/tree/tree.h"

#include <algorithm>
#include <array>
#include <vector>

namespace dtt::tree {

dtt::tree::Tree build_quadtree(const dtt::core::ConstParticlesView &particles,
                               const BuildParams &params) {
    Tree tree;
    if (particles.count == 0)
        return tree;

    tree.root = 0;

    // allocate nodes
    std::size_t capacity = 2 * (particles.count / params.max_leaf_size + 1);
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

    // build tree using stack
    std::vector<std::size_t> stack;
    stack.reserve(capacity);
    stack.push_back(tree.root);

    while (!stack.empty()) {
        std::size_t node_idx = stack.back();
        stack.pop_back();
        // Make sure we never reallocate while holding references to nodes.
        while (tree.node_count + 4 > capacity) {
            tree.nodes.extend(capacity << 1);
            capacity <<= 1;
        }
        Node &node = tree.nodes[node_idx];
        if (node.count <= params.max_leaf_size || node.depth >= params.max_depth) {
            double mass = 0.0, com_x = 0.0, com_y = 0.0;
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

        // partition particles into 4 quadrants
        double mid_x = 0.5 * (node.box.min_x + node.box.max_x);
        double mid_y = 0.5 * (node.box.min_y + node.box.max_y);
        // define child boxes
        std::array<BoundingBox, 4> child_boxes;
        child_boxes[0] = BoundingBox{node.box.min_x, node.box.min_y, mid_x, mid_y}; // SW
        child_boxes[1] = BoundingBox{mid_x, node.box.min_y, node.box.max_x, mid_y}; // SE
        child_boxes[2] = BoundingBox{node.box.min_x, mid_y, mid_x, node.box.max_y}; // NW
        child_boxes[3] = BoundingBox{mid_x, mid_y, node.box.max_x, node.box.max_y}; // NE
        // center of mass positions
        std::array<double, 12> child_coms{0}; // com_x, com_y, mass

        std::vector<std::size_t> quadrants[4];
        #pragma omp parallel for if (node.count >= 1000)
        for (std::size_t i = 0; i < node.count; ++i) {
            std::size_t p_idx = tree.indices[node.first + i];
            int quad = (particles.x[p_idx] >= mid_x) | ((particles.y[p_idx] >= mid_y) << 1);
            quadrants[quad].push_back(p_idx);
            child_coms[3 * quad + 0] += particles.mass[p_idx] * particles.x[p_idx];
            child_coms[3 * quad + 1] += particles.mass[p_idx] * particles.y[p_idx];
            child_coms[3 * quad + 2] += particles.mass[p_idx];
        }

        // create child nodes
        std::size_t child_offset = 0;
        double parent_mass = 0.0, parent_com_x = 0.0, parent_com_y = 0.0;

        #pragma unroll
        for (int q = 0; q < 4; ++q) {
            // check for empty quadrant and max capacity
            if (quadrants[q].empty()) {
                node.children[q] = -1;
                continue;
            }

            // write back permuted indices
            std::size_t child_first = node.first + child_offset;
            std::size_t child_count = quadrants[q].size();
            for (std::size_t p_idx : quadrants[q]) {
                tree.indices[node.first + child_offset++] = p_idx;
            }

            // create child node
            std::size_t child_tree_idx = tree.node_count++;
            Node &child = tree.nodes[child_tree_idx];
            child.box = child_boxes[q];
            child.depth = node.depth + 1;
            child.first = child_first;
            child.count = child_count;
            child.mass = child_coms[3 * q + 2];
            if (child.mass > 0.0) {
                child.com_x = child_coms[3 * q] / child.mass;
                child.com_y = child_coms[3 * q + 1] / child.mass;
            }

            node.children[q] = child_tree_idx;
            stack.push_back(child_tree_idx);
            parent_mass += child.mass;
            parent_com_x += child.mass * child.com_x;
            parent_com_y += child.mass * child.com_y;
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
