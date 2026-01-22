#include "dtt/tree/build.h"
#include "dtt/tree/bounding-box.h"
#include "dtt/tree/tree.h"

#include <algorithm>
#include <array>
#include <vector>

namespace dtt::tree {

dtt::tree::Tree build_quadtree_base(const dtt::core::ConstParticlesView &particles,
                                    const BuildParams &params) {
    Tree tree;
    if (particles.count == 0)
        return tree;

    tree.root = 0;

    // init particles indices
    tree.index_count = particles.count;
    tree.indices = dtt::core::Memory<std::size_t>::allocate(particles.count);
    for (std::size_t i = 0; i < particles.count; ++i)
        tree.indices[i] = i;

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

    // create root node
    Node root{.box = box, .depth = 0, .first = 0, .count = particles.count};
    tree.nodes[0] = root;
    tree.node_count = 1;

    // build tree using stack
    std::vector<std::size_t> stack;
    stack.reserve(128);
    stack.push_back(tree.root);

    while (!stack.empty()) {
        std::size_t node_idx = stack.back();
        stack.pop_back();
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
            if (tree.node_count >= capacity) {
                tree.nodes.extend(capacity << 1);
                capacity <<= 1;
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

dtt::tree::Tree build_quadtree(const dtt::core::ConstParticlesView &particles,
                               const BuildParams &params) {
    Tree tree;
    if (particles.count == 0)
        return tree;

    tree.root = 0;

    // compute root bounding box
    auto [min_x, max_x, min_y, max_y] = util::minmax_box(particles.x, particles.y, particles.count);
    BoundingBox box{.min_x = min_x, .min_y = min_y, .max_x = max_x, .max_y = max_y};
    box.make_square();
    box.pad(params.padding);
    tree.root_box = box;

    // init particles indices using morton codes
    tree.index_count = particles.count;
    tree.indices = dtt::core::Memory<std::size_t>::allocate(particles.count);
    dtt::core::Memory<uint64_t> morton_codes =
        dtt::core::Memory<uint64_t>::allocate(particles.count);

    util::morton_encode(morton_codes.data(), tree.indices.data(), particles.x, particles.y,
                        particles.count, box.min_x, box.min_y, box.max_x, box.max_y);
    util::radix_sort(morton_codes.data(), tree.indices.data(), particles.count);

    // allocate nodes
    std::size_t capacity = 2 * (particles.count / params.max_leaf_size + 1);
    tree.nodes = dtt::core::Memory<Node>::allocate(capacity);

    // create root node
    Node root{.box = box, .depth = 0, .first = 0, .count = particles.count};
    tree.nodes[0] = root;
    tree.node_count = 1;

    // build tree using memory-efficient stack
    std::vector<std::size_t> stack;
    stack.reserve(capacity);
    stack.push_back(tree.root);

    // workspace for child partitions
    dtt::core::Memory<std::size_t> scratch =
        dtt::core::Memory<std::size_t>::allocate(particles.count);

    
    while (!stack.empty()) {
        auto ni = stack.back(); stack.pop_back();
        Node &node = tree.nodes[ni];
        if (node.count <= params.max_leaf_size || node.depth >= params.max_depth) {
            double mass=0, cx=0, cy=0;
            #pragma omp parallel for reduction(+:mass,cx,cy) if (node.count>=500)
            for (std::size_t i=0;i<node.count;++i){
                auto pid = tree.indices[node.first+i];
                double m = particles.mass[pid];
                mass += m; cx += m*particles.x[pid]; cy += m*particles.y[pid];
            }
            node.mass = mass;
            if (mass>0){ node.com_x = cx/mass; node.com_y = cy/mass; }
            continue;
        }
        double mx = 0.5*(node.box.min_x+node.box.max_x);
        double my = 0.5*(node.box.min_y+node.box.max_y);
        std::size_t counts[4] = {0,0,0,0};
        // count
        #pragma omp parallel for reduction(+:counts[:4]) if (node.count>=500)
        for (std::size_t i=0;i<node.count;++i){
            auto pid = tree.indices[node.first+i];
            int q = (particles.x[pid] >= mx) | ((particles.y[pid] >= my) << 1);
            counts[q]++;
        }
        // offsets
        std::size_t offs[4]; offs[0]=0;
        for (int q=1;q<4;++q) offs[q]=offs[q-1]+counts[q-1];
        // scatter + accumulate child mass/com
        double child_mass[4]={0,0,0,0}, child_cx[4]={0,0,0,0}, child_cy[4]={0,0,0,0};
        #pragma omp parallel for reduction(+:child_mass[:4],child_cx[:4],child_cy[:4]) if (node.count>=500)
        for (std::size_t i=0;i<node.count;++i){
            auto pid = tree.indices[node.first+i];
            int q = (particles.x[pid] >= mx) | ((particles.y[pid] >= my) << 1);
            scratch[node.first + offs[q]++] = pid;
            double m = particles.mass[pid];
            child_mass[q] += m;
            child_cx[q] += m*particles.x[pid];
            child_cy[q] += m*particles.y[pid];
        }
        // copy back
        std::memcpy(tree.indices.data()+node.first, scratch.data()+node.first,
                    node.count*sizeof(std::size_t));
        // create children
        std::size_t start = node.first;
        double parent_mass=0, parent_cx=0, parent_cy=0;
        #pragma unroll
        for(int q=0;q<4;++q){
            if (counts[q]==0){ node.children[q]=-1; continue; }
            if (tree.node_count >= capacity){ tree.nodes.extend(capacity<<1); capacity<<=1; }
            std::size_t child_idx = tree.node_count++;
            Node &ch = tree.nodes[child_idx];
            ch.box = {/* optional tight bbox, or keep quadrant box */ node.box.min_x, node.box.min_y,
                    node.box.max_x, node.box.max_y}; // keep simple
            ch.depth = node.depth+1;
            ch.first = start;
            ch.count = counts[q];
            ch.mass = child_mass[q];
            if (child_mass[q]>0){ ch.com_x = child_cx[q]/child_mass[q]; ch.com_y = child_cy[q]/child_mass[q]; }
            node.children[q] = static_cast<int>(child_idx);
            stack.push_back(child_idx);
            parent_mass += ch.mass;
            parent_cx += ch.mass * ch.com_x;
            parent_cy += ch.mass * ch.com_y;
            start += counts[q];
        }
        node.mass = parent_mass;
        if (parent_mass>0){ node.com_x = parent_cx/parent_mass; node.com_y = parent_cy/parent_mass; }
    }

    return tree;
}

dtt::tree::Tree build_quadtree_omp(const dtt::core::ConstParticlesView &particles,
                                   const BuildParams &params) {
    (void)particles;
    (void)params;
    return Tree{};
}

} // namespace dtt::tree
