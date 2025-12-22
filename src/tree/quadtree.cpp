#include "dtt/tree/quadtree.h"

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>

namespace dtt::tree {

BoundingBox compute_root_bbox(const sim::Particles &particles, double padding) {
    BoundingBox bbox;
    bbox.min_x = std::numeric_limits<double>::infinity();
    bbox.min_y = std::numeric_limits<double>::infinity();
    bbox.max_x = -std::numeric_limits<double>::infinity();
    bbox.max_y = -std::numeric_limits<double>::infinity();

    for (std::size_t i = 0; i < particles.x.size(); ++i) {
        bbox.expand(particles.x[i], particles.y[i]);
    }
    bbox.pad(padding);
    bbox.make_square();
    return bbox;
}

Tree build_quadtree(const sim::Particles &particles, const BuildParams &params) {
    const std::size_t n = particles.x.size();
    Tree tree;
    if (n == 0) return tree;

    // compute root bounding box
    BoundingBox root_bbox = compute_root_bbox(particles, params.padding);
    tree.root_bbox = root_bbox;

    // compute morton keys
    std::vector<std::pair<uint64_t, std::size_t>> morton_keys;
    morton_keys.reserve(n);
    const uint32_t levels = static_cast<uint32_t>(std::min(params.max_depth, 21));
    
    for (std::size_t i = 0; i < n; ++i) {
        double nx, ny;
        normalize_point(root_bbox, particles.x[i], particles.y[i], nx, ny);
        morton_keys.push_back({morton_from_normalized(nx, ny, levels), i});
    }

    // sort by morton code
    std::sort(morton_keys.begin(), morton_keys.end(), [](const auto &a, const auto &b) {
        return a.first < b.first;
    });

    // fill in permutation
    tree.indices.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        tree.indices[i] = morton_keys[i].second;
    }

    
    // build tree nodes
    tree.nodes.clear();
    tree.nodes.reserve(2 * n); // rough estimate
    
    Node root{};
    root.bbox = root_bbox;
    std::fill(std::begin(root.children), std::end(root.children), -1);
    root.first = 0;
    root.count = static_cast<int>(n);
    root.mass = 0.0;
    root.com_x = root.com_y = 0.0;
    root.level = 0;

    tree.nodes.push_back(root);
    tree.root = 0;

    std::vector<int> stack;
    stack.push_back(tree.root);

    while (!stack.empty()) {
        int ni = stack.back();
        stack.pop_back();
        Node &node = tree.nodes[ni];

        if (node.count <= params.max_leaf_size || node.level >= params.max_depth) {
            // leaf node, compute mass and center of mass
            BoundingBox bbox;
            bbox.min_x = bbox.min_y = std::numeric_limits<double>::infinity();
            bbox.max_x = bbox.max_y = -std::numeric_limits<double>::infinity();
            double total_mass = 0.0, com_x = 0.0, com_y = 0.0;
            for (int k = 0; k < node.count; ++k) {
                std::size_t idx = tree.indices[node.first + k];
                bbox.expand(particles.x[idx], particles.y[idx]);
                total_mass += particles.mass[idx];
                com_x += particles.mass[idx] * particles.x[idx];
                com_y += particles.mass[idx] * particles.y[idx];
            }
            node.bbox = bbox;
            node.mass = total_mass;
            node.com_x = total_mass > 0.0 ? com_x / total_mass : 0.0;
            node.com_y = total_mass > 0.0 ? com_y / total_mass : 0.0;
            continue;
        }

        // split into quadrants by midpoint
        const double mx = 0.5 * (node.bbox.min_x + node.bbox.max_x);
        const double my = 0.5 * (node.bbox.min_y + node.bbox.max_y);

        std::array<std::vector<std::size_t>, 4> buckets;
        for (int k = 0; k < node.count; ++k) {
            std::size_t idx = tree.indices[node.first + k];
            // put into buckets: 0:SW,1:SE,2:NW,3:NE
            const int quad = (particles.x[idx] >= mx) | ((particles.y[idx] >= my) << 1);
            buckets[quad].push_back(idx);
        }

        // write back permuted indices
        int write = node.first;
        int child_first[4];
        int child_count[4];
        for (int q = 0; q < 4; ++q) {
            child_first[q] = write;
            child_count[q] = static_cast<int>(buckets[q].size());
            for (std::size_t idx : buckets[q]) tree.indices[write++] = idx;
        }

        // create child nodes
        for (int q = 0; q < 4; ++q) {
            if (child_count[q] == 0) {
                node.children[q] = -1;
                continue;
            }
            BoundingBox child_bbox;
            child_bbox.min_x = child_bbox.min_y = std::numeric_limits<double>::infinity();
            child_bbox.max_x = child_bbox.max_y = -std::numeric_limits<double>::infinity();
            double mass = 0.0, com_x = 0.0, com_y = 0.0;
            for (int k = 0; k < child_count[q]; ++k) {
                std::size_t idx = tree.indices[child_first[q] + k];
                child_bbox.expand(particles.x[idx], particles.y[idx]);
                mass += particles.mass[idx];
                com_x += particles.mass[idx] * particles.x[idx];
                com_y += particles.mass[idx] * particles.y[idx];
            }
            Node child{};
            child.bbox = child_bbox;
            std::fill(std::begin(child.children), std::end(child.children), -1);
            child.first = child_first[q];
            child.count = child_count[q];
            child.mass = mass;
            child.com_x = mass > 0.0 ? com_x / mass : 0.0;
            child.com_y = mass > 0.0 ? com_y / mass : 0.0;
            child.level = node.level + 1;
            tree.nodes.push_back(child);
            node.children[q] = static_cast<int>(tree.nodes.size() - 1);
            stack.push_back(node.children[q]);
        }
        // update parent mass/COM/bbox from children
        double total_mass = 0.0, com_x = 0.0, com_y = 0.0;
        BoundingBox tight;
        tight.min_x = tight.min_y = std::numeric_limits<double>::infinity();
        tight.max_x = tight.max_y = -std::numeric_limits<double>::infinity();
        for (int q = 0; q < 4; ++q) {
            int ci = node.children[q];
            if (ci == -1) continue;
            const Node &ch = tree.nodes[ci];
            total_mass += ch.mass;
            com_x += ch.mass * ch.com_x;
            com_y += ch.mass * ch.com_y;
            tight.expand(ch.bbox.min_x, ch.bbox.min_y);
            tight.expand(ch.bbox.max_x, ch.bbox.max_y);
        }
        node.mass = total_mass;
        node.com_x = total_mass > 0.0 ? com_x / total_mass : 0.0;
        node.com_y = total_mass > 0.0 ? com_y / total_mass : 0.0;
        node.bbox = tight;
    }
    return tree;
}

bool accept(const Node &a, const Node &b, const MAC &mac) {
    if (&a == &b) return false;
    const double dx = a.com_x - b.com_x;
    const double dy = a.com_y - b.com_y;
    const double dist_sq = dx * dx + dy * dy + 1e-18;
    const double size = std::max(a.bbox.max_extent(), b.bbox.max_extent());
    if (mac.theta <= 0.0) return false;
    return (size * size) < (mac.theta * mac.theta * dist_sq);
}

}
