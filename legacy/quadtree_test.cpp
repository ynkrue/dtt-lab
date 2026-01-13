#include "dtt/sim/forces.h"
#include "dtt/tree/quadtree.h"

#include <gtest/gtest.h>

namespace {

using dtt::sim::compute_forces_naive;
using dtt::sim::compute_forces_tree;
using dtt::sim::ForceField;
using dtt::sim::ForceParams;
using dtt::sim::Particles;
using dtt::tree::BuildParams;
using dtt::tree::Tree;

Particles make_test_particles() {
    Particles p;
    p.x = {0.0, 1.0, 0.0, 1.0};
    p.y = {0.0, 0.0, 1.0, 1.0};
    p.vx = {0.0, 0.0, 0.0, 0.0};
    p.vy = {0.0, 0.0, 0.0, 0.0};
    p.mass = {1.0, 1.0, 1.0, 1.0};
    return p;
}

TEST(QuadtreeTest, BuildProducesValidLeaves) {
    Particles p = make_test_particles();
    BuildParams params;
    params.max_leaf_size = 2;
    Tree t = dtt::tree::build_quadtree(p, params);

    ASSERT_NE(t.root, -1);
    ASSERT_FALSE(t.nodes.empty());
    EXPECT_EQ(t.nodes[t.root].count, static_cast<int>(p.x.size()));

    // Check every leaf has size <= max_leaf_size.
    for (const auto &n : t.nodes) {
        if (n.is_leaf()) {
            EXPECT_LE(n.count, params.max_leaf_size);
        }
    }
    // Mass is propagated to root.
    EXPECT_NEAR(t.nodes[t.root].mass, 4.0, 1e-12);
}

TEST(QuadtreeTest, ForcesMatchNaiveWhenNoAccept) {
    Particles p = make_test_particles();
    BuildParams params;
    params.max_leaf_size = 1; // force full refinement
    Tree t = dtt::tree::build_quadtree(p, params);

    ForceParams fp{.softening = 1e-6, .cutoff = std::nullopt};
    dtt::tree::MAC mac{.theta = 0.0}; // never accept, only leaf-leaf

    ForceField f_tree = compute_forces_tree(p, t, fp, mac);
    ForceField f_naive = compute_forces_naive(p, fp);

    ASSERT_EQ(f_tree.size(), f_naive.size());
    for (std::size_t i = 0; i < f_tree.size(); ++i) {
        EXPECT_NEAR(f_tree[i][0], f_naive[i][0], 1e-9);
        EXPECT_NEAR(f_tree[i][1], f_naive[i][1], 1e-9);
    }
}

} // namespace
