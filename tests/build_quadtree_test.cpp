#include "dtt/core/particles.h"
#include "dtt/tree/build.h"

#include <gtest/gtest.h>

#include <functional>
#include <vector>

using dtt::core::ConstParticlesView;
using dtt::core::ParticlesBuffer;
using dtt::tree::BuildParams;
using dtt::tree::Tree;
using BuildFn = Tree (*)(const ConstParticlesView &, const BuildParams &);

struct BuilderCase {
    const char *name;
    BuildFn fn;
};

class BuildTreeTest : public ::testing::TestWithParam<BuilderCase> {};

namespace {

ParticlesBuffer make_particles(const std::vector<double> &xs, const std::vector<double> &ys,
                               const std::vector<double> &masses) {
    const std::size_t n = xs.size();
    ParticlesBuffer buf = ParticlesBuffer::make(n);
    for (std::size_t i = 0; i < n; ++i) {
        buf.x[i] = xs[i];
        buf.y[i] = ys[i];
        buf.mass[i] = masses[i];
        buf.vx[i] = 0.0;
        buf.vy[i] = 0.0;
    }
    return buf;
}

TEST_P(BuildTreeTest, HandlesEmptyInput) {
    ParticlesBuffer buf = ParticlesBuffer::make(0);
    Tree t = GetParam().fn(buf.const_view(), BuildParams{});
    EXPECT_EQ(t.index_count, 0u);
    EXPECT_EQ(t.node_count, 0u);
    EXPECT_EQ(t.root, -1);
}

TEST_P(BuildTreeTest, ComputesRootBoundingBox) {
    ParticlesBuffer buf =
        make_particles({0.0, 1.0, 0.0, 1.0}, {0.0, 0.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 1.0});
    BuildParams params;
    params.padding = 0.0;
    Tree t = GetParam().fn(buf.const_view(), params);
    ASSERT_NE(t.root, -1);
    EXPECT_DOUBLE_EQ(t.root_box.min_x, 0.0);
    EXPECT_DOUBLE_EQ(t.root_box.min_y, 0.0);
    EXPECT_DOUBLE_EQ(t.root_box.max_x, 1.0);
    EXPECT_DOUBLE_EQ(t.root_box.max_y, 1.0);
}

TEST_P(BuildTreeTest, SplitsIntoLeavesRespectingMaxLeafSize) {
    ParticlesBuffer buf =
        make_particles({0.1, 0.9, 0.1, 0.9}, {0.1, 0.1, 0.9, 0.9}, {1.0, 1.0, 1.0, 1.0});
    BuildParams params;
    params.max_leaf_size = 1;
    params.max_depth = 8;
    params.padding = 0.0;

    Tree t = GetParam().fn(buf.const_view(), params);
    ASSERT_NE(t.root, -1);
    const auto &root = t.nodes[t.root];
    // All 4 children should exist and each contain exactly one particle.
    for (int c = 0; c < 4; ++c) {
        ASSERT_NE(root.children[c], -1) << "child " << c << " missing";
        const auto &ch = t.nodes[root.children[c]];
        EXPECT_EQ(ch.count, 1);
        EXPECT_GE(ch.depth, root.depth + 1);
        // Force the child range to be within indices buffer.
        ASSERT_LT(static_cast<std::size_t>(ch.first + ch.count), t.index_count + 1);
    }
}

TEST_P(BuildTreeTest, RespectsMaxDepth) {
    ParticlesBuffer buf =
        make_particles({0.1, 0.9, 0.1, 0.9}, {0.1, 0.1, 0.9, 0.9}, {1.0, 1.0, 1.0, 1.0});
    BuildParams params;
    params.max_leaf_size = 1;
    params.max_depth = 0; // disallow any split
    params.padding = 0.0;

    Tree t = GetParam().fn(buf.const_view(), params);
    ASSERT_NE(t.root, -1);
    const auto &root = t.nodes[t.root];
    EXPECT_EQ(root.children[0], -1);
    EXPECT_EQ(root.children[1], -1);
    EXPECT_EQ(root.children[2], -1);
    EXPECT_EQ(root.children[3], -1);
    EXPECT_EQ(root.count, static_cast<int>(buf.count));
}

TEST_P(BuildTreeTest, ComputesMassAndCenterOfMass) {
    ParticlesBuffer buf = make_particles({0.0, 1.0}, {0.0, 0.0}, {2.0, 4.0});
    BuildParams params;
    params.max_leaf_size = 1;
    params.padding = 0.0;

    Tree t = GetParam().fn(buf.const_view(), params);
    ASSERT_NE(t.root, -1);
    const auto &root = t.nodes[t.root];
    // Root mass and COM should reflect both particles.
    EXPECT_DOUBLE_EQ(root.mass, 6.0);
    EXPECT_DOUBLE_EQ(root.com_x, (2.0 * 0.0 + 4.0 * 1.0) / 6.0);
    EXPECT_DOUBLE_EQ(root.com_y, 0.0);
}

std::string BuilderName(const ::testing::TestParamInfo<BuilderCase> &info) {
    return info.param.name;
}

INSTANTIATE_TEST_SUITE_P(Builders, BuildTreeTest,
                         ::testing::Values(BuilderCase{"morton", &dtt::tree::build_quadtree}
                                           // Add more builders (e.g., {"omp",
                                           // &dtt::tree::build_quadtree_omp}) once implemented
                                           ),
                         BuilderName);

} // namespace
