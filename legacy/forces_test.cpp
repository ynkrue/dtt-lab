#include "dtt/sim/forces.h"
#include "dtt/sim/particles.h"

#include <gtest/gtest.h>

namespace dtt::sim {

TEST(ForcesTest, TwoParticleSymmetry) {
    Particles p;
    p.x = {0.0, 1.0};
    p.y = {0.0, 0.0};
    p.vx = {0.0, 0.0};
    p.vy = {0.0, 0.0};
    p.mass = {1.0, 1.0};

    ForceParams params{.softening = 1e-9, .cutoff = std::nullopt, .gravity = 1.0};
    ForceField f = compute_forces_naive(p, params);

    ASSERT_EQ(f.size(), 2u);
    EXPECT_NEAR(f[0][0], 1.0, 1e-9);
    EXPECT_NEAR(f[0][1], 0.0, 1e-12);
    EXPECT_NEAR(f[1][0], -1.0, 1e-9);
    EXPECT_NEAR(f[1][1], 0.0, 1e-12);
}

TEST(ForcesTest, RespectsCutoff) {
    Particles p;
    p.x = {0.0, 10.0};
    p.y = {0.0, 0.0};
    p.vx = {0.0, 0.0};
    p.vy = {0.0, 0.0};
    p.mass = {1.0, 1.0};

    ForceParams params{.softening = 1e-6, .cutoff = 1.0, .gravity = 1.0};
    ForceField f = compute_forces_naive(p, params);

    EXPECT_NEAR(f[0][0], 0.0, 1e-12);
    EXPECT_NEAR(f[1][0], 0.0, 1e-12);
}

TEST(ForcesTest, EulerStepUpdatesPositionAndVelocity) {
    Particles p;
    p.x = {0.0};
    p.y = {0.0};
    p.vx = {0.0};
    p.vy = {0.0};
    p.mass = {1.0};

    ForceField f(1);
    f[0] = {2.0, 0.0}; // ax = 2.0

    const double dt = 0.5;
    euler_step(p, f, dt);

    EXPECT_NEAR(p.vx[0], 1.0, 1e-12); // vx = 0 + 2*0.5
    EXPECT_NEAR(p.x[0], 0.5, 1e-12);  // x = 0 + vx_new * dt
    EXPECT_NEAR(p.vy[0], 0.0, 1e-12);
    EXPECT_NEAR(p.y[0], 0.0, 1e-12);
}

} // namespace dtt::sim
