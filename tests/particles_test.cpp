#include "dtt/sim/particles.h"

#include <gtest/gtest.h>

namespace dtt::sim {

TEST(ParticlesTest, ValidateSizes) {
    Particles p;
    p.x = {0.0, 1.0};
    p.y = {0.0, 1.0};
    p.vx = {0.0, 0.0};
    p.vy = {0.0, 0.0};
    p.mass = {1.0, 2.0};
    EXPECT_TRUE(validate_particles(p));

    p.mass.pop_back();
    EXPECT_FALSE(validate_particles(p));
}

TEST(ParticlesTest, RandomGeneratorProducesSizes) {
    const std::size_t n = 5;
    const double min_x = -1.0, max_x = 2.0;
    const double min_y = -3.0, max_y = 4.0;
    const double mass_min = 0.5, mass_max = 2.0;

    Particles p =
        create_rnd_particles(n, /*seed=*/42, min_x, max_x, min_y, max_y, mass_min, mass_max);
    EXPECT_EQ(p.x.size(), n);
    EXPECT_EQ(p.y.size(), n);
    EXPECT_EQ(p.vx.size(), n);
    EXPECT_EQ(p.vy.size(), n);
    EXPECT_EQ(p.mass.size(), n);

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_GE(p.x[i], min_x);
        EXPECT_LE(p.x[i], max_x);
        EXPECT_GE(p.y[i], min_y);
        EXPECT_LE(p.y[i], max_y);
        EXPECT_GE(p.mass[i], mass_min);
        EXPECT_LE(p.mass[i], mass_max);
        EXPECT_DOUBLE_EQ(p.vx[i], 0.0);
        EXPECT_DOUBLE_EQ(p.vy[i], 0.0);
    }
}

} // namespace dtt::sim
