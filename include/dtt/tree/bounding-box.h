#pragma once

namespace dtt::tree
{
    /**
    * @file tree.h
    * @brief Definitions related to tree data structures.
    * Defines structures and functions for managing spatial partitioning trees.
    * @author Yannik Rüfenacht
    */

    /**
     * @brief Axis-aligned bounding box in 2D space.
     * Used for spatial partitioning and interaction modeling.
     * Base structure for quadtree nodes.
     */
    struct BoundingBox {
        // corners
        double min_x{};
        double min_y{};
        double max_x{};
        double max_y{};

        // utility functions
        void pad(double eps);
        double width() const;
        double height() const;
        bool valid() const;
        void clamp(double &x, double &y) const;
        void make_square();
        double distance_sq_lower_bound(const BoundingBox &other) const;
    };
} // namespace dtt::tree