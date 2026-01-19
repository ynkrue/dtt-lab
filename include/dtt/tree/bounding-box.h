#pragma once

namespace dtt::tree {
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

    /// put padding around the box
    void pad(double eps);
    /// get shape of the box
    double width() const;
    double height() const;
    /// check for valid corners
    bool valid() const;
    /// force point to be inside the box
    void clamp(double &x, double &y) const;
    /// expand box to include point
    void expand(double x, double y);
    /// make square by expanding the smaller side
    void make_square();
    /// compute squared distance between two boxes
    double distance_sq_lower_bound(const BoundingBox &other) const;
};
} // namespace dtt::tree