#ifndef CONSTRAINTGENERATOR_H
#define CONSTRAINTGENERATOR_H

#include <vector>
#include <array>
#include <string>

// ============================================================
// ConstraintGenerator
// ------------------------------------------------------------
// Applies geometric constraints to an existing mesh by
// marking cells as active/inactive.
// Operates ONLY on cell centroids.
// ============================================================

class ConstraintGenerator {
public:
    struct Polygon2D {
        std::vector<std::array<double, 2>> vertices; // (x,y)
        bool invert = false; // true = mask outside
    };

    struct ZRange {
        double zMin;  // depth-positive down
        double zMax;
    };

    struct Cell {
        int id;
        std::array<double, 3> centroid; // (x,y,z)
        bool active;
    };

    // Construction
    ConstraintGenerator() = default;

    // Load constraints
    void setPolygonConstraint(const Polygon2D& poly);
    void setZRangeConstraint(const ZRange& zrange);

    // Apply constraints
    void apply(std::vector<Cell>& cells) const;

private:
    bool pointInPolygon(double x, double y,
                        const std::vector<std::array<double,2>>& poly) const;

private:
    bool hasPolygon_ = false;
    bool hasZRange_  = false;

    Polygon2D polygon_;
    ZRange zrange_;
};

#endif // CONSTRAINTGENERATOR_H
