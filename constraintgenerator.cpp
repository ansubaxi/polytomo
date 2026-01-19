#include "constraintgenerator.h"

// ============================================================
// Public API
// ============================================================

void ConstraintGenerator::setPolygonConstraint(const Polygon2D& poly) {
    polygon_ = poly;
    hasPolygon_ = true;
}

void ConstraintGenerator::setZRangeConstraint(const ZRange& zrange) {
    zrange_ = zrange;
    hasZRange_ = true;
}

void ConstraintGenerator::apply(std::vector<Cell>& cells) const {
    for (auto& c : cells) {
        if (!c.active) continue;

        // ----- Polygon constraint -----
        if (hasPolygon_) {
            bool inside = pointInPolygon(
                c.centroid[0], c.centroid[1], polygon_.vertices);

            if (!polygon_.invert && !inside)
                c.active = false;

            if (polygon_.invert && inside)
                c.active = false;
        }

        // ----- Z-range constraint -----
        if (hasZRange_) {
            double depth = -c.centroid[2]; // positive down
            if (depth < zrange_.zMin || depth > zrange_.zMax)
                c.active = false;
        }
    }
}

// ============================================================
// Geometry utilities
// ============================================================

// Ray casting algorithm (robust, fast)
bool ConstraintGenerator::pointInPolygon(
    double x, double y,
    const std::vector<std::array<double,2>>& poly) const
{
    bool inside = false;
    size_t n = poly.size();

    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        double xi = poly[i][0], yi = poly[i][1];
        double xj = poly[j][0], yj = poly[j][1];

        bool intersect =
            ((yi > y) != (yj > y)) &&
            (x < (xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi);

        if (intersect)
            inside = !inside;
    }

    return inside;
}
