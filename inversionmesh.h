#ifndef INVERSIONMESH_H
#define INVERSIONMESH_H

#include <vector>
#include <array>

// ============================================================
// InversionMesh
// ------------------------------------------------------------
// Solver-facing representation of an arbitrary polyhedral mesh.
// Completely independent of:
//  - grid structure
//  - prism / tetra / hex choice
//  - visualization
// ============================================================

class InversionMesh {
public:
    struct Element {
        int id = -1;

        // Geometry
        double volume = 0.0;
        std::array<double, 3> centroid{0.0, 0.0, 0.0};

        // Topology
        std::vector<int> neighbors;
        std::vector<double> sharedFaceArea;

        // Constraints
        bool active = true;
    };

    InversionMesh() = default;

    // Construction
    void clear();
    void addElement(const Element& e);

    // Access
    size_t size() const;
    const Element& element(size_t i) const;
    Element& element(size_t i);

    const std::vector<Element>& elements() const;

private:
    std::vector<Element> elements_;
};

#endif // INVERSIONMESH_H
