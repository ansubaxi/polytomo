#include "meshadapter.h"

// ============================================================
// Convert prisms → inversion elements
// ============================================================

InversionMesh MeshAdapter::buildInversionMesh(
        const ModelGenerator::MeshData& mesh) {

    InversionMesh invMesh;
    invMesh.clear();

    invMesh.clear();
    invMesh = InversionMesh();

    for (const auto& prism : mesh.prisms) {
        InversionMesh::Element e;

        e.id = prism.id;
        e.volume = prism.volume;
        e.centroid = prism.centroid;
        e.active = prism.active;

        // Topology
        e.neighbors = prism.neighbors;
        e.sharedFaceArea = prism.sharedFaceArea;

        invMesh.addElement(e);
    }

    return invMesh;
}
