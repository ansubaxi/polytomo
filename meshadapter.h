#ifndef MESHADAPTER_H
#define MESHADAPTER_H

#include "modelgenerator.h"
#include "inversionmesh.h"

// ============================================================
// MeshAdapter
// ------------------------------------------------------------
// Converts ModelGenerator mesh (prisms + vertices)
// into solver-facing InversionMesh elements.
// ============================================================

class MeshAdapter {
public:
    static InversionMesh buildInversionMesh(
        const ModelGenerator::MeshData& mesh);
};

#endif // MESHADAPTER_H
