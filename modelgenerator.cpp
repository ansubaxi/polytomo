#include "modelgenerator.h"

#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QRegExp>
#include <unordered_map>
#include <algorithm>
#include <cmath>

// ============================================================
// Constructor / Access
// ============================================================

ModelGenerator::ModelGenerator(QObject* parent)
    : QObject(parent) {}

ModelGenerator::~ModelGenerator() {}

void ModelGenerator::setParameters(const Parameters& params) {
    params_ = params;
}

const ModelGenerator::Parameters& ModelGenerator::parameters() const {
    return params_;
}

const ModelGenerator::MeshData& ModelGenerator::mesh() const {
    return mesh_;
}

// ============================================================
// Compute extents for UI display (before full generation)
// ============================================================

bool ModelGenerator::computeExtentFromParameters() {
    if (!readInputFile()) {
        return false;
    }
    
    computeExtent();
    buildSurfaceGrid();
    
    return true;
}

// ============================================================
// Main pipeline
// ============================================================

void ModelGenerator::process() {
    emit logMessage("Starting triangular-prism mesh generation");

    if (!readInputFile()) {
        emit finished(false, "Failed to read input file");
        return;
    }

    computeExtent();
    buildSurfaceGrid();
    buildLayeredVertices();
    buildPrisms();
    computePrismGeometry();
    buildAdjacency();

    if (!validateMesh()) {
        emit finished(false, "Mesh validation failed");
        return;
    }

    if (!writeMeshFiles()) {
        emit finished(false, "Failed to write output files");
        return;
    }

    // FIXED: Changed to 2-param signature, embed maxElevation in message
    emit finished(true, QString("Mesh generation completed successfully. Max elevation: %1 m").arg(mesh_.zMax));
}

// ============================================================
// Input parsing
// ============================================================

bool ModelGenerator::readInputFile() {
    QFile file(params_.inputFile);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return false;

    QTextStream in(&file);

    mesh_.x.clear();
    mesh_.y.clear();
    mesh_.elevation.clear();
    
    // Track lines to skip empty space
    int lineIndex = 0;

    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        if (line.isEmpty()) continue;

        QStringList t = line.split(QRegExp("[,\\s]+"), Qt::SkipEmptyParts);
        int maxCol = std::max({params_.eastingColumn,
                               params_.northingColumn,
                               params_.topoColumn});
        
        // Ensure line has enough columns
        if (t.size() <= maxCol) continue;
        
        // --- FIX START: Check ALL active columns, not just the first one ---
        bool okX, okY, okZ;
        double valX = t[params_.eastingColumn].toDouble(&okX);
        double valY = t[params_.northingColumn].toDouble(&okY);
        double valZ = t[params_.topoColumn].toDouble(&okZ);

        // If ANY of the required columns are not numbers, assumes it's a header or bad line
        if (!okX || !okY || !okZ) {
            if (lineIndex < 5) { // Only log this for the first few lines to avoid log spam
                emit logMessage(QString("Skipping non-numeric line %1: %2...").arg(lineIndex).arg(line.left(20)));
            }
            lineIndex++;
            continue; 
        }
        // --- FIX END ---

        mesh_.x.push_back(valX);
        mesh_.y.push_back(valY);
        mesh_.elevation.push_back(valZ);
        lineIndex++;
    }

    emit logMessage(QString("Loaded %1 valid input points").arg(mesh_.x.size()));
    return !mesh_.x.empty();
}

// ============================================================
// Extents
// ============================================================

void ModelGenerator::computeExtent() {
    if (mesh_.x.empty()) return;

    mesh_.xMin = mesh_.xMax = mesh_.x[0];
    mesh_.yMin = mesh_.yMax = mesh_.y[0];
    mesh_.zMin = mesh_.zMax = mesh_.elevation[0];

    for (size_t i = 1; i < mesh_.x.size(); ++i) {
        mesh_.xMin = std::min(mesh_.xMin, mesh_.x[i]);
        mesh_.xMax = std::max(mesh_.xMax, mesh_.x[i]);
        mesh_.yMin = std::min(mesh_.yMin, mesh_.y[i]);
        mesh_.yMax = std::max(mesh_.yMax, mesh_.y[i]);
        mesh_.zMin = std::min(mesh_.zMin, mesh_.elevation[i]);
        mesh_.zMax = std::max(mesh_.zMax, mesh_.elevation[i]);
    }
}

// ============================================================
// Surface grid
// ============================================================

void ModelGenerator::buildSurfaceGrid() {
    mesh_.nx = static_cast<int>((mesh_.xMax - mesh_.xMin) / params_.dx) + 1;
    mesh_.ny = static_cast<int>((mesh_.yMax - mesh_.yMin) / params_.dy) + 1;
    mesh_.nSurfaceVertices = mesh_.nx * mesh_.ny;
    mesh_.nLayers = static_cast<int>(params_.dzLayers.size());
}

// ============================================================
// Layered vertices
// ============================================================

void ModelGenerator::buildLayeredVertices() {
    mesh_.vertices.clear();
    mesh_.vertices.reserve(mesh_.nSurfaceVertices * (mesh_.nLayers + 1));

    auto nearestZ = [&](double xq, double yq) {
        double best = 1e30, z = mesh_.zMin;
        for (size_t i = 0; i < mesh_.x.size(); ++i) {
            double dx = mesh_.x[i] - xq;
            double dy = mesh_.y[i] - yq;
            double d = dx*dx + dy*dy;
            if (d < best) { best = d; z = mesh_.elevation[i]; }
        }
        return z;
    };

    // Surface layer (layer 0)
    for (int iy = 0; iy < mesh_.ny; ++iy) {
        for (int ix = 0; ix < mesh_.nx; ++ix) {
            double x = mesh_.xMin + ix * params_.dx;
            double y = mesh_.yMin + iy * params_.dy;
            double z = params_.topoAware ? nearestZ(x,y) : mesh_.zMax;
            mesh_.vertices.push_back({x,y,z});
        }
    }

    // Subsurface layers
    double depth = 0.0;
    for (int l = 0; l < mesh_.nLayers; ++l) {
        depth += params_.dzLayers[l];
        for (int i = 0; i < mesh_.nSurfaceVertices; ++i) {
            Vertex v = mesh_.vertices[i];
            v.z -= depth;
            mesh_.vertices.push_back(v);
        }
    }
}

// ============================================================
// Prism construction
// ============================================================

void ModelGenerator::buildPrisms() {
    mesh_.prisms.clear();
    int prismID = 0;

    auto vid = [&](int layer, int ix, int iy) {
        return layer * mesh_.nSurfaceVertices + iy * mesh_.nx + ix;
    };

    for (int l = 0; l < mesh_.nLayers; ++l) {
        for (int iy = 0; iy < mesh_.ny - 1; ++iy) {
            for (int ix = 0; ix < mesh_.nx - 1; ++ix) {

                int v00 = vid(l, ix, iy);
                int v10 = vid(l, ix+1, iy);
                int v01 = vid(l, ix, iy+1);
                int v11 = vid(l, ix+1, iy+1);

                int t00 = vid(l+1, ix, iy);
                int t10 = vid(l+1, ix+1, iy);
                int t01 = vid(l+1, ix, iy+1);
                int t11 = vid(l+1, ix+1, iy+1);

                // Triangle 1
                TriangularPrism p1;
                p1.id = prismID++;
                p1.layer = l;
                p1.vertexIDs = {v00, v10, v11, t00, t10, t11};
                mesh_.prisms.push_back(p1);

                // Triangle 2
                TriangularPrism p2;
                p2.id = prismID++;
                p2.layer = l;
                p2.vertexIDs = {v00, v11, v01, t00, t11, t01};
                mesh_.prisms.push_back(p2);
            }
        }
    }
}

// ============================================================
// Geometry
// ============================================================

void ModelGenerator::computePrismGeometry() {
    for (auto& p : mesh_.prisms) {
        computeCentroid(p);
        p.volume = computePrismVolume(p);
    }
}

void ModelGenerator::computeCentroid(TriangularPrism& p) const {
    double cx=0, cy=0, cz=0;
    for (int i=0;i<6;++i) {
        const Vertex& v = mesh_.vertices[p.vertexIDs[i]];
        cx+=v.x; cy+=v.y; cz+=v.z;
    }
    p.centroid = {cx/6, cy/6, cz/6};
}

double ModelGenerator::computePrismVolume(const TriangularPrism& p) const {
    const Vertex& a = mesh_.vertices[p.vertexIDs[0]];
    const Vertex& b = mesh_.vertices[p.vertexIDs[1]];
    const Vertex& c = mesh_.vertices[p.vertexIDs[2]];
    const Vertex& d = mesh_.vertices[p.vertexIDs[3]];

    double area = 0.5 * std::abs(
        (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x));
    double height = std::abs(a.z - d.z);
    return area * height;
}

// ============================================================
// Fast adjacency (NO STRINGS)
// ============================================================

struct FaceKey {
    int v[4];
    bool operator==(const FaceKey& o) const {
        for (int i=0;i<4;++i) if (v[i]!=o.v[i]) return false;
        return true;
    }
};

namespace std {
template<>
struct hash<FaceKey> {
    size_t operator()(const FaceKey& f) const {
        size_t h=0;
        for(int i=0;i<4;++i)
            h ^= std::hash<int>()(f.v[i]+0x9e3779b9+(h<<6)+(h>>2));
        return h;
    }
};
}

void ModelGenerator::buildAdjacency() {
    std::unordered_map<FaceKey,int> faceOwner;

    auto addFace = [&](std::array<int,4> verts, int pid) {
        std::sort(verts.begin(), verts.end());
        FaceKey key{{verts[0],verts[1],verts[2],verts[3]}};
        auto it = faceOwner.find(key);
        if (it==faceOwner.end()) faceOwner[key]=pid;
        else {
            mesh_.prisms[pid].neighbors.push_back(it->second);
            mesh_.prisms[it->second].neighbors.push_back(pid);
        }
    };

    for (const auto& p : mesh_.prisms) {
        addFace({p.vertexIDs[0],p.vertexIDs[1],p.vertexIDs[2],-1},p.id);
        addFace({p.vertexIDs[3],p.vertexIDs[4],p.vertexIDs[5],-1},p.id);
        addFace({p.vertexIDs[0],p.vertexIDs[1],p.vertexIDs[4],p.vertexIDs[3]},p.id);
        addFace({p.vertexIDs[1],p.vertexIDs[2],p.vertexIDs[5],p.vertexIDs[4]},p.id);
        addFace({p.vertexIDs[2],p.vertexIDs[0],p.vertexIDs[3],p.vertexIDs[5]},p.id);
    }
    
    // NEW: Compute shared face areas for TomoFast
    emit logMessage("Computing shared face areas...");
    for (auto& p : mesh_.prisms) {
        p.sharedFaceArea.resize(p.neighbors.size());
        
        for (size_t i = 0; i < p.neighbors.size(); ++i) {
            int neighborID = p.neighbors[i];
            p.sharedFaceArea[i] = computeSharedFaceArea(p, mesh_.prisms[neighborID]);
        }
    }
    emit logMessage(QString("Computed face areas for %1 adjacencies").arg(mesh_.prisms.size()));
}

// ============================================================
// NEW: Shared Face Area Computation (for TomoFast)
// ============================================================

double ModelGenerator::computeSharedFaceArea(
    const TriangularPrism& p1, 
    const TriangularPrism& p2) const 
{
    // Find shared vertices between two prisms
    std::vector<int> shared;
    for (int v1 : p1.vertexIDs) {
        for (int v2 : p2.vertexIDs) {
            if (v1 == v2) {
                shared.push_back(v1);
            }
        }
    }
    
    if (shared.size() == 3) {
        // Triangular face (top or bottom)
        const Vertex& a = mesh_.vertices[shared[0]];
        const Vertex& b = mesh_.vertices[shared[1]];
        const Vertex& c = mesh_.vertices[shared[2]];
        
        // Cross product for triangle area
        double dx1 = b.x - a.x, dy1 = b.y - a.y, dz1 = b.z - a.z;
        double dx2 = c.x - a.x, dy2 = c.y - a.y, dz2 = c.z - a.z;
        
        double cx = dy1*dz2 - dz1*dy2;
        double cy = dz1*dx2 - dx1*dz2;
        double cz = dx1*dy2 - dy1*dx2;
        
        return 0.5 * std::sqrt(cx*cx + cy*cy + cz*cz);
    }
    else if (shared.size() == 4) {
        // Quadrilateral face (side face) - split into 2 triangles
        const Vertex& a = mesh_.vertices[shared[0]];
        const Vertex& b = mesh_.vertices[shared[1]];
        const Vertex& c = mesh_.vertices[shared[2]];
        const Vertex& d = mesh_.vertices[shared[3]];
        
        // Triangle 1: a, b, c
        double dx1 = b.x - a.x, dy1 = b.y - a.y, dz1 = b.z - a.z;
        double dx2 = c.x - a.x, dy2 = c.y - a.y, dz2 = c.z - a.z;
        double cx1 = dy1*dz2 - dz1*dy2;
        double cy1 = dz1*dx2 - dx1*dz2;
        double cz1 = dx1*dy2 - dy1*dx2;
        double area1 = 0.5 * std::sqrt(cx1*cx1 + cy1*cy1 + cz1*cz1);
        
        // Triangle 2: a, c, d
        dx1 = c.x - a.x; dy1 = c.y - a.y; dz1 = c.z - a.z;
        dx2 = d.x - a.x; dy2 = d.y - a.y; dz2 = d.z - a.z;
        double cx2 = dy1*dz2 - dz1*dy2;
        double cy2 = dz1*dx2 - dx1*dz2;
        double cz2 = dx1*dy2 - dy1*dx2;
        double area2 = 0.5 * std::sqrt(cx2*cx2 + cy2*cy2 + cz2*cz2);
        
        return area1 + area2;
    }
    
    return 0.0; // No shared face (shouldn't happen)
}

// ============================================================
// Output
// ============================================================

bool ModelGenerator::writeMeshFiles() const {
    // FIX: Using correct Enum names from header
    if (params_.outputFormat == Parameters::POLYTOMO_PMESH || params_.outputFormat == Parameters::BOTH)
        if (!writePMesh(params_.outputFile)) return false;

    if (params_.outputFormat == Parameters::VTK_VTU || params_.outputFormat == Parameters::BOTH) {
        QString vtuFile = params_.outputFile;
        if (!vtuFile.endsWith(".vtu")) {
            vtuFile += ".vtu";
        }
        if (!writeVTK(vtuFile)) return false;
    }

    return true;
}

bool ModelGenerator::writePMesh(const QString& filename) const {
    QFile f(filename);
    if (!f.open(QIODevice::WriteOnly|QIODevice::Text)) return false;
    QTextStream out(&f);

    out << mesh_.vertices.size() << " " << mesh_.prisms.size() << "\n";
    for (const auto& v : mesh_.vertices)
        out << v.x<<" "<<v.y<<" "<<v.z<<"\n";

    for (const auto& p : mesh_.prisms) {
        for (int id : p.vertexIDs) out << id << " ";
        out << "\n";
    }
    return true;
}

bool ModelGenerator::writeVTK(const QString& filename) const {
    QFile f(filename);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
    QTextStream o(&f);

    o << "<?xml version=\"1.0\"?>\n";
    o << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    o << "<UnstructuredGrid>\n";
    o << "<Piece NumberOfPoints=\"" << mesh_.vertices.size()
      << "\" NumberOfCells=\"" << mesh_.prisms.size() << "\">\n";

    // ------------------ Points ------------------
    o << "<Points>\n";
    o << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& v : mesh_.vertices)
        o << v.x << " " << v.y << " " << v.z << "\n";
    o << "</DataArray>\n";
    o << "</Points>\n";

    // ------------------ Cells ------------------
    o << "<Cells>\n";

    // Connectivity
    o << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& p : mesh_.prisms) {
        for (int i = 0; i < 6; ++i)
            o << p.vertexIDs[i] << " ";
        o << "\n";
    }
    o << "</DataArray>\n";

    // Offsets
    o << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (size_t i = 0; i < mesh_.prisms.size(); ++i) {
        offset += 6;
        o << offset << "\n";
    }
    o << "</DataArray>\n";

    // Types (VTK_WEDGE = 13)
    o << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < mesh_.prisms.size(); ++i)
        o << "13\n";
    o << "</DataArray>\n";

    o << "</Cells>\n";

    // ------------------ Cell Data ------------------
    o << "<CellData Scalars=\"layer_id\">\n";

    // Layer index
    o << "<DataArray type=\"Int32\" Name=\"layer_id\" format=\"ascii\">\n";
    for (const auto& p : mesh_.prisms)
        o << p.layer << "\n";
    o << "</DataArray>\n";

    // Depth (calculated from actual vertex positions)
    o << "<DataArray type=\"Float64\" Name=\"depth\" format=\"ascii\">\n";
    for (const auto& p : mesh_.prisms) {
        // Get Z coordinates of all 6 vertices
        double z_sum = 0.0;
        double z_min = mesh_.vertices[p.vertexIDs[0]].z;
        double z_max = z_min;
        
        for (int vid : p.vertexIDs) {
            double z = mesh_.vertices[vid].z;
            z_sum += z;
            if (z < z_min) z_min = z;
            if (z > z_max) z_max = z;
        }
        
        double z_avg = z_sum / 6.0;
        
        // For topography-aware: depth below local surface
        // For flat: simple negative of average Z
        double depth;
        if (params_.topoAware) {
            // Top of cell is at z_max, calculate depth below it
            depth = z_max - z_avg;
        } else {
            // Flat mesh: Z is negative below surface
            depth = -z_avg;
        }
        
        // Safety: ensure non-negative
        if (depth < 0) depth = 0;
        
        o << depth << "\n";
    }
    o << "</DataArray>\n";

    // Active flag
    o << "<DataArray type=\"Int32\" Name=\"active\" format=\"ascii\">\n";
    for (const auto& p : mesh_.prisms)
        o << (p.active ? 1 : 0) << "\n";
    o << "</DataArray>\n";

    o << "</CellData>\n";

    o << "</Piece>\n";
    o << "</UnstructuredGrid>\n";
    o << "</VTKFile>\n";

    return true;
}

// ============================================================
// Validation
// ============================================================

bool ModelGenerator::validateMesh() const {
    for (const auto& p : mesh_.prisms)
        if (p.volume<=0) return false;
    return true;
}

// ============================================================
// Estimation Helpers
// ============================================================

long long ModelGenerator::getEstimatedVertices() const {
    if (mesh_.nx == 0 || mesh_.ny == 0 || mesh_.nLayers == 0) {
        // If mesh not built, estimate from parameters
        if (!mesh_.x.empty() && !params_.dzLayers.empty()) {
            double rangeX = mesh_.xMax - mesh_.xMin;
            double rangeY = mesh_.yMax - mesh_.yMin;
            int nx = static_cast<int>(rangeX / params_.dx) + 1;
            int ny = static_cast<int>(rangeY / params_.dy) + 1;
            int nz = static_cast<int>(params_.dzLayers.size());
            return (long long)(nx + 1) * (ny + 1) * (nz + 1);
        }
        return 0;
    }
    return (long long)(mesh_.nx + 1) * (mesh_.ny + 1) * (mesh_.nLayers + 1);
}

long long ModelGenerator::getEstimatedPrisms() const {
    if (mesh_.nx == 0 || mesh_.ny == 0 || mesh_.nLayers == 0) {
        // If mesh not built, estimate from parameters
        if (!mesh_.x.empty() && !params_.dzLayers.empty()) {
            double rangeX = mesh_.xMax - mesh_.xMin;
            double rangeY = mesh_.yMax - mesh_.yMin;
            int nx = static_cast<int>(rangeX / params_.dx) + 1;
            int ny = static_cast<int>(rangeY / params_.dy) + 1;
            int nz = static_cast<int>(params_.dzLayers.size());
            return (long long)nx * ny * nz * 2;
        }
        return 0;
    }
    return (long long)mesh_.nx * mesh_.ny * mesh_.nLayers * 2;
}

double ModelGenerator::getEstimatedMemoryMB() const {
    long long v = getEstimatedVertices();
    long long p = getEstimatedPrisms();
    
    // Rough estimate:
    // Vertex: 3 doubles (24 bytes)
    // Prism: 6 ints (24 bytes) + neighbors (~32 bytes) + sharedFaceArea (~32 bytes) + overhead
    double bytes = v * 24.0 + p * 100.0; 
    return bytes / (1024.0 * 1024.0);
}
