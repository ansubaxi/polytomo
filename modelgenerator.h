#ifndef MODELGENERATOR_H
#define MODELGENERATOR_H

#include <QObject>
#include <QString>
#include <vector>
#include <array>

class ModelGenerator : public QObject {
    Q_OBJECT

public:
    // Parameters Structure
    struct Parameters {
        QString inputFile;
        QString outputFile;
        double dx = 10.0;
        double dy = 10.0;
        double gridOffset = 0.0;
        double depthBelowSurface = 1000.0;
        std::vector<double> dzLayers;
        int eastingColumn = 0;
        int northingColumn = 1;
        int topoColumn = 2;
        bool topoAware = true;
        
        enum OutputFormat {
            POLYTOMO_PMESH,
            VTK_VTU,
            BOTH
        } outputFormat = POLYTOMO_PMESH;
    };

    // Mesh Data Structures
    struct Vertex {
        double x, y, z;
    };

    struct TriangularPrism {
        int id;
        int layer;
        std::array<int, 6> vertexIDs;
        double volume;
        std::array<double, 3> centroid;
        std::vector<int> neighbors;
        std::vector<double> sharedFaceArea;  // FIXED: Added for TomoFast integration
        bool active = true;
    };

    struct MeshData {
        std::vector<double> x, y, elevation;
        double xMin = 0, xMax = 0;
        double yMin = 0, yMax = 0;
        double zMin = 0, zMax = 0;
        int nx = 0, ny = 0, nLayers = 0;
        int nSurfaceVertices = 0;  // FIXED: Added for vertex indexing
        std::vector<Vertex> vertices;
        std::vector<TriangularPrism> prisms;
    };

    explicit ModelGenerator(QObject* parent = nullptr);
    ~ModelGenerator();

    void setParameters(const Parameters& params);
    const Parameters& parameters() const;  // FIXED: Added declaration
    const MeshData& mesh() const;
    
    // Call after setParameters to compute extents for UI display
    bool computeExtentFromParameters();
    
    // Getters for mainwindow
    int getNx() const { return mesh_.nx; }
    int getNy() const { return mesh_.ny; }
    int getNz() const { return mesh_.nLayers; }
    double getXMin() const { return mesh_.xMin; }
    double getXMax() const { return mesh_.xMax; }
    double getYMin() const { return mesh_.yMin; }
    double getYMax() const { return mesh_.yMax; }
    double getMinElevation() const { return mesh_.zMin; }
    double getMaxElevation() const { return mesh_.zMax; }
    
    // Estimation helpers (called before process())
    long long getEstimatedVertices() const;
    long long getEstimatedPrisms() const;
    double getEstimatedMemoryMB() const;

public slots:
    void process();

signals:
    void progressUpdated(int percent);
    void logMessage(const QString& msg);
    void finished(bool success, const QString& message);  // FIXED: Consistent 2-param signature

private:
    bool readInputFile();
    void computeExtent();
    void buildSurfaceGrid();
    void buildLayeredVertices();
    void buildPrisms();
    void computePrismGeometry();
    void computeCentroid(TriangularPrism& p) const;
    double computePrismVolume(const TriangularPrism& p) const;
    void buildAdjacency();
    double computeSharedFaceArea(const TriangularPrism& p1, const TriangularPrism& p2) const;  // NEW: For TomoFast
    bool writeMeshFiles() const;
    bool writePMesh(const QString& filename) const;
    bool writeVTK(const QString& filename) const;
    bool validateMesh() const;

    Parameters params_;
    MeshData mesh_;
};

#endif // MODELGENERATOR_H
