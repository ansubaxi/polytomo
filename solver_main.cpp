/**
 * solver_main.cpp
 * 
 * MPI-enabled command-line solver for PolyTomo mesh generation
 * 
 * Usage:
 *   Single process:  PolyTomoSolver.exe parameters.json
 *   Multi-process:   mpiexec -n 128 PolyTomoSolver.exe parameters.json
 */

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

// We'll use a simplified version of ModelGenerator without Qt dependencies
// This is a standalone MPI solver

struct Parameters {
    std::string inputFile;
    std::string outputFile;
    double dx;
    double dy;
    double gridOffset;
    double depthBelowSurface;
    std::vector<double> dzLayers;
    bool topoAware;
    int eastingColumn;
    int northingColumn;
    int topoColumn;
    
    // Output format: 0=pmesh, 1=vtu, 2=both
    int outputFormat;
    
    // Calculated
    int nx, ny, nz;
    double xMin, xMax, yMin, yMax;
    double maxElevation;
};

struct Vertex {
    double x, y, z;
    long long globalID;
};

struct TriangularPrism {
    long long vertices[6];
    long long cellID;
    int layer_id;           // Layer index (0, 1, 2, ...)
    double depth;           // Depth at center of cell (meters below surface)
};

struct Point3D {
    double x, y, z;
};

// Function declarations
bool readParametersJSON(const std::string& filename, Parameters& params);
bool readInputFile(const std::string& filename, std::vector<Point3D>& dataPoints, Parameters& params);
void calculateGridExtent(const std::vector<Point3D>& dataPoints, Parameters& params);
void calculateGridDimensions(Parameters& params);
double interpolateElevation(double x, double y, const std::vector<Point3D>& dataPoints, const Parameters& params);
void calculateMPIDomain(int rank, int size, int ny_total, int& start_y, int& end_y, int& count_y);
bool generateMeshMPI(int rank, int size, const Parameters& params, const std::vector<Point3D>& dataPoints);
bool writePolyTomoFormatMPI(const std::string& filename, int rank, int size,
                            const std::vector<Vertex>& localVertices,
                            const std::vector<TriangularPrism>& localPrisms,
                            long long globalVertexStartID,
                            long long globalCellStartID,
                            long long totalVertices,
                            long long totalCells);
bool writeVTKFormatMPI(const std::string& filename, int rank, int size,
                       const std::vector<Vertex>& localVertices,
                       const std::vector<TriangularPrism>& localPrisms,
                       long long globalVertexStartID,
                       long long globalCellStartID,
                       long long totalVertices,
                       long long totalCells);

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Check arguments
    if (argc < 2) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <parameters.json>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }
    
    std::string paramFile = argv[1];
    
    if (rank == 0) {
        std::cout << "===============================================" << std::endl;
        std::cout << "PolyTomo MPI Mesh Generator" << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << "MPI processes: " << size << std::endl;
        std::cout << "Parameter file: " << paramFile << std::endl;
        std::cout << "===============================================" << std::endl;
    }
    
    // Read parameters (all ranks read independently for simplicity)
    Parameters params;
    if (!readParametersJSON(paramFile, params)) {
        if (rank == 0) {
            std::cerr << "ERROR: Failed to read parameters file" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }
    
    // Read input data (all ranks read - could optimize with broadcast)
    std::vector<Point3D> dataPoints;
    if (!readInputFile(params.inputFile, dataPoints, params)) {
        if (rank == 0) {
            std::cerr << "ERROR: Failed to read input file: " << params.inputFile << std::endl;
        }
        MPI_Finalize();
        return 1;
    }
    
    if (rank == 0) {
        std::cout << "Input data: " << dataPoints.size() << " points" << std::endl;
        std::cout << "Grid extent: X=[" << params.xMin << ", " << params.xMax << "]" << std::endl;
        std::cout << "             Y=[" << params.yMin << ", " << params.yMax << "]" << std::endl;
        std::cout << "Grid size: " << params.nx << " x " << params.ny << " x " << params.nz << std::endl;
        std::cout << "Max elevation: " << params.maxElevation << " m" << std::endl;
        std::cout << "===============================================" << std::endl;
    }
    
    // Generate mesh in parallel
    bool success = generateMeshMPI(rank, size, params, dataPoints);
    
    if (rank == 0) {
        if (success) {
            std::cout << "===============================================" << std::endl;
            std::cout << "SUCCESS: Mesh generation completed!" << std::endl;
            std::cout << "Output file: " << params.outputFile << std::endl;
            std::cout << "===============================================" << std::endl;
        } else {
            std::cout << "ERROR: Mesh generation failed!" << std::endl;
        }
    }
    
    MPI_Finalize();
    return success ? 0 : 1;
}

bool readParametersJSON(const std::string& filename, Parameters& params) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    // Read entire file into string
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string content = buffer.str();
    
    // Remove whitespace/newlines for easier parsing
    content.erase(std::remove(content.begin(), content.end(), '\n'), content.end());
    content.erase(std::remove(content.begin(), content.end(), '\r'), content.end());
    content.erase(std::remove(content.begin(), content.end(), '\t'), content.end());
    
    // Helper functions to extract values
    auto extractString = [&](const std::string& key) -> std::string {
        size_t pos = content.find("\"" + key + "\":");
        if (pos == std::string::npos) return "";
        pos = content.find('"', pos + key.length() + 3);
        if (pos == std::string::npos) return "";
        size_t end = content.find('"', pos + 1);
        return content.substr(pos + 1, end - pos - 1);
    };
    
    auto extractDouble = [&](const std::string& key) -> double {
        size_t pos = content.find("\"" + key + "\":");
        if (pos == std::string::npos) return 0.0;
        pos += key.length() + 3;
        while (pos < content.length() && (content[pos] == ' ' || content[pos] == ':')) pos++;
        return std::stod(content.substr(pos));
    };
    
    auto extractInt = [&](const std::string& key) -> int {
        size_t pos = content.find("\"" + key + "\":");
        if (pos == std::string::npos) return 0;
        pos += key.length() + 3;
        while (pos < content.length() && (content[pos] == ' ' || content[pos] == ':')) pos++;
        return std::stoi(content.substr(pos));
    };
    
    auto extractBool = [&](const std::string& key) -> bool {
        size_t pos = content.find("\"" + key + "\":");
        if (pos == std::string::npos) return false;
        return content.find("true", pos) != std::string::npos;
    };
    
    // Extract all parameters
    params.inputFile = extractString("inputFile");
    params.outputFile = extractString("outputFile");
    params.dx = extractDouble("dx");
    params.dy = extractDouble("dy");
    params.gridOffset = extractDouble("gridOffset");
    params.depthBelowSurface = extractDouble("depthBelowSurface");
    params.topoAware = extractBool("topoAware");
    params.eastingColumn = extractInt("eastingColumn");
    params.northingColumn = extractInt("northingColumn");
    params.topoColumn = extractInt("topoColumn");
    params.outputFormat = extractInt("outputFormat");
    
    // Set default output format if not specified
    if (params.outputFormat < 0 || params.outputFormat > 2) {
        params.outputFormat = 0;
    }
    
    // === CRITICAL FIX: Properly parse layers array ===
    params.dzLayers.clear();
    size_t layersPos = content.find("\"layers\":");
    if (layersPos != std::string::npos) {
        size_t bracketStart = content.find('[', layersPos);
        size_t bracketEnd = content.find(']', bracketStart);
        
        if (bracketStart != std::string::npos && bracketEnd != std::string::npos) {
            std::string layersStr = content.substr(bracketStart + 1, bracketEnd - bracketStart - 1);
            // Remove any remaining whitespace
            layersStr.erase(std::remove(layersStr.begin(), layersStr.end(), ' '), layersStr.end());
            
            // Parse comma-separated numbers
            std::stringstream ss(layersStr);
            std::string token;
            while (std::getline(ss, token, ',')) {
                if (!token.empty()) {
                    try {
                        double dz = std::stod(token);
                        params.dzLayers.push_back(dz);
                    } catch (...) {
                        // Skip invalid values
                    }
                }
            }
        }
    }
    
    // === CRITICAL FIX: Fallback if no layers defined ===
    if (params.dzLayers.empty() && params.depthBelowSurface > 0) {
        std::cerr << "WARNING: No layers in JSON! Creating 1 default layer." << std::endl;
        params.dzLayers.push_back(params.depthBelowSurface);
    }
    
    params.nz = params.dzLayers.size();
    
    // Debug output
    std::cout << "Parsed " << params.dzLayers.size() << " layers: ";
    for (size_t i = 0; i < std::min(size_t(5), params.dzLayers.size()); ++i) {
        std::cout << params.dzLayers[i] << "m ";
    }
    if (params.dzLayers.size() > 5) std::cout << "...";
    std::cout << std::endl;
    
    return true;
}

bool readInputFile(const std::string& filename, std::vector<Point3D>& dataPoints, Parameters& params) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<double> values;
        
        while (std::getline(ss, token, ',')) {
            try {
                values.push_back(std::stod(token));
            } catch (...) {
                continue;
            }
        }
        
        if (values.size() >= 3) {
            Point3D p;
            p.x = values[params.eastingColumn];
            p.y = values[params.northingColumn];
            p.z = values[params.topoColumn];
            dataPoints.push_back(p);
        }
    }
    
    if (dataPoints.empty()) {
        return false;
    }
    
    // Calculate grid extent and dimensions
    calculateGridExtent(dataPoints, params);
    calculateGridDimensions(params);
    
    return true;
}

void calculateGridExtent(const std::vector<Point3D>& dataPoints, Parameters& params) {
    params.xMin = params.xMax = dataPoints[0].x;
    params.yMin = params.yMax = dataPoints[0].y;
    params.maxElevation = dataPoints[0].z;
    
    for (const auto& p : dataPoints) {
        if (p.x < params.xMin) params.xMin = p.x;
        if (p.x > params.xMax) params.xMax = p.x;
        if (p.y < params.yMin) params.yMin = p.y;
        if (p.y > params.yMax) params.yMax = p.y;
        if (p.z > params.maxElevation) params.maxElevation = p.z;
    }
}

void calculateGridDimensions(Parameters& params) {
    double dx_extent = params.xMax - params.xMin;
    double dy_extent = params.yMax - params.yMin;
    
    // Add offset to avoid edge points
    double offset_x = dx_extent * params.gridOffset;
    double offset_y = dy_extent * params.gridOffset;
    
    params.xMin -= offset_x;
    params.xMax += offset_x;
    params.yMin -= offset_y;
    params.yMax += offset_y;
    
    dx_extent = params.xMax - params.xMin;
    dy_extent = params.yMax - params.yMin;
    
    params.nx = static_cast<int>(std::ceil(dx_extent / params.dx));
    params.ny = static_cast<int>(std::ceil(dy_extent / params.dy));
}

double interpolateElevation(double x, double y, const std::vector<Point3D>& dataPoints, const Parameters& params) {
    // Simple nearest neighbor interpolation
    // For production, use IDW or kriging
    
    double minDist = 1e30;
    double elevation = params.maxElevation;
    
    for (const auto& p : dataPoints) {
        double dx = x - p.x;
        double dy = y - p.y;
        double dist = std::sqrt(dx*dx + dy*dy);
        
        if (dist < minDist) {
            minDist = dist;
            elevation = p.z;
        }
    }
    
    return elevation;
}

void calculateMPIDomain(int rank, int size, int ny_total, int& start_y, int& end_y, int& count_y) {
    int rows_per_rank = ny_total / size;
    int remainder = ny_total % size;
    
    start_y = rank * rows_per_rank + std::min(rank, remainder);
    count_y = rows_per_rank + (rank < remainder ? 1 : 0);
    end_y = start_y + count_y;
}

bool generateMeshMPI(int rank, int size, const Parameters& params, const std::vector<Point3D>& dataPoints) {
    // Calculate domain decomposition
    int local_start_y, local_end_y, local_ny_count;
    calculateMPIDomain(rank, size, params.ny, local_start_y, local_end_y, local_ny_count);
    
    if (rank == 0) {
        std::cout << "Starting parallel mesh generation..." << std::endl;
    }
    std::cout << "Rank " << rank << " processing Y-rows " << local_start_y 
              << " to " << local_end_y << " (" << local_ny_count << " rows)" << std::endl;
    
    // Calculate cumulative depths
    std::vector<double> cumulativeDepths;
    cumulativeDepths.push_back(0.0);
    for (const auto& dz : params.dzLayers) {
        cumulativeDepths.push_back(cumulativeDepths.back() + dz);
    }
    
    // Calculate number of vertices this rank will generate
    int rows_of_vertices = local_ny_count;
    if (rank == size - 1) rows_of_vertices++; // Last rank includes final edge
    
    long long vertices_per_depth_layer = (long long)(params.nx + 1) * rows_of_vertices;
    long long local_num_vertices = vertices_per_depth_layer * (params.nz + 1);
    
    // Calculate global vertex start ID using MPI_Exscan
    long long globalVertexStartID = 0;
    MPI_Exscan(&local_num_vertices, &globalVertexStartID, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) globalVertexStartID = 0;
    
    // Calculate number of cells this rank will generate
    long long local_num_cells = (long long)params.nx * local_ny_count * params.nz * 2;
    
    // Calculate global cell start ID
    long long globalCellStartID = 0;
    MPI_Exscan(&local_num_cells, &globalCellStartID, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) globalCellStartID = 0;
    
    std::cout << "Rank " << rank << " generating " << local_num_vertices << " vertices, " 
              << local_num_cells << " cells" << std::endl;
    
    // Generate local vertices
    std::vector<Vertex> localVertices;
    localVertices.reserve(local_num_vertices);
    
    long long currentVID = globalVertexStartID;
    
    for (int k = 0; k <= params.nz; ++k) {
        double depth = (k < cumulativeDepths.size()) ? cumulativeDepths[k] : params.depthBelowSurface;
        
        for (int j = 0; j < rows_of_vertices; ++j) {
            int gridJ = local_start_y + j;
            double y = params.yMin + gridJ * params.dy;
            
            for (int i = 0; i <= params.nx; ++i) {
                double x = params.xMin + i * params.dx;
                double elevation = interpolateElevation(x, y, dataPoints, params);
                
                Vertex v;
                v.x = x;
                v.y = y;
                v.z = params.topoAware ? (elevation - depth) : -depth;
                v.globalID = currentVID++;
                
                localVertices.push_back(v);
            }
        }
    }
    
    std::cout << "Rank " << rank << " generated vertices" << std::endl;
    
    // Generate local prisms
    std::vector<TriangularPrism> localPrisms;
    localPrisms.reserve(local_num_cells);
    
    long long currentCID = globalCellStartID;
    
    // Helper function to calculate global vertex ID
    auto mathGlobalID = [&](int i, int global_j, int k) -> long long {
        long long nodes_per_layer = (long long)(params.ny + 1) * (params.nx + 1);
        long long nodes_per_row = params.nx + 1;
        return (long long)k * nodes_per_layer + (long long)global_j * nodes_per_row + i;
    };
    
    for (int k = 0; k < params.nz; ++k) {
        for (int j = 0; j < local_ny_count; ++j) {
            int global_j = local_start_y + j;
            
            for (int i = 0; i < params.nx; ++i) {
                // Get cell center for elevation lookup
                double cell_x = params.xMin + (i + 0.5) * params.dx;
                double cell_y = params.yMin + (global_j + 0.5) * params.dy;
                double surface_elev = interpolateElevation(cell_x, cell_y, dataPoints, params);
                
                // Calculate actual Z coordinates at layer boundaries
                double depth_top = cumulativeDepths[k];
                double depth_bottom = (k + 1 < cumulativeDepths.size()) ? 
                                     cumulativeDepths[k + 1] : params.depthBelowSurface;
                
                double z_top = params.topoAware ? (surface_elev - depth_top) : -depth_top;
                double z_bottom = params.topoAware ? (surface_elev - depth_bottom) : -depth_bottom;
                double z_center = (z_top + z_bottom) / 2.0;
                
                // CRITICAL FIX: Depth below MAXIMUM elevation (not local surface!)
                // This ensures proper vertical color gradient in visualization
                double cell_depth = params.maxElevation - z_center;
                
                // Calculate 8 corner vertices of hexahedron
                long long v0 = mathGlobalID(i, global_j, k);
                long long v1 = mathGlobalID(i + 1, global_j, k);
                long long v2 = mathGlobalID(i + 1, global_j + 1, k);
                long long v3 = mathGlobalID(i, global_j + 1, k);
                
                long long v4 = mathGlobalID(i, global_j, k + 1);
                long long v5 = mathGlobalID(i + 1, global_j, k + 1);
                long long v6 = mathGlobalID(i + 1, global_j + 1, k + 1);
                long long v7 = mathGlobalID(i, global_j + 1, k + 1);
                
                // Split hexahedron into 2 triangular prisms
                // Prism 1: v0, v1, v2 (bottom) + v4, v5, v6 (top)
                TriangularPrism p1;
                p1.vertices[0] = v0;
                p1.vertices[1] = v1;
                p1.vertices[2] = v2;
                p1.vertices[3] = v4;
                p1.vertices[4] = v5;
                p1.vertices[5] = v6;
                p1.cellID = currentCID++;
                p1.layer_id = k;
                p1.depth = cell_depth;  // Depth from Z-coordinates
                localPrisms.push_back(p1);
                
                // Prism 2: v0, v2, v3 (bottom) + v4, v6, v7 (top)
                TriangularPrism p2;
                p2.vertices[0] = v0;
                p2.vertices[1] = v2;
                p2.vertices[2] = v3;
                p2.vertices[3] = v4;
                p2.vertices[4] = v6;
                p2.vertices[5] = v7;
                p2.cellID = currentCID++;
                p2.layer_id = k;
                p2.depth = cell_depth;  // Depth from Z-coordinates
                localPrisms.push_back(p2);
            }
        }
    }
    
    std::cout << "Rank " << rank << " generated prisms" << std::endl;
    
    // Calculate total vertices and cells
    long long totalVertices = (long long)(params.nx + 1) * (params.ny + 1) * (params.nz + 1);
    long long totalCells = (long long)params.nx * params.ny * params.nz * 2;
    
    // Write mesh file(s) based on output format
    bool success = true;
    
    // outputFormat: 0=pmesh, 1=vtu, 2=both
    
    // 1. Handle PMESH format
    if (params.outputFormat == 0 || params.outputFormat == 2) {
        std::string pmeshFile = params.outputFile;
        
        // Ensure proper extension (check END of filename)
        if (pmeshFile.length() < 6 || pmeshFile.substr(pmeshFile.length() - 6) != ".pmesh") {
            // Strip existing extension if any
            size_t lastDot = pmeshFile.find_last_of('.');
            if (lastDot != std::string::npos) {
                pmeshFile = pmeshFile.substr(0, lastDot);
            }
            pmeshFile += ".pmesh";
        }
        
        if (rank == 0) std::cout << "Writing PMESH to: " << pmeshFile << std::endl;
        
        if (!writePolyTomoFormatMPI(pmeshFile, rank, size, 
                                   localVertices, localPrisms,
                                   globalVertexStartID, globalCellStartID,
                                   totalVertices, totalCells)) {
            success = false;
            if (rank == 0) std::cerr << "Failed to write PMESH file." << std::endl;
        }
    }
    
    // 2. Handle VTU format
    if (params.outputFormat == 1 || params.outputFormat == 2) {
        std::string vtuFile = params.outputFile;
        
        // Ensure proper extension (check END of filename)
        if (vtuFile.length() < 4 || vtuFile.substr(vtuFile.length() - 4) != ".vtu") {
            // Strip existing extension if any
            size_t lastDot = vtuFile.find_last_of('.');
            if (lastDot != std::string::npos) {
                vtuFile = vtuFile.substr(0, lastDot);
            }
            vtuFile += ".vtu";
        }
        
        if (rank == 0) std::cout << "Writing VTU to: " << vtuFile << std::endl;
        
        if (!writeVTKFormatMPI(vtuFile, rank, size, 
                              localVertices, localPrisms,
                              globalVertexStartID, globalCellStartID,
                              totalVertices, totalCells)) {
            success = false;
            if (rank == 0) std::cerr << "Failed to write VTU file." << std::endl;
        }
    }
    
    return success;
}

bool writePolyTomoFormatMPI(const std::string& filename, int rank, int size,
                            const std::vector<Vertex>& localVertices,
                            const std::vector<TriangularPrism>& localPrisms,
                            long long globalVertexStartID,
                            long long globalCellStartID,
                            long long totalVertices,
                            long long totalCells)
{
    MPI_File fh;
    
    // Delete file if exists (Rank 0 only)
    if (rank == 0) {
        std::remove(filename.c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Open file for parallel writing
    int err = MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                            MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    if (err != MPI_SUCCESS) {
        if (rank == 0) {
            std::cerr << "MPI Error: Could not open file for writing: " << filename << std::endl;
        }
        return false;
    }
    
    // CRITICAL FIX: Truncate file and seek to start
    MPI_File_set_size(fh, 0);
    MPI_File_seek(fh, 0, MPI_SEEK_SET);
    
    // ===== STEP 1: WRITE HEADER (Rank 0 only) =====
    std::string header;
    MPI_Offset headerOffset = 0;
    
    if (rank == 0) {
        std::ostringstream oss;
        oss << "# PolyTomo Polyhedral Mesh v1.0 (MPI)\n";
        oss << "# Generated by PolyTomo MPI Mesh Generator\n";
        oss << "# Total vertices: " << totalVertices << "\n";
        oss << "# Total cells: " << totalCells << "\n";
        oss << "VERTICES " << totalVertices << "\n";
        header = oss.str();
        
        MPI_File_write(fh, header.c_str(), header.size(), MPI_CHAR, MPI_STATUS_IGNORE);
        headerOffset = header.size();
        
        std::cout << "Writing header..." << std::endl;
    }
    
    // Broadcast header size to all ranks
    MPI_Bcast(&headerOffset, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
    
    // ===== STEP 2: WRITE VERTICES (Parallel) =====
    // Prepare local vertex buffer
    std::ostringstream vBuffer;
    vBuffer.imbue(std::locale::classic());  // CRITICAL: Force dot decimal separator for coordinates!
    vBuffer << std::fixed << std::setprecision(6);
    for (const auto& v : localVertices) {
        vBuffer << v.x << " " << v.y << " " << v.z << "\n";
    }
    std::string localVBuffer = vBuffer.str();
    
    // Calculate write offset for this rank using MPI_Exscan
    MPI_Offset vLen = localVBuffer.size();
    MPI_Offset vOffset = 0;
    MPI_Exscan(&vLen, &vOffset, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) vOffset = 0;
    
    // Write vertices at calculated offset
    MPI_File_write_at_all(fh, headerOffset + vOffset, localVBuffer.c_str(), 
                          vLen, MPI_CHAR, MPI_STATUS_IGNORE);
    
    if (rank == 0) {
        std::cout << "Writing vertices..." << std::endl;
    }
    
    // ===== STEP 3: WRITE CELLS HEADER (Rank 0 only) =====
    // Calculate total size of vertices section
    MPI_Offset totalVerticesBytes = 0;
    MPI_Allreduce(&vLen, &totalVerticesBytes, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    MPI_Offset cellsHeaderStart = headerOffset + totalVerticesBytes;
    
    std::string cellsHeader;
    MPI_Offset cellsHeaderLen = 0;
    
    if (rank == 0) {
        std::ostringstream oss;
        oss << "\nCELLS " << totalCells << "\n";
        cellsHeader = oss.str();
        cellsHeaderLen = cellsHeader.size();
        
        MPI_File_write_at(fh, cellsHeaderStart, cellsHeader.c_str(), 
                         cellsHeaderLen, MPI_CHAR, MPI_STATUS_IGNORE);
        
        std::cout << "Writing cells..." << std::endl;
    }
    
    // Broadcast cells header size
    MPI_Bcast(&cellsHeaderLen, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
    
    // ===== STEP 4: WRITE CELLS (Parallel) =====
    // Prepare local cells buffer
    std::ostringstream cBuffer;
    cBuffer.imbue(std::locale::classic());  // Safety: Ensure consistent formatting
    for (const auto& p : localPrisms) {
        cBuffer << "PRISM " << p.cellID << " 2\n";
        cBuffer << "  3 " << p.vertices[0] << " " << p.vertices[1] << " " << p.vertices[2] << "\n";
        cBuffer << "  3 " << p.vertices[3] << " " << p.vertices[4] << " " << p.vertices[5] << "\n";
    }
    std::string localCBuffer = cBuffer.str();
    
    // Calculate write offset for cells
    MPI_Offset cLen = localCBuffer.size();
    MPI_Offset cOffset = 0;
    MPI_Exscan(&cLen, &cOffset, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) cOffset = 0;
    
    // Write cells at calculated offset
    MPI_File_write_at_all(fh, cellsHeaderStart + cellsHeaderLen + cOffset, 
                          localCBuffer.c_str(), cLen, MPI_CHAR, MPI_STATUS_IGNORE);
    
    // ===== STEP 5: WRITE FOOTER (Rank 0 only) =====
    // Calculate total size of cells section
    MPI_Offset totalCellsBytes = 0;
    MPI_Allreduce(&cLen, &totalCellsBytes, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    MPI_Offset footerStart = cellsHeaderStart + cellsHeaderLen + totalCellsBytes;
    
    if (rank == 0) {
        std::string footer = "\nEND\n";
        MPI_File_write_at(fh, footerStart, footer.c_str(), footer.size(), 
                         MPI_CHAR, MPI_STATUS_IGNORE);
        
        std::cout << "Writing footer..." << std::endl;
    }
    
    // Close file
    MPI_File_close(&fh);
    
    if (rank == 0) {
        std::cout << "File writing completed successfully!" << std::endl;
    }
    
    return true;
}

bool writeVTKFormatMPI(const std::string& filename, int rank, int size,
                       const std::vector<Vertex>& localVertices,
                       const std::vector<TriangularPrism>& localPrisms,
                       long long globalVertexStartID,
                       long long globalCellStartID,
                       long long totalVertices,
                       long long totalCells) {
    // Open file with MPI-IO
    MPI_File fh;
    int result = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), 
                              MPI_MODE_CREATE | MPI_MODE_WRONLY,
                              MPI_INFO_NULL, &fh);
    
    if (result != MPI_SUCCESS) {
        if (rank == 0) {
            std::cerr << "ERROR: Could not open file for writing: " << filename << std::endl;
        }
        return false;
    }
    
    // Truncate file to ensure clean start
    MPI_File_set_size(fh, 0);
    
    // CRITICAL FIX: Explicitly seek to start of file after truncation
    // Without this, file pointer might not be at position 0
    MPI_File_seek(fh, 0, MPI_SEEK_SET);
    
    // ===== STEP 1: WRITE HEADER (Rank 0 only) =====
    MPI_Offset headerLen = 0;
    if (rank == 0) {
        std::cout << "Writing VTU header..." << std::endl;
        
        std::ostringstream header;
        header.imbue(std::locale::classic());  // CRITICAL: Force dot decimal separator
        header << "<?xml version=\"1.0\"?>\n";
        header << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
        header << "  <UnstructuredGrid>\n";
        header << "    <Piece NumberOfPoints=\"" << totalVertices 
               << "\" NumberOfCells=\"" << totalCells << "\">\n";
        header << "      <Points>\n";
        header << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        
        std::string headerStr = header.str();
        headerLen = headerStr.size();
        
        // Diagnostic: Print first 100 chars of header
        std::cout << "Header first 100 chars: " << headerStr.substr(0, 100) << std::endl;
        std::cout << "Header length: " << headerLen << " bytes" << std::endl;
        
        MPI_File_write(fh, headerStr.c_str(), headerLen, MPI_CHAR, MPI_STATUS_IGNORE);
        
        // Diagnostic: Verify write position
        MPI_Offset currentPos;
        MPI_File_get_position(fh, &currentPos);
        std::cout << "File position after header write: " << currentPos << " bytes" << std::endl;
    }
    
    // Broadcast header length
    MPI_Bcast(&headerLen, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
    
    // ===== STEP 2: WRITE POINTS (Parallel) =====
    std::ostringstream vBuffer;
    vBuffer.imbue(std::locale::classic());  // CRITICAL: Force dot decimal separator for coordinates!
    for (const auto& v : localVertices) {
        vBuffer << "          " << std::fixed << std::setprecision(6) 
                << v.x << " " << v.y << " " << v.z << "\n";
    }
    std::string localVBuffer = vBuffer.str();
    
    // Calculate write offset
    MPI_Offset vLen = localVBuffer.size();
    MPI_Offset vOffset = 0;
    MPI_Exscan(&vLen, &vOffset, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) vOffset = 0;
    
    // Write vertices
    MPI_File_write_at_all(fh, headerLen + vOffset, localVBuffer.c_str(), 
                          vLen, MPI_CHAR, MPI_STATUS_IGNORE);
    
    // ===== STEP 3: WRITE POINTS FOOTER (Rank 0 only) =====
    MPI_Offset totalVerticesBytes = 0;
    MPI_Allreduce(&vLen, &totalVerticesBytes, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    MPI_Offset pointsFooterStart = headerLen + totalVerticesBytes;
    
    MPI_Offset pointsFooterLen = 0;
    if (rank == 0) {
        std::string pointsFooter = "        </DataArray>\n";
        pointsFooter += "      </Points>\n";
        pointsFooter += "      <Cells>\n";
        pointsFooter += "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
        
        pointsFooterLen = pointsFooter.size();
        MPI_File_write_at(fh, pointsFooterStart, pointsFooter.c_str(), 
                         pointsFooterLen, MPI_CHAR, MPI_STATUS_IGNORE);
        
        std::cout << "Writing VTU cells..." << std::endl;
    }
    
    MPI_Bcast(&pointsFooterLen, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
    
    // ===== STEP 4: WRITE CONNECTIVITY (Parallel) =====
    std::ostringstream connBuffer;
    connBuffer.imbue(std::locale::classic());  // Safety: Ensure consistent formatting
    for (const auto& p : localPrisms) {
        connBuffer << "          ";
        for (int i = 0; i < 6; i++) {
            connBuffer << p.vertices[i];
            if (i < 5) connBuffer << " ";
        }
        connBuffer << "\n";
    }
    std::string localConnBuffer = connBuffer.str();
    
    MPI_Offset connLen = localConnBuffer.size();
    MPI_Offset connOffset = 0;
    MPI_Exscan(&connLen, &connOffset, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) connOffset = 0;
    
    MPI_Offset connStart = pointsFooterStart + pointsFooterLen;
    MPI_File_write_at_all(fh, connStart + connOffset, localConnBuffer.c_str(), 
                          connLen, MPI_CHAR, MPI_STATUS_IGNORE);
    
    // ===== STEP 5: WRITE OFFSETS (Rank 0 computes, all write) =====
    MPI_Offset totalConnBytes = 0;
    MPI_Allreduce(&connLen, &totalConnBytes, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    MPI_Offset offsetsHeaderStart = connStart + totalConnBytes;
    
    MPI_Offset offsetsHeaderLen = 0;
    if (rank == 0) {
        std::string offsetsHeader = "        </DataArray>\n";
        offsetsHeader += "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
        
        offsetsHeaderLen = offsetsHeader.size();
        MPI_File_write_at(fh, offsetsHeaderStart, offsetsHeader.c_str(), 
                         offsetsHeaderLen, MPI_CHAR, MPI_STATUS_IGNORE);
    }
    
    MPI_Bcast(&offsetsHeaderLen, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
    
    // Write offsets (each prism has 6 vertices, so offset increases by 6)
    std::ostringstream offsetsBuffer;
    offsetsBuffer.imbue(std::locale::classic());  // Safety: Ensure consistent formatting
    for (size_t i = 0; i < localPrisms.size(); i++) {
        long long offset = (globalCellStartID + i + 1) * 6;
        offsetsBuffer << "          " << offset << "\n";
    }
    std::string localOffsetsBuffer = offsetsBuffer.str();
    
    MPI_Offset offsetsLen = localOffsetsBuffer.size();
    MPI_Offset offsetsOffset = 0;
    MPI_Exscan(&offsetsLen, &offsetsOffset, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) offsetsOffset = 0;
    
    MPI_Offset offsetsStart = offsetsHeaderStart + offsetsHeaderLen;
    MPI_File_write_at_all(fh, offsetsStart + offsetsOffset, localOffsetsBuffer.c_str(), 
                          offsetsLen, MPI_CHAR, MPI_STATUS_IGNORE);
    
    // ===== STEP 6: WRITE CELL TYPES (Rank 0 computes, all write) =====
    MPI_Offset totalOffsetsBytes = 0;
    MPI_Allreduce(&offsetsLen, &totalOffsetsBytes, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    MPI_Offset typesHeaderStart = offsetsStart + totalOffsetsBytes;
    
    MPI_Offset typesHeaderLen = 0;
    if (rank == 0) {
        std::string typesHeader = "        </DataArray>\n";
        typesHeader += "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        
        typesHeaderLen = typesHeader.size();
        MPI_File_write_at(fh, typesHeaderStart, typesHeader.c_str(), 
                         typesHeaderLen, MPI_CHAR, MPI_STATUS_IGNORE);
    }
    
    MPI_Bcast(&typesHeaderLen, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
    
    // Write cell types (13 = VTK_WEDGE = triangular prism)
    std::ostringstream typesBuffer;
    typesBuffer.imbue(std::locale::classic());  // Safety: Ensure consistent formatting
    for (size_t i = 0; i < localPrisms.size(); i++) {
        typesBuffer << "          13\n";
    }
    std::string localTypesBuffer = typesBuffer.str();
    
    MPI_Offset typesLen = localTypesBuffer.size();
    MPI_Offset typesOffset = 0;
    MPI_Exscan(&typesLen, &typesOffset, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) typesOffset = 0;
    
    MPI_Offset typesStart = typesHeaderStart + typesHeaderLen;
    MPI_File_write_at_all(fh, typesStart + typesOffset, localTypesBuffer.c_str(), 
                          typesLen, MPI_CHAR, MPI_STATUS_IGNORE);
    
    // ===== STEP 7: WRITE CELL DATA HEADER (Rank 0 only) =====
    MPI_Offset totalTypesBytes = 0;
    MPI_Allreduce(&typesLen, &totalTypesBytes, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    MPI_Offset cellDataStart = typesStart + totalTypesBytes;
    
    MPI_Offset cellDataHeaderLen = 0;
    if (rank == 0) {
        std::string cellDataHeader = "        </DataArray>\n";
        cellDataHeader += "      </Cells>\n";
        cellDataHeader += "      <CellData>\n";
        cellDataHeader += "        <DataArray type=\"Float64\" Name=\"depth\" format=\"ascii\">\n";
        
        cellDataHeaderLen = cellDataHeader.size();
        MPI_File_write_at(fh, cellDataStart, cellDataHeader.c_str(), cellDataHeaderLen, 
                         MPI_CHAR, MPI_STATUS_IGNORE);
    }
    
    // Broadcast cell data header length
    MPI_Bcast(&cellDataHeaderLen, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
    
    // ===== STEP 8: WRITE DEPTH DATA (Parallel) =====
    std::ostringstream depthBuffer;
    depthBuffer.imbue(std::locale::classic());
    depthBuffer << std::fixed << std::setprecision(6);
    for (const auto& p : localPrisms) {
        depthBuffer << "          " << p.depth << "\n";
    }
    std::string localDepthBuffer = depthBuffer.str();
    
    MPI_Offset depthLen = localDepthBuffer.size();
    MPI_Offset depthOffset = 0;
    MPI_Exscan(&depthLen, &depthOffset, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) depthOffset = 0;
    
    MPI_Offset depthDataStart = cellDataStart + cellDataHeaderLen;
    MPI_File_write_at_all(fh, depthDataStart + depthOffset, localDepthBuffer.c_str(), 
                          depthLen, MPI_CHAR, MPI_STATUS_IGNORE);
    
    // ===== STEP 9: WRITE LAYER_ID HEADER AND DATA (Rank 0 header + Parallel data) =====
    MPI_Offset totalDepthBytes = 0;
    MPI_Allreduce(&depthLen, &totalDepthBytes, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    MPI_Offset layerHeaderStart = depthDataStart + totalDepthBytes;
    
    MPI_Offset layerHeaderLen = 0;
    if (rank == 0) {
        std::string layerHeader = "        </DataArray>\n";
        layerHeader += "        <DataArray type=\"Int32\" Name=\"layer_id\" format=\"ascii\">\n";
        
        layerHeaderLen = layerHeader.size();
        MPI_File_write_at(fh, layerHeaderStart, layerHeader.c_str(), layerHeaderLen, 
                         MPI_CHAR, MPI_STATUS_IGNORE);
    }
    
    // Broadcast layer header length
    MPI_Bcast(&layerHeaderLen, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
    
    // Write layer_id data
    std::ostringstream layerBuffer;
    layerBuffer.imbue(std::locale::classic());
    for (const auto& p : localPrisms) {
        layerBuffer << "          " << p.layer_id << "\n";
    }
    std::string localLayerBuffer = layerBuffer.str();
    
    MPI_Offset layerLen = localLayerBuffer.size();
    MPI_Offset layerOffset = 0;
    MPI_Exscan(&layerLen, &layerOffset, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) layerOffset = 0;
    
    MPI_Offset layerDataStart = layerHeaderStart + layerHeaderLen;
    MPI_File_write_at_all(fh, layerDataStart + layerOffset, localLayerBuffer.c_str(), 
                          layerLen, MPI_CHAR, MPI_STATUS_IGNORE);
    
    // ===== STEP 10: WRITE FOOTER (Rank 0 only) =====
    MPI_Offset totalLayerBytes = 0;
    MPI_Allreduce(&layerLen, &totalLayerBytes, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
    MPI_Offset footerStart = layerDataStart + totalLayerBytes;
    
    if (rank == 0) {
        std::string footer = "        </DataArray>\n";
        footer += "      </CellData>\n";
        footer += "    </Piece>\n";
        footer += "  </UnstructuredGrid>\n";
        footer += "</VTKFile>\n";
        
        MPI_File_write_at(fh, footerStart, footer.c_str(), footer.size(), 
                         MPI_CHAR, MPI_STATUS_IGNORE);
        
        std::cout << "VTU file writing completed successfully!" << std::endl;
    }
    
    // Close file
    MPI_File_close(&fh);
    
    return true;
}
