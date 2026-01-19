#ifndef VISUALIZERX_H
#define VISUALIZERX_H

#include <vector>
#include <string>
#include <iostream>

/**
 * Class for handling visualization tasks including surface visualization,
 * Z-exaggeration transform, and handling colormaps.
 */
class VisualizerX {
public:
    // New Member Variables
    float zExaggeration;
    std::vector<std::string> colorMap;

    // Constructor
    VisualizerX();

    // Method declarations
    void setZExaggeration(float exaggeration);
    void visualizeSurface(const std::vector<float>& data);
    void applyColormap(const std::string& mapName);

private:
    void loadDefaultColormap();
};

#endif // VISUALIZERX_H