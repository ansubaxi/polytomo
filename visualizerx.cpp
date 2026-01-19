// Corrected visualizerx.cpp file

#include <iostream>
#include <vector>
#include <algorithm>

// Function to apply Z-exaggeration transform
void applyZExaggeration(std::vector<float>& heights, float exaggerationFactor) {
    for (auto& height : heights) {
        height *= exaggerationFactor;
    }
}

// Function to filter surface layers
std::vector<float> filterSurfaceLayers(const std::vector<float>& data, float threshold) {
    std::vector<float> filtered;
    std::copy_if(data.begin(), data.end(), std::back_inserter(filtered), [&](float value) {
        return value > threshold;
    });
    return filtered;
}

// Function to handle colormap data array
std::vector<int> generateColormap(size_t size) {
    std::vector<int> colormap(size);
    for (size_t i = 0; i < size; ++i) {
        colormap[i] = static_cast<int>(255.0f * (static_cast<float>(i) / size));
    }
    return colormap;
}

int main() {
    // Example usage
    std::vector<float> heights = {1.0, 2.0, 3.0};
    applyZExaggeration(heights, 2.0);
    
    std::vector<float> surfaceData = {1.5, 0.5, 3.5, 4.0};
    auto filteredData = filterSurfaceLayers(surfaceData, 2.0);
    
    auto colormap = generateColormap(256);
    
    // Print results
    std::cout << "Heights after Z-exaggeration:" << std::endl;
    for (const auto& height : heights) {
        std::cout << height << ' ';
    }
    std::cout << std::endl;
    
    std::cout << "Filtered surface data:" << std::endl;
    for (const auto& value : filteredData) {
        std::cout << value << ' ';
    }
    std::cout << std::endl;
    
    std::cout << "Colormap data:" << std::endl;
    for (const auto& color : colormap) {
        std::cout << color << ' ';
    }
    std::cout << std::endl;
    return 0;
}