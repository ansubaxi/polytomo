#ifndef VISUALIZERX_H
#define VISUALIZERX_H

#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkColorTransferFunction.h>
#include <vtkLookupTable.h>

class VisualizerX {
public:
    VisualizerX();
    ~VisualizerX();

    void SetZExaggeration(double factor);
    void VisualizeSurface(vtkSmartPointer<vtkPolyData> surfaceData);
    void SetColorMap(const std::string& colorMapName);
    
private:
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor> interactor;
    vtkSmartPointer<vtkPolyDataMapper> mapper;
    vtkSmartPointer<vtkActor> actor;
    vtkSmartPointer<vtkColorTransferFunction> colorTransferFunction;
    vtkSmartPointer<vtkLookupTable> lookupTable;
};

#endif // VISUALIZERX_H
