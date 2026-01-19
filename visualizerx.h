#ifndef VISUALIZERX_H
#define VISUALIZERX_H

#include <QMainWindow>
#include <QFileInfo>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkActor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkCubeAxesActor.h>

class QSlider;
class QCheckBox;
class QComboBox;
class QPushButton;
class QLabel;
class QVTKOpenGLNativeWidget;
class QLineEdit;

class VisualizerX : public QMainWindow {
    Q_OBJECT

public:
    // NEW: Default constructor for standalone use
    explicit VisualizerX(QWidget* parent = nullptr);
    
    // Keep this for backward compatibility with mainwindow.cpp
    explicit VisualizerX(const QString& vtuFile, QWidget* parent = nullptr);
    
    ~VisualizerX();

protected:
    void showEvent(QShowEvent* event) override;

private slots:
    void onBrowseFile();  // NEW: Browse for mesh file
    void loadData();
    void onZExaggerationChanged(int value);
    void onDownsampleChanged(int value);      // NEW: Downsampling slider
    void onSurfaceOnlyToggled(bool checked);  // NEW: Surface extraction
    void onWireframeToggled(bool checked);
    void onEdgesToggled(bool checked);
    void onCubeAxesToggled(bool checked);
    void onColorSchemeChanged(int index);
    void onSaveScreenshot();
    void onResetCamera();

private:
    void setupUI();
    void setupVTK();
    void applyViridisColor();
    void applyPlasmaColor();
    void applyBlueRedColor();
    void applyLayerColor();
    void applySolidColor();
    void loadMeshFile(const QString& filepath);  // NEW: Actual loader
    void applyFilters();                         // NEW: Apply downsampling/surface extraction
    void updateCellCountLabel();                 // NEW: Update displayed cell count

    QString vtuFile_;
    bool dataLoaded_;
    double zExaggeration_;
    int downsampleRate_;       // NEW: Downsample every Nth cell
    bool surfaceOnly_;         // NEW: Show only surface layer
    
    // UI Elements
    QVTKOpenGLNativeWidget* vtkNativeWidget_;
    QLineEdit* filePathEdit_;       // NEW: Show current file
    QPushButton* browseBtn_;        // NEW: Browse button
    QPushButton* reloadBtn_;        // NEW: Reload button
    QSlider* zExaggerationSlider_;
    QSlider* downsampleSlider_;     // NEW: Downsampling slider
    QLabel* zExagLabel_;
    QLabel* downsampleLabel_;       // NEW: Downsample label
    QLabel* cellCountLabel_;        // NEW: Cell count display
    QCheckBox* surfaceOnlyCheck_;   // NEW: Surface-only checkbox
    QCheckBox* wireframeCheck_;
    QCheckBox* edgesCheck_;
    QCheckBox* cubeAxesCheck_;
    QComboBox* colorSchemeCombo_;
    QPushButton* screenshotBtn_;
    QLabel* statusLabel_;

    // VTK Objects
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;
    vtkSmartPointer<vtkRenderer> renderer_;
    vtkSmartPointer<vtkActor> meshActor_;
    vtkSmartPointer<vtkDataSetMapper> mapper_;
    vtkSmartPointer<vtkUnstructuredGrid> grid_;           // Original full grid
    vtkSmartPointer<vtkUnstructuredGrid> filteredGrid_;   // NEW: After filtering
    vtkSmartPointer<vtkLookupTable> lut_;
    vtkSmartPointer<vtkScalarBarActor> scalarBar_;
    vtkSmartPointer<vtkCubeAxesActor> cubeAxes_;
    vtkSmartPointer<vtkOrientationMarkerWidget> axesWidget_;
};

#endif // VISUALIZERX_H
