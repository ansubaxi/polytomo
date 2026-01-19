#include "visualizerx.h"

#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QSlider>
#include <QCheckBox>
#include <QComboBox>
#include <QPushButton>
#include <QFileDialog>
#include <QLabel>
#include <QLineEdit>
#include <QGroupBox>
#include <QMessageBox>
#include <QTimer>
#include <QDebug>
#include <QApplication>
#include <QVTKOpenGLNativeWidget.h>

// VTK Includes
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetMapper.h>
#include <vtkProperty.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkCellData.h>
#include <vtkLightKit.h>
#include <vtkCamera.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkTextProperty.h>
#include <vtkTextActor.h>
#include <vtkAxesActor.h>
#include <vtkOutputWindow.h>
#include <vtkCaptionActor2D.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>

// NEW: Default constructor for standalone use
VisualizerX::VisualizerX(QWidget* parent)
    : QMainWindow(parent),
      vtuFile_(""),
      dataLoaded_(false),
      zExaggeration_(1.0),
      downsampleRate_(1),
      surfaceOnly_(false)
{
    vtkOutputWindow::SetGlobalWarningDisplay(0);
    setWindowTitle("PolyTomo Visualizer");
    resize(1600, 1000);
    setupUI();
    setupVTK();
}

// Keep for backward compatibility
VisualizerX::VisualizerX(const QString& vtuFile, QWidget* parent)
    : QMainWindow(parent),
      vtuFile_(vtuFile),
      dataLoaded_(false),
      zExaggeration_(1.0),
      downsampleRate_(1),
      surfaceOnly_(false)
{
    vtkOutputWindow::SetGlobalWarningDisplay(0);
    setWindowTitle("PolyTomo Visualizer - " + QFileInfo(vtuFile).fileName());
    resize(1600, 1000);
    setupUI();
    setupVTK();
}

VisualizerX::~VisualizerX() {
    if (axesWidget_) {
        axesWidget_->SetInteractor(nullptr);
        axesWidget_->EnabledOff();
    }
}

void VisualizerX::showEvent(QShowEvent* event) {
    QMainWindow::showEvent(event);
    if (!dataLoaded_ && !vtuFile_.isEmpty()) {
        QTimer::singleShot(100, this, &VisualizerX::loadData);
    }
}

void VisualizerX::setupUI() {
    QWidget* central = new QWidget(this);
    setCentralWidget(central);
    QHBoxLayout* mainLayout = new QHBoxLayout(central);
    mainLayout->setContentsMargins(0,0,0,0);

    // --- Detect theme from application palette ---
    bool isDark = (qApp->palette().color(QPalette::Window).lightness() < 128);
    
    // --- Sidebar Controls ---
    QWidget* sidebar = new QWidget();
    sidebar->setFixedWidth(300);
    
    // Dynamic theme based on main window
    QString sidebarBg = isDark ? "#2b2b2b" : "#f5f5f5";
    QString groupBorder = isDark ? "#555" : "#ccc";
    QString textColor = isDark ? "#ffffff" : "#000000";
    QString labelColor = isDark ? "#cccccc" : "#333333";
    QString buttonBg = isDark ? "#3d3d3d" : "#e0e0e0";
    QString buttonHover = isDark ? "#4d4d4d" : "#d0d0d0";
    QString inputBg = isDark ? "#3d3d3d" : "#ffffff";
    QString inputBorder = isDark ? "#555" : "#ccc";
    
    sidebar->setStyleSheet(QString(
        "QWidget { background-color: %1; }"
        "QGroupBox { "
        "    color: %2; "
        "    border: 1px solid %3; "
        "    border-radius: 4px; "
        "    margin-top: 8px; "
        "    padding-top: 8px; "
        "    font-weight: bold; "
        "}"
        "QGroupBox::title { "
        "    subcontrol-origin: margin; "
        "    subcontrol-position: top left; "
        "    left: 10px; "
        "    padding: 0 5px; "
        "}"
        "QLabel { color: %4; }"
        "QPushButton { "
        "    background-color: %5; "
        "    color: %2; "
        "    border: none; "
        "    padding: 6px 12px; "
        "    border-radius: 3px; "
        "}"
        "QPushButton:hover { background-color: %6; }"
        "QLineEdit, QComboBox { "
        "    background-color: %7; "
        "    color: %2; "
        "    border: 1px solid %8; "
        "    padding: 4px; "
        "    border-radius: 3px; "
        "}"
    ).arg(sidebarBg, textColor, groupBorder, labelColor, 
          buttonBg, buttonHover, inputBg, inputBorder));

    QVBoxLayout* sideLayout = new QVBoxLayout(sidebar);
    sideLayout->setSpacing(5);
    sideLayout->setContentsMargins(10,10,10,10);
    
    QLabel* title = new QLabel("PolyTomo Visualizer");
    title->setStyleSheet("font-size: 16px; font-weight: bold; padding: 10px 0;");
    title->setAlignment(Qt::AlignCenter);
    sideLayout->addWidget(title);
    sideLayout->addSpacing(10);  // Add spacing after title
    
    // File Loading Group
    QGroupBox* grpFile = new QGroupBox("Mesh File");
    QVBoxLayout* lFile = new QVBoxLayout(grpFile);
    
    filePathEdit_ = new QLineEdit();
    filePathEdit_->setReadOnly(true);
    filePathEdit_->setPlaceholderText("No file loaded");
    if (!vtuFile_.isEmpty()) {
        filePathEdit_->setText(QFileInfo(vtuFile_).fileName());
    }
    lFile->addWidget(filePathEdit_);
    
    QHBoxLayout* lFileBtns = new QHBoxLayout();
    browseBtn_ = new QPushButton("Browse...");
    connect(browseBtn_, &QPushButton::clicked, this, &VisualizerX::onBrowseFile);
    lFileBtns->addWidget(browseBtn_);
    
    reloadBtn_ = new QPushButton("Reload");
    reloadBtn_->setEnabled(!vtuFile_.isEmpty());
    connect(reloadBtn_, &QPushButton::clicked, this, &VisualizerX::loadData);
    lFileBtns->addWidget(reloadBtn_);
    
    lFile->addLayout(lFileBtns);
    sideLayout->addWidget(grpFile);
    sideLayout->addSpacing(15);  // Add spacing after group

    // 1. Geometry Group
    QGroupBox* grpGeo = new QGroupBox("Geometry");
    QVBoxLayout* lGeo = new QVBoxLayout(grpGeo);
    
    lGeo->addWidget(new QLabel("Z Exaggeration:"));
    QHBoxLayout* lZ = new QHBoxLayout();
    zExaggerationSlider_ = new QSlider(Qt::Horizontal);
    zExaggerationSlider_->setRange(1, 50);
    zExaggerationSlider_->setValue(10);
    zExaggerationSlider_->setTickPosition(QSlider::TicksBelow);
    zExaggerationSlider_->setTickInterval(5);
    connect(zExaggerationSlider_, &QSlider::valueChanged, this, &VisualizerX::onZExaggerationChanged);
    zExagLabel_ = new QLabel("1.0x");
    zExagLabel_->setStyleSheet("color: #4CAF50; font-weight: bold; min-width: 50px;");
    lZ->addWidget(zExaggerationSlider_);
    lZ->addWidget(zExagLabel_);
    lGeo->addLayout(lZ);
    
    // NEW: Downsampling slider
    lGeo->addWidget(new QLabel("Downsample (every Nth):"));
    QHBoxLayout* lDownsample = new QHBoxLayout();
    downsampleSlider_ = new QSlider(Qt::Horizontal);
    downsampleSlider_->setRange(1, 20);
    downsampleSlider_->setValue(1);
    downsampleSlider_->setTickPosition(QSlider::TicksBelow);
    downsampleSlider_->setTickInterval(2);
    connect(downsampleSlider_, &QSlider::valueChanged, this, &VisualizerX::onDownsampleChanged);
    downsampleLabel_ = new QLabel("1 (all)");
    downsampleLabel_->setStyleSheet("color: #2196F3; font-weight: bold; min-width: 80px;");
    lDownsample->addWidget(downsampleSlider_);
    lDownsample->addWidget(downsampleLabel_);
    lGeo->addLayout(lDownsample);
    
    // NEW: Cell count label
    cellCountLabel_ = new QLabel("Cells: 0");
    cellCountLabel_->setStyleSheet("color: #FF9800; font-size: 11px;");
    lGeo->addWidget(cellCountLabel_);
    
    sideLayout->addWidget(grpGeo);
    sideLayout->addSpacing(15);  // Add spacing after group

    // 2. Display Group
    QGroupBox* grpDisp = new QGroupBox("Display Options");
    QVBoxLayout* lDisp = new QVBoxLayout(grpDisp);
    
    // NEW: Surface-only checkbox
    surfaceOnlyCheck_ = new QCheckBox("Surface Prisms Only");
    surfaceOnlyCheck_->setToolTip("Show only top layer to verify topography");
    connect(surfaceOnlyCheck_, &QCheckBox::toggled, this, &VisualizerX::onSurfaceOnlyToggled);
    lDisp->addWidget(surfaceOnlyCheck_);
    
    wireframeCheck_ = new QCheckBox("Wireframe Mode");
    connect(wireframeCheck_, &QCheckBox::toggled, this, &VisualizerX::onWireframeToggled);
    lDisp->addWidget(wireframeCheck_);

    edgesCheck_ = new QCheckBox("Show Cell Edges");
    connect(edgesCheck_, &QCheckBox::toggled, this, &VisualizerX::onEdgesToggled);
    lDisp->addWidget(edgesCheck_);

    cubeAxesCheck_ = new QCheckBox("Show Axes Ruler");
    cubeAxesCheck_->setChecked(true);
    connect(cubeAxesCheck_, &QCheckBox::toggled, this, &VisualizerX::onCubeAxesToggled);
    lDisp->addWidget(cubeAxesCheck_);
    sideLayout->addWidget(grpDisp);
    sideLayout->addSpacing(15);  // Add spacing after group

    // 3. Coloring Group
    QGroupBox* grpColor = new QGroupBox("Color Scheme");
    QVBoxLayout* lColor = new QVBoxLayout(grpColor);
    
    colorSchemeCombo_ = new QComboBox();
    colorSchemeCombo_->addItems({
        "Depth (Viridis)", 
        "Depth (Plasma)", 
        "Depth (Blue→Red)", 
        "Layer Index", 
        "Uniform Gray"
    });
    // No inline style - will use sidebar style
    connect(colorSchemeCombo_, QOverload<int>::of(&QComboBox::currentIndexChanged), 
            this, &VisualizerX::onColorSchemeChanged);
    lColor->addWidget(colorSchemeCombo_);
    sideLayout->addWidget(grpColor);
    sideLayout->addSpacing(20);  // Extra spacing before actions

    // 4. Actions
    QPushButton* btnReset = new QPushButton("Reset Camera");
    connect(btnReset, &QPushButton::clicked, this, &VisualizerX::onResetCamera);
    sideLayout->addWidget(btnReset);

    screenshotBtn_ = new QPushButton("Take Screenshot (3x)");
    connect(screenshotBtn_, &QPushButton::clicked, this, &VisualizerX::onSaveScreenshot);
    sideLayout->addWidget(screenshotBtn_);

    sideLayout->addStretch();
    
    statusLabel_ = new QLabel("Ready. Use Browse to load a mesh file.");
    if (!vtuFile_.isEmpty()) {
        statusLabel_->setText("Loading...");
    }
    statusLabel_->setStyleSheet("color: #888; font-size: 11px;");
    statusLabel_->setWordWrap(true);
    sideLayout->addWidget(statusLabel_);

    mainLayout->addWidget(sidebar);
    
    // --- VTK Widget ---
    vtkNativeWidget_ = new QVTKOpenGLNativeWidget();
    mainLayout->addWidget(vtkNativeWidget_, 1);
}

void VisualizerX::setupVTK() {
    renderWindow_ = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    vtkNativeWidget_->setRenderWindow(renderWindow_);

    renderer_ = vtkSmartPointer<vtkRenderer>::New();
    renderer_->SetBackground(0.15, 0.15, 0.18);  // Dark blue-gray
    renderWindow_->AddRenderer(renderer_);

    // Lighting
    vtkSmartPointer<vtkLightKit> lightKit = vtkSmartPointer<vtkLightKit>::New();
    lightKit->AddLightsToRenderer(renderer_);

    // Mapper and Actor
    mapper_ = vtkSmartPointer<vtkDataSetMapper>::New();
    meshActor_ = vtkSmartPointer<vtkActor>::New();
    meshActor_->SetMapper(mapper_);
    renderer_->AddActor(meshActor_);

    // Scalar Bar
    scalarBar_ = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar_->SetTitle("Depth (m)");
    scalarBar_->SetNumberOfLabels(5);
    scalarBar_->GetLabelTextProperty()->SetColor(1, 1, 1);
    scalarBar_->GetTitleTextProperty()->SetColor(1, 1, 1);
    scalarBar_->GetLabelTextProperty()->SetFontSize(12);
    scalarBar_->GetTitleTextProperty()->SetFontSize(14);
    scalarBar_->SetPosition(0.88, 0.1);
    scalarBar_->SetWidth(0.1);
    scalarBar_->SetHeight(0.8);
    renderer_->AddViewProp(scalarBar_);

    // Cube Axes
    cubeAxes_ = vtkSmartPointer<vtkCubeAxesActor>::New();
    cubeAxes_->SetCamera(renderer_->GetActiveCamera());
    cubeAxes_->GetTitleTextProperty(0)->SetColor(1.0, 1.0, 1.0);
    cubeAxes_->GetLabelTextProperty(0)->SetColor(1.0, 1.0, 1.0);
    cubeAxes_->GetTitleTextProperty(1)->SetColor(1.0, 1.0, 1.0);
    cubeAxes_->GetLabelTextProperty(1)->SetColor(1.0, 1.0, 1.0);
    cubeAxes_->GetTitleTextProperty(2)->SetColor(1.0, 1.0, 1.0);
    cubeAxes_->GetLabelTextProperty(2)->SetColor(1.0, 1.0, 1.0);
    cubeAxes_->DrawXGridlinesOn();
    cubeAxes_->DrawYGridlinesOn();
    cubeAxes_->DrawZGridlinesOn();
    cubeAxes_->SetXTitle("Easting (m)");
    cubeAxes_->SetYTitle("Northing (m)");
    cubeAxes_->SetZTitle("Elevation (m)");
    cubeAxes_->SetFlyModeToOuterEdges();
    renderer_->AddActor(cubeAxes_);

    // Orientation Axes
    vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
    axes->SetTotalLength(1.5, 1.5, 1.5);
    axes->SetShaftTypeToLine();
    axes->GetXAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
    axes->GetYAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
    axes->GetZAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();

    axesWidget_ = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    axesWidget_->SetOrientationMarker(axes);
    axesWidget_->SetInteractor(renderWindow_->GetInteractor());
    axesWidget_->SetViewport(0.0, 0.0, 0.2, 0.2);
    axesWidget_->SetEnabled(1);
    axesWidget_->InteractiveOff();
}

void VisualizerX::onBrowseFile() {
    QString file = QFileDialog::getOpenFileName(
        this,
        "Open VTU Mesh File",
        QString(),
        "VTU Files (*.vtu);;All Files (*)"
    );
    
    if (!file.isEmpty()) {
        vtuFile_ = file;
        filePathEdit_->setText(QFileInfo(file).fileName());
        setWindowTitle("PolyTomo Visualizer - " + QFileInfo(file).fileName());
        reloadBtn_->setEnabled(true);
        loadData();
    }
}

void VisualizerX::loadData() {
    if (vtuFile_.isEmpty()) {
        QMessageBox::warning(this, "No File", "Please select a mesh file first.");
        return;
    }

    loadMeshFile(vtuFile_);
}

void VisualizerX::loadMeshFile(const QString& filepath) {
    if (!QFile::exists(filepath)) {
        QMessageBox::critical(this, "Error", "File not found: " + filepath);
        statusLabel_->setText("ERROR: File not found");
        return;
    }

    statusLabel_->setText("Loading mesh...");
    QApplication::processEvents();

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = 
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filepath.toStdString().c_str());
    reader->Update();

    grid_ = reader->GetOutput();
    
    if (!grid_ || grid_->GetNumberOfCells() == 0) {
        QMessageBox::critical(this, "Error", "Failed to load mesh or mesh is empty");
        statusLabel_->setText("ERROR: Failed to load mesh");
        return;
    }

    dataLoaded_ = true;
    
    // Apply filters (downsampling/surface extraction)
    applyFilters();
    
    // Initial color scheme
    applyViridisColor();
    
    // Update cube axes bounds
    double bounds[6];
    filteredGrid_->GetBounds(bounds);
    cubeAxes_->SetBounds(bounds);
    
    // Reset camera
    renderer_->ResetCamera();
    renderWindow_->Render();
    
    statusLabel_->setText("Mesh loaded successfully");
}

// NEW: Apply downsampling and surface extraction filters
void VisualizerX::applyFilters() {
    if (!grid_ || grid_->GetNumberOfCells() == 0) {
        return;
    }
    
    vtkSmartPointer<vtkUnstructuredGrid> workingGrid = grid_;
    
    // STEP 1: Surface extraction (if enabled)
    if (surfaceOnly_) {
        // Extract cells where layer_id == 0 (top layer)
        vtkDataArray* layerArray = workingGrid->GetCellData()->GetArray("layer_id");
        
        if (layerArray) {
            vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
            threshold->SetInputData(workingGrid);
            threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "layer_id");
            threshold->SetLowerThreshold(0);
            threshold->SetUpperThreshold(0);
            threshold->Update();
            
            workingGrid = threshold->GetOutput();
            statusLabel_->setText("Surface layer extracted");
        } else {
            statusLabel_->setText("Warning: No layer_id array found");
        }
    }
    
    // STEP 2: Downsampling (if rate > 1)
    if (downsampleRate_ > 1) {
        vtkIdType numCells = workingGrid->GetNumberOfCells();
        vtkSmartPointer<vtkUnstructuredGrid> downsampledGrid = 
            vtkSmartPointer<vtkUnstructuredGrid>::New();
        
        vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
        downsampledGrid->SetPoints(newPoints);
        
        // Copy cell and point data arrays structure
        downsampledGrid->GetCellData()->CopyAllocate(workingGrid->GetCellData());
        
        // Sample every Nth cell
        vtkIdType sampledCount = 0;
        for (vtkIdType i = 0; i < numCells; i += downsampleRate_) {
            vtkCell* cell = workingGrid->GetCell(i);
            vtkIdList* pointIds = cell->GetPointIds();
            
            vtkSmartPointer<vtkIdList> newPointIds = vtkSmartPointer<vtkIdList>::New();
            
            for (vtkIdType j = 0; j < pointIds->GetNumberOfIds(); ++j) {
                vtkIdType oldPtId = pointIds->GetId(j);
                double pt[3];
                workingGrid->GetPoint(oldPtId, pt);
                vtkIdType newPtId = newPoints->InsertNextPoint(pt);
                newPointIds->InsertNextId(newPtId);
            }
            
            downsampledGrid->InsertNextCell(cell->GetCellType(), newPointIds);
            downsampledGrid->GetCellData()->CopyData(workingGrid->GetCellData(), i, sampledCount);
            sampledCount++;
        }
        
        workingGrid = downsampledGrid;
    }
    
    filteredGrid_ = workingGrid;
    
    // Update mapper
    mapper_->SetInputData(filteredGrid_);
    mapper_->Update();
    
    // Update cell count label
    updateCellCountLabel();
}

void VisualizerX::updateCellCountLabel() {
    if (!filteredGrid_) {
        cellCountLabel_->setText("Cells: 0");
        return;
    }
    
    vtkIdType displayedCells = filteredGrid_->GetNumberOfCells();
    vtkIdType totalCells = grid_ ? grid_->GetNumberOfCells() : 0;
    
    QString text;
    if (displayedCells < totalCells) {
        text = QString("Cells: %1 / %2 (%3%)")
            .arg(displayedCells)
            .arg(totalCells)
            .arg(int(100.0 * displayedCells / totalCells));
    } else {
        text = QString("Cells: %1").arg(displayedCells);
    }
    
    cellCountLabel_->setText(text);
}

void VisualizerX::onZExaggerationChanged(int value) {
    zExaggeration_ = value / 10.0;
    zExagLabel_->setText(QString("%1x").arg(zExaggeration_, 0, 'f', 1));
    
    if (!meshActor_) return;
    
    meshActor_->SetScale(1.0, 1.0, zExaggeration_);
    renderWindow_->Render();
}

// NEW: Downsampling slider handler
void VisualizerX::onDownsampleChanged(int value) {
    downsampleRate_ = value;
    
    if (value == 1) {
        downsampleLabel_->setText("1 (all)");
    } else {
        downsampleLabel_->setText(QString("%1").arg(value));
    }
    
    if (dataLoaded_) {
        applyFilters();
        onColorSchemeChanged(colorSchemeCombo_->currentIndex());
        renderWindow_->Render();
    }
}

// NEW: Surface-only checkbox handler
void VisualizerX::onSurfaceOnlyToggled(bool checked) {
    surfaceOnly_ = checked;
    
    if (dataLoaded_) {
        applyFilters();
        onColorSchemeChanged(colorSchemeCombo_->currentIndex());
        renderWindow_->Render();
    }
}

void VisualizerX::onWireframeToggled(bool checked) {
    if (!meshActor_) return;
    meshActor_->GetProperty()->SetRepresentationToSurface();
    if (checked) {
        meshActor_->GetProperty()->SetRepresentationToWireframe();
    }
    renderWindow_->Render();
}

void VisualizerX::onEdgesToggled(bool checked) {
    if (!meshActor_) return;
    meshActor_->GetProperty()->SetEdgeVisibility(checked);
    meshActor_->GetProperty()->SetEdgeColor(0.2, 0.2, 0.2);
    renderWindow_->Render();
}

void VisualizerX::onCubeAxesToggled(bool checked) {
    if (!cubeAxes_) return;
    cubeAxes_->SetVisibility(checked);
    renderWindow_->Render();
}

void VisualizerX::onColorSchemeChanged(int index) {
    if (!dataLoaded_) return;
    
    switch (index) {
        case 0: applyViridisColor(); break;
        case 1: applyPlasmaColor(); break;
        case 2: applyBlueRedColor(); break;
        case 3: applyLayerColor(); break;
        case 4: applySolidColor(); break;
    }
    
    renderWindow_->Render();
}

void VisualizerX::applyViridisColor() {
    if (!filteredGrid_ || !filteredGrid_->GetCellData()->HasArray("depth")) {
        applySolidColor();
        return;
    }
    
    vtkDataArray* depthArray = filteredGrid_->GetCellData()->GetArray("depth");
    double range[2];
    depthArray->GetRange(range);
    
    lut_ = vtkSmartPointer<vtkLookupTable>::New();
    lut_->SetNumberOfTableValues(256);
    lut_->SetHueRange(0.667, 0.0);  // Blue to Red
    lut_->SetSaturationRange(1.0, 1.0);
    lut_->SetValueRange(0.5, 1.0);
    lut_->Build();
    
    // Viridis-like colors
    for (int i = 0; i < 256; ++i) {
        double t = i / 255.0;
        double r, g, b;
        
        // Viridis color map approximation
        if (t < 0.25) {
            double u = t / 0.25;
            r = 0.267 * (1-u) + 0.282 * u;
            g = 0.005 * (1-u) + 0.141 * u;
            b = 0.329 * (1-u) + 0.490 * u;
        } else if (t < 0.5) {
            double u = (t - 0.25) / 0.25;
            r = 0.282 * (1-u) + 0.255 * u;
            g = 0.141 * (1-u) + 0.384 * u;
            b = 0.490 * (1-u) + 0.553 * u;
        } else if (t < 0.75) {
            double u = (t - 0.5) / 0.25;
            r = 0.255 * (1-u) + 0.478 * u;
            g = 0.384 * (1-u) + 0.647 * u;
            b = 0.553 * (1-u) + 0.408 * u;
        } else {
            double u = (t - 0.75) / 0.25;
            r = 0.478 * (1-u) + 0.993 * u;
            g = 0.647 * (1-u) + 0.906 * u;
            b = 0.408 * (1-u) + 0.145 * u;
        }
        
        lut_->SetTableValue(i, r, g, b, 1.0);
    }
    
    mapper_->SetScalarModeToUseCellData();
    mapper_->SelectColorArray("depth");
    mapper_->SetLookupTable(lut_);
    mapper_->SetScalarRange(range);
    mapper_->ScalarVisibilityOn();
    
    scalarBar_->SetLookupTable(lut_);
    scalarBar_->SetTitle("Depth (m)");
    scalarBar_->SetVisibility(1);
}

void VisualizerX::applyPlasmaColor() {
    if (!filteredGrid_ || !filteredGrid_->GetCellData()->HasArray("depth")) {
        applySolidColor();
        return;
    }
    
    vtkDataArray* depthArray = filteredGrid_->GetCellData()->GetArray("depth");
    double range[2];
    depthArray->GetRange(range);
    
    lut_ = vtkSmartPointer<vtkLookupTable>::New();
    lut_->SetNumberOfTableValues(256);
    
    // Plasma colormap
    for (int i = 0; i < 256; ++i) {
        double t = i / 255.0;
        double r, g, b;
        
        if (t < 0.25) {
            double u = t / 0.25;
            r = 0.050 * (1-u) + 0.514 * u;
            g = 0.030 * (1-u) + 0.064 * u;
            b = 0.529 * (1-u) + 0.687 * u;
        } else if (t < 0.5) {
            double u = (t - 0.25) / 0.25;
            r = 0.514 * (1-u) + 0.788 * u;
            g = 0.064 * (1-u) + 0.109 * u;
            b = 0.687 * (1-u) + 0.530 * u;
        } else if (t < 0.75) {
            double u = (t - 0.5) / 0.25;
            r = 0.788 * (1-u) + 0.957 * u;
            g = 0.109 * (1-u) + 0.411 * u;
            b = 0.530 * (1-u) + 0.185 * u;
        } else {
            double u = (t - 0.75) / 0.25;
            r = 0.957 * (1-u) + 0.940 * u;
            g = 0.411 * (1-u) + 0.976 * u;
            b = 0.185 * (1-u) + 0.131 * u;
        }
        
        lut_->SetTableValue(i, r, g, b, 1.0);
    }
    
    lut_->Build();
    
    mapper_->SetScalarModeToUseCellData();
    mapper_->SelectColorArray("depth");
    mapper_->SetLookupTable(lut_);
    mapper_->SetScalarRange(range);
    mapper_->ScalarVisibilityOn();
    
    scalarBar_->SetLookupTable(lut_);
    scalarBar_->SetTitle("Depth (m)");
    scalarBar_->SetVisibility(1);
}

void VisualizerX::applyBlueRedColor() {
    if (!filteredGrid_ || !filteredGrid_->GetCellData()->HasArray("depth")) {
        applySolidColor();
        return;
    }
    
    vtkDataArray* depthArray = filteredGrid_->GetCellData()->GetArray("depth");
    double range[2];
    depthArray->GetRange(range);
    
    lut_ = vtkSmartPointer<vtkLookupTable>::New();
    lut_->SetNumberOfTableValues(256);
    lut_->SetHueRange(0.667, 0.0);  // Blue to Red
    lut_->Build();
    
    mapper_->SetScalarModeToUseCellData();
    mapper_->SelectColorArray("depth");
    mapper_->SetLookupTable(lut_);
    mapper_->SetScalarRange(range);
    mapper_->ScalarVisibilityOn();
    
    scalarBar_->SetLookupTable(lut_);
    scalarBar_->SetTitle("Depth (m)");
    scalarBar_->SetVisibility(1);
}

void VisualizerX::applyLayerColor() {
    if (!filteredGrid_ || !filteredGrid_->GetCellData()->HasArray("layer_id")) {
        applySolidColor();
        return;
    }
    
    vtkDataArray* layerArray = filteredGrid_->GetCellData()->GetArray("layer_id");
    double range[2];
    layerArray->GetRange(range);
    
    int numLayers = static_cast<int>(range[1] - range[0]) + 1;
    
    lut_ = vtkSmartPointer<vtkLookupTable>::New();
    lut_->SetNumberOfTableValues(numLayers);
    lut_->SetHueRange(0.0, 0.85);  // Full spectrum
    lut_->Build();
    
    mapper_->SetScalarModeToUseCellData();
    mapper_->SelectColorArray("layer_id");
    mapper_->SetLookupTable(lut_);
    mapper_->SetScalarRange(range);
    mapper_->ScalarVisibilityOn();
    
    scalarBar_->SetLookupTable(lut_);
    scalarBar_->SetTitle("Layer Index");
    scalarBar_->SetVisibility(1);
}

void VisualizerX::applySolidColor() {
    mapper_->ScalarVisibilityOff();
    meshActor_->GetProperty()->SetColor(0.7, 0.7, 0.7);
    scalarBar_->SetVisibility(0);
}

void VisualizerX::onSaveScreenshot() {
    if (!dataLoaded_) {
        QMessageBox::warning(this, "No Data", "Please load a mesh first");
        return;
    }
    
    QString filename = QFileDialog::getSaveFileName(
        this,
        "Save Screenshot",
        QFileInfo(vtuFile_).baseName() + "_screenshot.png",
        "PNG Images (*.png)"
    );
    
    if (filename.isEmpty()) return;
    
    vtkSmartPointer<vtkWindowToImageFilter> filter = 
        vtkSmartPointer<vtkWindowToImageFilter>::New();
    filter->SetInput(renderWindow_);
    filter->SetScale(3);  // 3x resolution
    filter->SetInputBufferTypeToRGBA();
    filter->ReadFrontBufferOff();
    filter->Update();
    
    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(filename.toStdString().c_str());
    writer->SetInputConnection(filter->GetOutputPort());
    writer->Write();
    
    statusLabel_->setText("Screenshot saved: " + QFileInfo(filename).fileName());
}

void VisualizerX::onResetCamera() {
    if (!renderer_) return;
    renderer_->ResetCamera();
    renderWindow_->Render();
}
