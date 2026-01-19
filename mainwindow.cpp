#include "mainwindow.h"
#include "modelgenerator.h"
#include "visualizerx.h"

// --- QT Headers ---
#include <QMenuBar>
#include <QMenu>
#include <QAction>
#include <QKeySequence>
#include <QApplication>
#include <QStyleFactory>
#include <QDateTime>
#include <QScrollBar>
#include <QFile>
#include <QTextStream>
#include <QRegExp>
#include <cmath>
#include <QHeaderView>
#include <QTimer>
#include <QScrollArea>
#include <QMessageBox>
#include <QFileInfo>
#include <QProcess>      // Added for MPI
#include <QProcessEnvironment> // For MPI environment
#include <QJsonDocument> // Added for MPI
#include <QJsonObject>   // Added for MPI
#include <QJsonArray>    // Added for MPI
#include <QDir>          // Added for MPI

// ==========================
// Constructor / Destructor
// ==========================
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent),
      generator(nullptr),
      visualizer(nullptr),
      generatorThread(nullptr),
      estimateUpdateTimer(nullptr),
      mpiProcess(nullptr), // Initialize MPI process
      totalDataPoints(0),
      isDarkMode(false)
{
    setWindowTitle("PolyTomo MeshGen v1.0 - Polyhedral Mesh Generator");
    setMinimumSize(1200, 800);

    // Initialize settings
    settings = new QSettings("PolyTomo", "MeshGen", this);

    // CRITICAL: Initialize Timer BEFORE setupUI() to prevent crashes
    estimateUpdateTimer = new QTimer(this);
    estimateUpdateTimer->setSingleShot(true);
    estimateUpdateTimer->setInterval(500); // 500ms delay for debouncing
    connect(estimateUpdateTimer, &QTimer::timeout, this, &MainWindow::updateMeshEstimate);

    // Build the UI
    setupUI();

    // Load settings (including theme preference)
    loadSettings();

    // Apply the saved theme
    applyCurrentTheme();

    // Initial log message
    logMessage("═══════════════════════════════════════════════════════", "HEADER");
    logMessage("PolyTomo MeshGen v1.0", "HEADER");
    logMessage("Polyhedral Potential Field Tomography", "HEADER");
    logMessage("Ready. Load topography file to begin.", "INFO");
}

MainWindow::~MainWindow() {
    saveSettings();
    if (generatorThread && generatorThread->isRunning()) {
        generatorThread->quit();
        generatorThread->wait();
    }
    // Clean up MPI process
    if (mpiProcess) {
        if (mpiProcess->state() == QProcess::Running) {
            mpiProcess->kill();
            mpiProcess->waitForFinished();
        }
        delete mpiProcess;
    }
    if (visualizer) {
        delete visualizer;
    }
}

// ==========================
// UI Setup
// ==========================
void MainWindow::setupUI() {
    createMenuBar();

    // Central widget with scroll area (Responsive Design)
    QScrollArea* scrollArea = new QScrollArea(this);
    scrollArea->setWidgetResizable(true);
    scrollArea->setFrameShape(QFrame::NoFrame);
    
    QWidget* scrollWidget = new QWidget();
    QVBoxLayout* mainLayout = new QVBoxLayout(scrollWidget);
    mainLayout->setSpacing(15);
    mainLayout->setContentsMargins(15, 15, 15, 15);

    // --- 1. Input Section ---
    createInputSection();
    mainLayout->addWidget(inputGroup);

    // --- 2. Grid Section ---
    createGridSection();
    mainLayout->addWidget(gridGroup);

    // --- 3. Depth Section ---
    createDepthSection();
    mainLayout->addWidget(depthGroup);

    // --- 4. Estimate Panel ---
    createEstimatePanel();
    mainLayout->addWidget(estimateGroup);

    // --- 5. Output Section ---
    createOutputSection();
    mainLayout->addWidget(outputGroup);

    // --- 6. MPI Section ---
    createMPISection();
    mainLayout->addWidget(mpiGroup);

    // --- 7. Actions ---
    createActionSection();
    mainLayout->addWidget(actionGroup);

    // Add stretch to push logs to the bottom
    mainLayout->addStretch();

    // --- 8. Log Section ---
    createLogSection();
    mainLayout->addWidget(logGroup);

    // Finalize Scroll Area
    scrollArea->setWidget(scrollWidget);
    setCentralWidget(scrollArea);
}

void MainWindow::createMenuBar() {
    QMenuBar* menuBar = new QMenuBar(this);
    setMenuBar(menuBar);

    QMenu* fileMenu = menuBar->addMenu("&File");
    QAction* openAction = fileMenu->addAction("&Open Topography...");
    openAction->setShortcut(QKeySequence::Open);
    connect(openAction, &QAction::triggered, this, &MainWindow::browseInputFile);

    fileMenu->addSeparator();
    QAction* exitAction = fileMenu->addAction("E&xit");
    exitAction->setShortcut(QKeySequence::Quit);
    connect(exitAction, &QAction::triggered, this, &QWidget::close);

    QMenu* viewMenu = menuBar->addMenu("&View");
    QAction* themeAction = viewMenu->addAction("Toggle &Dark/Light Mode");
    themeAction->setShortcut(QKeySequence("Ctrl+D"));
    connect(themeAction, &QAction::triggered, this, &MainWindow::toggleTheme);

    QMenu* helpMenu = menuBar->addMenu("&Help");
    helpMenu->addAction("&About", [this]() {
        QMessageBox::about(this, "About PolyTomo", 
            "PolyTomo MeshGen v1.0\nPolyhedral Potential Field Tomography");
    });
}

// ==========================
// Sections Creation
// ==========================
void MainWindow::createInputSection() {
    inputGroup = new QGroupBox("1. Input Topography Data");
    QGridLayout* layout = new QGridLayout();
    layout->setVerticalSpacing(12);
    
    // File Selection
    layout->addWidget(new QLabel("Topography File:"), 0, 0);
    inputFileEdit = new QLineEdit();
    inputFileEdit->setReadOnly(true);
    inputFileEdit->setPlaceholderText("Select CSV file with X, Y, Elevation data...");
    layout->addWidget(inputFileEdit, 0, 1);
    
    browseInputBtn = new QPushButton("Browse...");
    connect(browseInputBtn, &QPushButton::clicked, this, &MainWindow::browseInputFile);
    layout->addWidget(browseInputBtn, 0, 2);

    // Info Labels
    dataPointsLabel = new QLabel("Data points: Not loaded");
    layout->addWidget(dataPointsLabel, 1, 0, 1, 3);
    dataExtentLabel = new QLabel("Survey extent: Not calculated");
    layout->addWidget(dataExtentLabel, 2, 0, 1, 3);

    // Column Mapping
    layout->addWidget(new QLabel("Easting Column:"), 3, 0);
    eastingCombo = new QComboBox();
    layout->addWidget(eastingCombo, 3, 1, 1, 2);

    layout->addWidget(new QLabel("Northing Column:"), 4, 0);
    northingCombo = new QComboBox();
    layout->addWidget(northingCombo, 4, 1, 1, 2);

    layout->addWidget(new QLabel("Elevation Column:"), 5, 0);
    topoCombo = new QComboBox();
    layout->addWidget(topoCombo, 5, 1, 1, 2);

    // Fill default combos (block signals during initial setup)
    QStringList defaults = {"Column 0", "Column 1", "Column 2", "Column 3"};
    
    eastingCombo->blockSignals(true);
    northingCombo->blockSignals(true);
    topoCombo->blockSignals(true);
    
    eastingCombo->addItems(defaults); eastingCombo->setCurrentIndex(0);
    northingCombo->addItems(defaults); northingCombo->setCurrentIndex(1);
    topoCombo->addItems(defaults); topoCombo->setCurrentIndex(2);
    
    eastingCombo->blockSignals(false);
    northingCombo->blockSignals(false);
    topoCombo->blockSignals(false);
    
    // Connect column changes to rescan file
    connect(eastingCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &MainWindow::scanInputFile);
    connect(northingCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &MainWindow::scanInputFile);
    connect(topoCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &MainWindow::scanInputFile);

    inputGroup->setLayout(layout);
}

void MainWindow::createGridSection() {
    gridGroup = new QGroupBox("2. Prism Horizontal Dimensions");
    QGridLayout* layout = new QGridLayout();
    layout->setVerticalSpacing(12);

    layout->addWidget(new QLabel("Prism Size X (dx):"), 0, 0);
    dxSpinBox = new QDoubleSpinBox();
    dxSpinBox->setRange(0.1, 10000.0);
    dxSpinBox->setValue(10.0);
    dxSpinBox->setSuffix(" m");
    connect(dxSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), 
            this, &MainWindow::onPrismSizeChanged);
    layout->addWidget(dxSpinBox, 0, 1);

    layout->addWidget(new QLabel("Prism Size Y (dy):"), 1, 0);
    dySpinBox = new QDoubleSpinBox();
    dySpinBox->setRange(0.1, 10000.0);
    dySpinBox->setValue(10.0);
    dySpinBox->setSuffix(" m");
    connect(dySpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), 
            this, &MainWindow::onPrismSizeChanged);
    layout->addWidget(dySpinBox, 1, 1);

    layout->addWidget(new QLabel("Grid Offset:"), 2, 0);
    gridOffsetSpinBox = new QDoubleSpinBox();
    gridOffsetSpinBox->setRange(0.0, 1.0);
    gridOffsetSpinBox->setValue(0.1);
    gridOffsetSpinBox->setSingleStep(0.05);
    connect(gridOffsetSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), 
            this, &MainWindow::onPrismSizeChanged);
    layout->addWidget(gridOffsetSpinBox, 2, 1);

    prismInfoLabel = new QLabel("Prism base: 10.0 m × 10.0 m");
    prismInfoLabel->setStyleSheet("color: #0078d7; font-style: italic;");
    layout->addWidget(prismInfoLabel, 3, 0, 1, 2);

    gridGroup->setLayout(layout);
}

void MainWindow::createDepthSection() {
    depthGroup = new QGroupBox("3. Vertical Layer Configuration");
    QVBoxLayout* mainLayout = new QVBoxLayout();
    mainLayout->setSpacing(10);

    // Topography Aware Checkbox
    topoAwareCheckBox = new QCheckBox("Topography Aware (drape mesh under surface)");
    topoAwareCheckBox->setChecked(true);
    mainLayout->addWidget(topoAwareCheckBox);

    // Total Depth
    QHBoxLayout* depthLayout = new QHBoxLayout();
    depthLayout->addWidget(new QLabel("Total Inversion Depth:"));
    totalDepthSpinBox = new QDoubleSpinBox();
    totalDepthSpinBox->setRange(10.0, 50000.0);
    totalDepthSpinBox->setValue(1000.0);
    totalDepthSpinBox->setSuffix(" m");
    connect(totalDepthSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), 
            this, &MainWindow::updateTotalDepth);
    depthLayout->addWidget(totalDepthSpinBox);
    mainLayout->addLayout(depthLayout);

    // Layer Type Selector
    QHBoxLayout* typeLayout = new QHBoxLayout();
    typeLayout->addWidget(new QLabel("Layer Type:"));
    cellTypeCombo = new QComboBox();
    cellTypeCombo->addItem("Uniform Layers");
    cellTypeCombo->addItem("Variable Layers");
    cellTypeCombo->setCurrentIndex(1);
    connect(cellTypeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &MainWindow::onLayerTableChanged);
    typeLayout->addWidget(cellTypeCombo);
    mainLayout->addLayout(typeLayout);

    // Uniform Settings (Hidden by default based on combo)
    QHBoxLayout* uniformLayout = new QHBoxLayout();
    uniformLayout->addWidget(new QLabel("Uniform Thickness:"));
    uniformDzSpinBox = new QDoubleSpinBox();
    uniformDzSpinBox->setRange(1.0, 1000.0);
    uniformDzSpinBox->setValue(50.0);
    uniformDzSpinBox->setSuffix(" m");
    connect(uniformDzSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), 
            this, &MainWindow::updateTotalDepth);
    uniformLayout->addWidget(uniformDzSpinBox);
    mainLayout->addLayout(uniformLayout);

    // Variable Layer Table
    layerTable = new QTableWidget(0, 3);
    layerTable->setHorizontalHeaderLabels({"Thickness (m)", "Count", "Cumulative Depth"});
    layerTable->horizontalHeader()->setStretchLastSection(true);
    layerTable->setMinimumHeight(150);
    connect(layerTable, &QTableWidget::cellChanged, this, &MainWindow::onLayerTableChanged);
    mainLayout->addWidget(layerTable);

    // Buttons
    QHBoxLayout* btnLayout = new QHBoxLayout();
    addLayerBtn = new QPushButton("Add Layer");
    connect(addLayerBtn, &QPushButton::clicked, this, &MainWindow::addVariableLayer);
    btnLayout->addWidget(addLayerBtn);

    removeLayerBtn = new QPushButton("Remove Layer");
    connect(removeLayerBtn, &QPushButton::clicked, this, &MainWindow::removeVariableLayer);
    btnLayout->addWidget(removeLayerBtn);
    mainLayout->addLayout(btnLayout);

    // Summary Labels (Must be created before adding layers to prevent crash)
    totalLayersLabel = new QLabel("Total layers: 0");
    mainLayout->addWidget(totalLayersLabel);

    calculatedDepthLabel = new QLabel("Calculated depth: 0.0 m");
    mainLayout->addWidget(calculatedDepthLabel);

    // Init Logic
    addVariableLayer(); // Add default layer
    addVariableLayer();
    
    // Initial State update
    onLayerTableChanged(); // Hides/Shows widgets based on Combo

    depthGroup->setLayout(mainLayout);
}

void MainWindow::createEstimatePanel() {
    estimateGroup = new QGroupBox("4. Mesh Estimate");
    QGridLayout* layout = new QGridLayout();
    layout->setVerticalSpacing(8);

    layout->addWidget(new QLabel("Grid Dimensions:"), 0, 0);
    gridDimensionsLabel = new QLabel("Not calculated");
    gridDimensionsLabel->setStyleSheet("font-weight: bold;");
    layout->addWidget(gridDimensionsLabel, 0, 1);

    layout->addWidget(new QLabel("Total Prisms:"), 1, 0);
    totalCellsLabel = new QLabel("0");
    layout->addWidget(totalCellsLabel, 1, 1);

    layout->addWidget(new QLabel("Estimated Memory:"), 2, 0);
    memoryEstimateLabel = new QLabel("0.0 MB");
    layout->addWidget(memoryEstimateLabel, 2, 1);

    warningLabel = new QLabel("");
    warningLabel->setWordWrap(true);
    layout->addWidget(warningLabel, 3, 0, 1, 2);

    estimateGroup->setLayout(layout);
}

void MainWindow::createOutputSection() {
    outputGroup = new QGroupBox("5. Output Format");
    QVBoxLayout* layout = new QVBoxLayout();

    QHBoxLayout* fileLayout = new QHBoxLayout();
    fileLayout->addWidget(new QLabel("Output File:"));
    outputFileEdit = new QLineEdit();
    fileLayout->addWidget(outputFileEdit);
    browseOutputBtn = new QPushButton("Browse...");
    connect(browseOutputBtn, &QPushButton::clicked, this, &MainWindow::browseOutputFile);
    fileLayout->addWidget(browseOutputBtn);
    layout->addLayout(fileLayout);

    // Format Radio Buttons
    formatButtonGroup = new QButtonGroup(this);
    
    formatPMeshRadio = new QRadioButton("PolyTomo Format (.pmesh)");
    formatPMeshRadio->setChecked(true);
    formatButtonGroup->addButton(formatPMeshRadio);
    layout->addWidget(formatPMeshRadio);

    formatVTKRadio = new QRadioButton("VTK Format (.vtu)");
    formatButtonGroup->addButton(formatVTKRadio);
    layout->addWidget(formatVTKRadio);

    formatBothRadio = new QRadioButton("Both Formats");
    formatButtonGroup->addButton(formatBothRadio);
    layout->addWidget(formatBothRadio);

    outputGroup->setLayout(layout);
}

void MainWindow::createMPISection() {
    mpiGroup = new QGroupBox("6. Parallelization");
    QVBoxLayout* layout = new QVBoxLayout();

    useMPICheckbox = new QCheckBox("Enable MPI Parallelization");
    layout->addWidget(useMPICheckbox);

    QHBoxLayout* spinLayout = new QHBoxLayout();
    spinLayout->addWidget(new QLabel("Processes:"));
    mpiProcessesSpinBox = new QSpinBox();
    mpiProcessesSpinBox->setRange(1, 256);
    mpiProcessesSpinBox->setValue(4);
    mpiProcessesSpinBox->setEnabled(false); 
    connect(useMPICheckbox, &QCheckBox::toggled, mpiProcessesSpinBox, &QSpinBox::setEnabled);
    spinLayout->addWidget(mpiProcessesSpinBox);
    spinLayout->addStretch();
    layout->addLayout(spinLayout);

    mpiWarningLabel = new QLabel("Note: MPI requires mpiexec to be in system path.");
    mpiWarningLabel->setStyleSheet("color: orange; font-style: italic;");
    layout->addWidget(mpiWarningLabel);

    mpiGroup->setLayout(layout);
}

void MainWindow::createActionSection() {
    actionGroup = new QGroupBox("7. Generate");
    QVBoxLayout* layout = new QVBoxLayout();

    QHBoxLayout* btnLayout = new QHBoxLayout();
    generateBtn = new QPushButton("Generate Mesh (Serial)");
    generateBtn->setMinimumHeight(45);
    generateBtn->setStyleSheet("font-weight: bold; font-size: 14px;");
    connect(generateBtn, &QPushButton::clicked, this, &MainWindow::generateModel);
    
    // --- MPI Button Addition ---
    generateMPIBtn = new QPushButton("Generate Mesh (MPI Parallel)");
    generateMPIBtn->setMinimumHeight(45);
    generateMPIBtn->setStyleSheet(
        "QPushButton {"
        "  background-color: #28a745;"
        "  color: #ffffff !important;"
        "  font-weight: bold;"
        "  font-size: 14px;"
        "  border: 2px solid #28a745;"
        "  border-radius: 4px;"
        "  padding: 8px 16px;"
        "}"
        "QPushButton:hover {"
        "  background-color: #218838;"
        "  color: #ffffff !important;"
        "}"
        "QPushButton:disabled {"
        "  background-color: #6c757d;"
        "  color: #ffffff !important;"
        "}"
    );
    // Initially disabled if inputs invalid, updated in validateInputs
    generateMPIBtn->setEnabled(false);
    connect(generateMPIBtn, &QPushButton::clicked, this, &MainWindow::generateModelMPI);
    // ---------------------------

    visualizeBtn = new QPushButton("Visualize Mesh");
    visualizeBtn->setMinimumHeight(45);
    visualizeBtn->setEnabled(false);
    connect(visualizeBtn, &QPushButton::clicked, this, &MainWindow::visualizeModel);

    btnLayout->addWidget(generateBtn);
    btnLayout->addWidget(generateMPIBtn); // Add to layout
    btnLayout->addWidget(visualizeBtn);

    layout->addLayout(btnLayout);

    progressBar = new QProgressBar();
    progressBar->setRange(0, 100);
    progressBar->setValue(0);
    layout->addWidget(progressBar);

    actionGroup->setLayout(layout);
}

void MainWindow::createLogSection() {
    logGroup = new QGroupBox("Log");
    QVBoxLayout* layout = new QVBoxLayout();

    logTextEdit = new QTextEdit();
    logTextEdit->setReadOnly(true);
    logTextEdit->setMaximumHeight(150);
    // Use standard font for log
    logTextEdit->setStyleSheet("font-family: Consolas, Monospace; font-size: 11px;");
    layout->addWidget(logTextEdit);

    QHBoxLayout* btnLayout = new QHBoxLayout();
    clearLogBtn = new QPushButton("Clear Log");
    connect(clearLogBtn, &QPushButton::clicked, logTextEdit, &QTextEdit::clear);
    btnLayout->addWidget(clearLogBtn);

    saveLogBtn = new QPushButton("Save Log...");
    connect(saveLogBtn, &QPushButton::clicked, [this]() {
        QString filename = QFileDialog::getSaveFileName(this, "Save Log", "", "Text Files (*.txt)");
        if (!filename.isEmpty()) {
            QFile file(filename);
            if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
                QTextStream out(&file);
                out << logTextEdit->toPlainText();
                file.close();
                logMessage("Log saved to: " + filename, "SUCCESS");
            }
        }
    });
    btnLayout->addWidget(saveLogBtn);
    btnLayout->addStretch();
    layout->addLayout(btnLayout);

    logGroup->setLayout(layout);
}

// ==========================
// Logic & Slots
// ==========================

void MainWindow::browseInputFile() {
    QString filename = QFileDialog::getOpenFileName(
        this, "Select Input", "", "CSV Files (*.csv *.txt);;All Files (*)");

    if (!filename.isEmpty()) {
        inputFileEdit->setText(filename);
        currentInputFile = filename;

        // Reset UI immediately to indicate loading
        dataPointsLabel->setText("Loading...");
        dataExtentLabel->setText("Calculating extent...");
        inputFileValid = false;

        QFile file(filename);
        if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            QTextStream in(&file);
            QString headerLine = in.readLine();
            QStringList headers = headerLine.split(QRegExp("[,\\s]+"), Qt::SkipEmptyParts);

            if (!headers.isEmpty()) {
                bool isNum;
                headers[0].toDouble(&isNum);
                if (!isNum) {
                    eastingCombo->clear();
                    northingCombo->clear();
                    topoCombo->clear();
                    eastingCombo->addItems(headers);
                    northingCombo->addItems(headers);
                    topoCombo->addItems(headers);
                    
                    if (headers.size() > 2) {
                        eastingCombo->setCurrentIndex(0);
                        northingCombo->setCurrentIndex(1);
                        topoCombo->setCurrentIndex(2);
                    }
                }
            }
            file.close();
            onInputFileLoaded();
        }
        validateInputs();
    }
}

void MainWindow::onInputFileLoaded() {
    QFile file(currentInputFile);
    if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream in(&file);
        int count = 0;
        while (!in.atEnd()) {
            if (!in.readLine().trimmed().isEmpty()) count++;
        }
        if (count > 0) count--; 
        
        totalDataPoints = count;
        dataPointsLabel->setText(QString("Data points: %1").arg(count));
        logMessage(QString("Loaded input file with %1 points").arg(count), "INFO");
        inputFileValid = true;
        
        // Scan file for extents
        scanInputFile();
        
        if(estimateUpdateTimer) estimateUpdateTimer->start();
    }
}

void MainWindow::scanInputFile() {
    extentsValid = false;
    
    if (currentInputFile.isEmpty()) return;
    
    int eastCol = eastingCombo->currentIndex();
    int northCol = northingCombo->currentIndex();
    int topoCol = topoCombo->currentIndex();
    
    // Safety check for invalid indices (e.g., -1 when combo is empty)
    if (eastCol < 0 || northCol < 0 || topoCol < 0) {
        dataExtentLabel->setText("Select columns first");
        return;
    }
    
    QFile file(currentInputFile);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        logMessage("Failed to open input file for scanning", "ERROR");
        return;
    }
    
    QTextStream in(&file);
    
    std::vector<double> x, y, z;
    
    int lineNum = 0;
    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        lineNum++;
        
        if (line.isEmpty()) continue;
        
        QStringList tokens = line.split(QRegExp("[,\\s]+"), Qt::SkipEmptyParts);
        int maxCol = std::max({eastCol, northCol, topoCol});
        
        // Safety check: ensure we have enough columns
        if (tokens.size() <= maxCol) {
            continue;  // Skip invalid lines silently
        }
        
        // Safety check: validate indices
        if (eastCol >= tokens.size() || northCol >= tokens.size() || topoCol >= tokens.size()) {
            continue;
        }
        
        // Try to parse as doubles, skip if invalid
        bool okX, okY, okZ;
        double vx = tokens[eastCol].toDouble(&okX);
        double vy = tokens[northCol].toDouble(&okY);
        double vz = tokens[topoCol].toDouble(&okZ);
        
        if (!okX || !okY || !okZ) {
            continue;  // Skip non-numeric lines (like headers)
        }
        
        x.push_back(vx);
        y.push_back(vy);
        z.push_back(vz);
    }
    
    if (x.empty()) {
        logMessage("No valid data points found in file", "ERROR");
        dataExtentLabel->setText("No valid data found");
        return;
    }
    
    // Calculate extents
    cachedXMin = cachedXMax = x[0];
    cachedYMin = cachedYMax = y[0];
    cachedZMin = cachedZMax = z[0];
    
    for (size_t i = 1; i < x.size(); ++i) {
        cachedXMin = std::min(cachedXMin, x[i]);
        cachedXMax = std::max(cachedXMax, x[i]);
        cachedYMin = std::min(cachedYMin, y[i]);
        cachedYMax = std::max(cachedYMax, y[i]);
        cachedZMin = std::min(cachedZMin, z[i]);
        cachedZMax = std::max(cachedZMax, z[i]);
    }
    
    extentsValid = true;
    logMessage(QString("Scanned %1 valid data points - extents computed").arg(x.size()), "INFO");
}

void MainWindow::browseOutputFile() {
    QString filename = QFileDialog::getSaveFileName(this, "Save Output", "", "PolyTomo Mesh (*.pmesh);;All Files (*)");
    if (!filename.isEmpty()) {
        outputFileEdit->setText(filename);
        currentOutputFile = filename;
        validateInputs();
    }
}

void MainWindow::onPrismSizeChanged() {
    double dx = dxSpinBox->value();
    double dy = dySpinBox->value();
    prismInfoLabel->setText(QString("Prism base: %1 m × %2 m").arg(dx).arg(dy));
    if(estimateUpdateTimer) estimateUpdateTimer->start();
}

void MainWindow::onLayerTableChanged() {
    bool isVariable = (cellTypeCombo->currentIndex() == 1);
    layerTable->setVisible(isVariable);
    addLayerBtn->setVisible(isVariable);
    removeLayerBtn->setVisible(isVariable);
    uniformDzSpinBox->setVisible(!isVariable);
    updateTotalDepth();
}

void MainWindow::addVariableLayer() {
    int row = layerTable->rowCount();
    layerTable->insertRow(row);
    double defThick = (row == 0) ? 25.0 : 50.0;
    layerTable->setItem(row, 0, new QTableWidgetItem(QString::number(defThick)));
    layerTable->setItem(row, 1, new QTableWidgetItem("4"));
    QTableWidgetItem* cumItem = new QTableWidgetItem("0.0");
    cumItem->setFlags(cumItem->flags() & ~Qt::ItemIsEditable);
    cumItem->setBackground(QColor(230, 230, 230));
    layerTable->setItem(row, 2, cumItem);
    updateTotalDepth();
}

void MainWindow::removeVariableLayer() {
    int row = layerTable->currentRow();
    if (row >= 0) {
        layerTable->removeRow(row);
        updateTotalDepth();
    }
}

void MainWindow::updateTotalDepth() {
    int nz = 0;
    double calcDepth = 0.0;
    if (cellTypeCombo->currentIndex() == 0) {
        double dz = uniformDzSpinBox->value();
        double total = totalDepthSpinBox->value();
        nz = static_cast<int>(std::ceil(total / dz));
        calcDepth = nz * dz;
    } else {
        for (int i = 0; i < layerTable->rowCount(); ++i) {
            QTableWidgetItem* tItem = layerTable->item(i, 0);
            QTableWidgetItem* cItem = layerTable->item(i, 1);
            if (tItem && cItem) {
                double t = tItem->text().toDouble();
                int c = cItem->text().toInt();
                calcDepth += t * c;
                nz += c;
                if (layerTable->item(i, 2))
                    layerTable->item(i, 2)->setText(QString::number(calcDepth, 'f', 1));
            }
        }
    }
    if (totalLayersLabel) totalLayersLabel->setText(QString("Total Layers: %1").arg(nz));
    if (calculatedDepthLabel) calculatedDepthLabel->setText(QString("Calculated Depth: %1 m").arg(calcDepth, 0, 'f', 1));
    if(estimateUpdateTimer) estimateUpdateTimer->start();
}

void MainWindow::updateMeshEstimate() {
    if (!extentsValid) {
        dataExtentLabel->setText("Load data file first");
        gridDimensionsLabel->setText("Grid: ? x ? x ?");
        return;
    }
    
    // Get parameters from UI
    double dx = dxSpinBox->value();
    double dy = dySpinBox->value();
    
    // Build dzLayers
    std::vector<double> dzLayers;
    if (cellTypeCombo->currentIndex() == 0) {
        double dz = uniformDzSpinBox->value();
        double totalDepth = totalDepthSpinBox->value();
        int n = static_cast<int>(std::ceil(totalDepth / dz));
        dzLayers.assign(n, dz);
    } else {
        for (int i = 0; i < layerTable->rowCount(); ++i) {
            double t = layerTable->item(i, 0)->text().toDouble();
            int c = layerTable->item(i, 1)->text().toInt();
            for (int k = 0; k < c; ++k)
                dzLayers.push_back(t);
        }
    }
    
    // Calculate grid dimensions using cached extents
    int nx = static_cast<int>((cachedXMax - cachedXMin) / dx) + 1;
    int ny = static_cast<int>((cachedYMax - cachedYMin) / dy) + 1;
    int nz = static_cast<int>(dzLayers.size());
    
    // Update extent display
    QString extentInfo = QString(
        "X: [%1, %2] m\n"
        "Y: [%3, %4] m\n"
        "Z: [%5, %6] m")
        .arg(cachedXMin, 0, 'f', 1)
        .arg(cachedXMax, 0, 'f', 1)
        .arg(cachedYMin, 0, 'f', 1)
        .arg(cachedYMax, 0, 'f', 1)
        .arg(cachedZMin, 0, 'f', 1)
        .arg(cachedZMax, 0, 'f', 1);
    dataExtentLabel->setText(extentInfo);
    
    // Calculate estimates
    long long vertices = (long long)(nx + 1) * (ny + 1) * (nz + 1);
    long long prisms = (long long)nx * ny * nz * 2;
    double memoryMB = (vertices * 24.0 + prisms * 100.0) / (1024.0 * 1024.0);
    
    updateEstimateDisplay(nx, ny, nz, vertices, prisms, memoryMB);
}

void MainWindow::onEstimatesReceived(int nx, int ny, int nz, long long vertices, long long prisms, double memoryMB) {
    Q_UNUSED(vertices);
    updateEstimateDisplay(nx, ny, nz, vertices, prisms, memoryMB);
}

void MainWindow::updateEstimateDisplay(int nx, int ny, int nz, long long vertices, long long prisms, double memoryMB) {
    if (nx > 0) {
        estimatedNx = nx; estimatedNy = ny; estimatedNz = nz;
        estimatedVertices = vertices; estimatedPrisms = prisms; estimatedMemoryMB = memoryMB;

        gridDimensionsLabel->setText(QString("%1 × %2 × %3").arg(nx).arg(ny).arg(nz));
        totalCellsLabel->setText(QString("%L1").arg(prisms));
        memoryEstimateLabel->setText(QString("%L1 MB").arg(memoryMB, 0, 'f', 1));
        
        checkMeshSize(); 
        
        meshEstimateValid = true;
    } else {
        meshEstimateValid = false;
    }
}

void MainWindow::validateInputs() {
    bool inputsValid = !currentInputFile.isEmpty() && !currentOutputFile.isEmpty();
    generateBtn->setEnabled(inputsValid);
    if (generateMPIBtn) generateMPIBtn->setEnabled(inputsValid);
}

void MainWindow::checkMeshSize() {
    if (estimatedPrisms > 1000000) {
        warningLabel->setText("⚠ Large mesh! (>1M cells). Generation may be slow.");
        warningLabel->setStyleSheet("color: orange; font-weight: bold;");
    } else {
        warningLabel->setText("✓ Grid size looks good.");
        warningLabel->setStyleSheet("color: green;");
    }
}

void MainWindow::generateModel() {
    if (currentInputFile.isEmpty() || currentOutputFile.isEmpty()) {
        QMessageBox::warning(this, "Missing Input", "Please select input and output files.");
        return;
    }
    generateBtn->setEnabled(false);
    if(generateMPIBtn) generateMPIBtn->setEnabled(false);
    progressBar->setValue(0);
    
    if (!generator) generator = new ModelGenerator();
    ModelGenerator::Parameters params;
    params.inputFile = currentInputFile;
    params.outputFile = currentOutputFile;
    params.dx = dxSpinBox->value();
    params.dy = dySpinBox->value();
    params.gridOffset = gridOffsetSpinBox->value();
    params.depthBelowSurface = totalDepthSpinBox->value();
    params.topoAware = topoAwareCheckBox->isChecked();
    params.eastingColumn = eastingCombo->currentIndex();
    params.northingColumn = northingCombo->currentIndex();
    params.topoColumn = topoCombo->currentIndex();
    if (cellTypeCombo->currentIndex() == 0) {
        double dz = uniformDzSpinBox->value();
        int n = static_cast<int>(std::ceil(params.depthBelowSurface / dz));
        for(int i=0; i<n; ++i) params.dzLayers.push_back(dz);
    } else {
        for (int i=0; i<layerTable->rowCount(); ++i) {
            double t = layerTable->item(i, 0)->text().toDouble();
            int c = layerTable->item(i, 1)->text().toInt();
            for(int k=0; k<c; ++k) params.dzLayers.push_back(t);
        }
    }
    if (formatPMeshRadio->isChecked()) params.outputFormat = ModelGenerator::Parameters::POLYTOMO_PMESH;
    else if (formatVTKRadio->isChecked()) params.outputFormat = ModelGenerator::Parameters::VTK_VTU;
    else params.outputFormat = ModelGenerator::Parameters::BOTH;
    
    generator->setParameters(params);
    if (!generatorThread) generatorThread = new QThread(this);
    
    connect(generator, &ModelGenerator::progressUpdated, this, &MainWindow::onGenerationProgress);
    connect(generator, &ModelGenerator::logMessage, this, [this](const QString& msg){ logMessage(msg, "INFO"); });
    connect(generator, &ModelGenerator::finished, this, &MainWindow::onGenerationComplete);
    
    generator->moveToThread(generatorThread);
    connect(generatorThread, &QThread::started, generator, &ModelGenerator::process);
    
    generatorThread->start();
    logMessage("Started Serial generation...", "INFO");
}

void MainWindow::onGenerationProgress(int value) {
    progressBar->setValue(value);
}

void MainWindow::onGenerationComplete(bool success, const QString& message) {
    generateBtn->setEnabled(true);
    if (generateMPIBtn) generateMPIBtn->setEnabled(true);
    
    if (generatorThread) {
        generatorThread->quit();
        generatorThread->wait();
    }
    if (success) {
        logMessage("Generation Complete!", "SUCCESS");
        visualizeBtn->setEnabled(true);
        QMessageBox::information(this, "Success", message);
    } else {
        logMessage("Generation Failed: " + message, "ERROR");
        QMessageBox::critical(this, "Error", message);
    }
}

// ============================================================================
// MPI SPECIFIC FUNCTIONS
// ============================================================================

bool MainWindow::writeParametersJSON(const QString& filename) {
    QJsonObject json;
    
    json["inputFile"] = inputFileEdit->text().replace("\\", "/");
    json["outputFile"] = outputFileEdit->text().replace("\\", "/");
    json["dx"] = dxSpinBox->value();
    json["dy"] = dySpinBox->value();
    json["gridOffset"] = gridOffsetSpinBox->value();
    json["depthBelowSurface"] = totalDepthSpinBox->value();
    json["topoAware"] = topoAwareCheckBox->isChecked() ? "true" : "false";
    json["eastingColumn"] = eastingCombo->currentIndex();
    json["northingColumn"] = northingCombo->currentIndex();
    json["topoColumn"] = topoCombo->currentIndex();
    
    // Output format: 0=pmesh, 1=vtu, 2=both
    int outputFormat = 0; // Default pmesh
    if (formatVTKRadio->isChecked()) {
        outputFormat = 1; // VTU only
    } else if (formatBothRadio->isChecked()) {
        outputFormat = 2; // Both formats
    }
    json["outputFormat"] = outputFormat;
    
    // Layers
    QJsonArray layers;
    if (cellTypeCombo->currentIndex() == 0) {
        double dz = uniformDzSpinBox->value();
        int n = static_cast<int>(std::ceil(json["depthBelowSurface"].toDouble() / dz));
        for(int i=0; i<n; ++i) layers.append(dz);
    } else {
        for (int i=0; i<layerTable->rowCount(); ++i) {
            double t = layerTable->item(i, 0)->text().toDouble();
            int c = layerTable->item(i, 1)->text().toInt();
            for(int k=0; k<c; ++k) layers.append(t);
        }
    }
    
    // === CRITICAL FIX: Add layers to JSON (moved outside if/else) ===
    json["layers"] = layers;
    logMessage(QString("Writing %1 layers to JSON").arg(layers.size()), "INFO");
    
    QFile file(filename);
    if (!file.open(QIODevice::WriteOnly)) {
        logMessage("ERROR: Could not write parameters file: " + filename, "ERROR");
        return false;
    }
    
    QJsonDocument doc(json);
    file.write(doc.toJson());
    file.close();
    return true;
}

QString MainWindow::getMPIExecutablePath() {
    // Basic search, extend as needed
    return "mpiexec"; 
}

QString MainWindow::getSolverExecutablePath() {
    // Look in same directory
    QString appDir = QCoreApplication::applicationDirPath();
    QString solver = appDir + "/PolyTomoSolver.exe";
    if (QFile::exists(solver)) return solver;
    return "PolyTomoSolver.exe"; // Hope it's in PATH
}

void MainWindow::generateModelMPI() {
    logMessage("═══════════════════════════════════════════════════════", "HEADER");
    logMessage("Starting MPI Parallel Mesh Generation", "HEADER");
    logMessage("═══════════════════════════════════════════════════════", "HEADER");
    
    if (inputFileEdit->text().isEmpty() || outputFileEdit->text().isEmpty()) {
        QMessageBox::warning(this, "Error", "Invalid inputs.");
        return;
    }
    
    // Write params
    QString tempDir = QDir::tempPath();
    QString paramFile = tempDir + "/polytomo_mpi_params.json";
    if (!writeParametersJSON(paramFile)) {
        logMessage("Failed to write parameters file", "ERROR");
        return;
    }
    
    QString mpiexec = getMPIExecutablePath();
    QString solver = getSolverExecutablePath();
    int procs = mpiProcessesSpinBox->value();
    
    logMessage("MPI Exec: " + mpiexec, "INFO");
    logMessage("Solver: " + solver, "INFO");
    logMessage("Procs: " + QString::number(procs), "INFO");
    
    // Check solver exists
    if (!QFile::exists(solver)) {
        logMessage("ERROR: PolyTomoSolver.exe not found!", "ERROR");
        QMessageBox::critical(this, "Error", "PolyTomoSolver.exe not found at: " + solver);
        return;
    }
    
    // Cleanup old process
    if (mpiProcess) {
        if (mpiProcess->state() == QProcess::Running) {
            mpiProcess->kill();
            mpiProcess->waitForFinished();
        }
        delete mpiProcess;
    }
    
    mpiProcess = new QProcess(this);
    
    // Inherit environment (Intel MPI configured in main.cpp)
    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    if (env.contains("I_MPI_ROOT")) {
        logMessage("Intel MPI environment detected", "SUCCESS");
    } else {
        logMessage("WARNING: Intel MPI environment not configured", "WARNING");
    }
    mpiProcess->setProcessEnvironment(env);
    
    connect(mpiProcess, &QProcess::readyReadStandardOutput, this, &MainWindow::onMPIProcessOutput);
    connect(mpiProcess, &QProcess::readyReadStandardError, this, &MainWindow::onMPIProcessError);
    connect(mpiProcess, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished), 
            this, &MainWindow::onMPIProcessFinished);
            
    QStringList args;
    args << "-n" << QString::number(procs) << solver << paramFile;
    
    // Disable buttons and update MPI button text
    generateBtn->setEnabled(false);
    generateMPIBtn->setEnabled(false);
    generateMPIBtn->setText("⏳ Generating (MPI)...");
    generateMPIBtn->setStyleSheet(
        "QPushButton {"
        "  background-color: #ffc107;"
        "  color: #000000 !important;"
        "  font-weight: bold;"
        "  font-size: 14px;"
        "}"
    );
    visualizeBtn->setEnabled(false);
    progressBar->setRange(0, 0); // Indeterminate mode
    progressBar->setVisible(true);
    
    logMessage("Executing: " + mpiexec + " " + args.join(" "), "INFO");
    logMessage("───────────────────────────────────────────────────────", "INFO");
    
    mpiTimer.start();
    mpiProcess->start(mpiexec, args);
    
    if (!mpiProcess->waitForStarted(5000)) {
        int exitCode = mpiProcess->exitCode();
        logMessage("Failed to start MPI process. Exit code: " + QString::number(exitCode), "ERROR");
        
        QString errorMsg = "Failed to start MPI process!\n\n";
        if (exitCode == -1073741515) {
            errorMsg += "Error: Required DLLs not found\n\n";
            errorMsg += "Solution: Restart your computer and try again.";
        }
        QMessageBox::critical(this, "Error", errorMsg);
        
        generateBtn->setEnabled(true);
        generateMPIBtn->setEnabled(true);
        generateMPIBtn->setText("Generate Mesh (MPI Parallel)");
        generateMPIBtn->setStyleSheet(
            "QPushButton { background-color: #28a745; color: #ffffff !important; font-weight: bold; font-size: 14px; }"
        );
        progressBar->setRange(0, 100);
        progressBar->setVisible(false);
    } else {
        logMessage("MPI process started successfully!", "SUCCESS");
    }
}

void MainWindow::onMPIProcessOutput() {
    if (!mpiProcess) return;
    QString out = mpiProcess->readAllStandardOutput();
    logMessage(out.trimmed(), "MPI");
}

void MainWindow::onMPIProcessError() {
    if (!mpiProcess) return;
    QString err = mpiProcess->readAllStandardError();
    if (!err.trimmed().isEmpty()) logMessage(err.trimmed(), "ERROR");
}

void MainWindow::onMPIProcessFinished(int exitCode, QProcess::ExitStatus exitStatus) {
    progressBar->setRange(0, 100);
    progressBar->setValue(100);
    generateBtn->setEnabled(true);
    generateMPIBtn->setEnabled(true);
    
    // Restore button style
    generateMPIBtn->setText("Generate Mesh (MPI Parallel)");
    generateMPIBtn->setStyleSheet(
        "QPushButton { background-color: #28a745; color: white; font-weight: bold; padding: 8px 16px; border-radius: 4px; }"
        "QPushButton:hover { background-color: #218838; }"
    );
    
    double sec = mpiTimer.elapsed() / 1000.0;
    
    if (exitStatus == QProcess::NormalExit && exitCode == 0) {
        logMessage("MPI Generation Success!", "SUCCESS");
        logMessage(QString("Time: %1 sec").arg(sec), "SUCCESS");
        visualizeBtn->setEnabled(true);
        currentOutputFile = outputFileEdit->text();
        QMessageBox::information(this, "Success", "MPI Generation Complete.");
    } else {
        logMessage("MPI Generation Failed. Exit Code: " + QString::number(exitCode), "ERROR");
        QMessageBox::critical(this, "Error", "MPI Solver Failed.");
    }
    
    progressBar->setVisible(false);
}

// ==========================
// Visualization
// ==========================
void MainWindow::visualizeModel() {
    // 1. Determine the expected VTU filename
    QString vtuFile = currentOutputFile;
    QFileInfo fi(vtuFile);
    
    // If output was set to .pmesh, the generator creates a .vtu file alongside it
    if (fi.suffix().toLower() == "pmesh") {
        vtuFile = fi.path() + "/" + fi.completeBaseName() + ".vtu";
    } 
    // If output has no extension or is already vtu
    else if (fi.suffix().toLower() != "vtu") {
        vtuFile += ".vtu";
    }

    // 2. Check if it exists
    if (!checkFileExists(vtuFile)) {
        QMessageBox::warning(this, "Visualization Error", 
            QString("Could not find VTU file:\n%1\n\n"
                    "If you generated only '.pmesh' format, visualization is not available.\n"
                    "Please select 'VTK' or 'Both' formats and regenerate.").arg(vtuFile));
        return;
    }
    
    // 3. Launch
    if (visualizer) delete visualizer;
    visualizer = new VisualizerX(vtuFile);
    connect(visualizer, &QObject::destroyed, this, &MainWindow::onVisualizationClosed);
    visualizer->show();
    logMessage("Opened Visualizer: " + vtuFile, "INFO");
}

void MainWindow::onVisualizationClosed() {
    visualizer = nullptr;
}

void MainWindow::saveSettings() {
    settings->beginGroup("MainWindow");
    settings->setValue("dx", dxSpinBox->value());
    settings->setValue("dy", dySpinBox->value());
    settings->setValue("totalDepth", totalDepthSpinBox->value());
    settings->setValue("darkMode", isDarkMode);
    settings->setValue("useMPI", useMPICheckbox->isChecked());
    settings->setValue("inputFile", currentInputFile);
    settings->setValue("outputFile", currentOutputFile);
    settings->endGroup();
}

void MainWindow::loadSettings() {
    settings->beginGroup("MainWindow");
    dxSpinBox->setValue(settings->value("dx", 10.0).toDouble());
    dySpinBox->setValue(settings->value("dy", 10.0).toDouble());
    totalDepthSpinBox->setValue(settings->value("totalDepth", 1000.0).toDouble());
    isDarkMode = settings->value("darkMode", false).toBool();
    useMPICheckbox->setChecked(settings->value("useMPI", false).toBool());
    QString inF = settings->value("inputFile").toString();
    if (!inF.isEmpty()) { inputFileEdit->setText(inF); currentInputFile = inF; }
    QString outF = settings->value("outputFile").toString();
    if (!outF.isEmpty()) { outputFileEdit->setText(outF); currentOutputFile = outF; }
    settings->endGroup();
}

bool MainWindow::checkFileExists(const QString& path) {
    return QFile::exists(path);
}

void MainWindow::logMessage(const QString& message, const QString& type) {
    if (!logTextEdit) return;
    QString timestamp = QDateTime::currentDateTime().toString("HH:mm:ss");
    QString color = "black";
    if (isDarkMode) color = "white"; // Basic adaptivity
    
    if (type == "ERROR") color = "#ff4444";
    else if (type == "SUCCESS") color = "#00cc00";
    else if (type == "HEADER") color = "#0078d7";
    else if (type == "MPI") color = "#ff8c00";
    
    QString html = QString("<span style='color:gray'>[%1]</span> <span style='color:%2'>%3</span>")
                   .arg(timestamp).arg(color).arg(message);
                   
    logTextEdit->append(html);
    logTextEdit->verticalScrollBar()->setValue(logTextEdit->verticalScrollBar()->maximum());
}

void MainWindow::toggleTheme() {
    isDarkMode = !isDarkMode;
    applyCurrentTheme();
    logMessage(isDarkMode ? "Switched to Dark Mode" : "Switched to Light Mode", "INFO");
}

void MainWindow::applyCurrentTheme() {
    if (isDarkMode) applyDarkTheme(); else applyLightTheme();
}

void MainWindow::applyLightTheme() {
    // 1. Set Palette (Base colors)
    QPalette p;
    p.setColor(QPalette::Window, QColor(240, 240, 240));
    p.setColor(QPalette::WindowText, Qt::black);
    p.setColor(QPalette::Base, Qt::white);
    p.setColor(QPalette::AlternateBase, QColor(245, 245, 245));
    p.setColor(QPalette::Text, Qt::black);
    p.setColor(QPalette::Button, QColor(240, 240, 240));
    p.setColor(QPalette::ButtonText, Qt::black);
    qApp->setPalette(p);

    // 2. Set Stylesheet
    qApp->setStyleSheet(R"(
        /* Global Reset */
        QWidget {
            background-color: #f0f0f0;
            color: black;
            font-family: 'Segoe UI', 'Roboto', sans-serif;
            font-size: 10pt;
        }
        
        QScrollArea { border: none; background-color: #f0f0f0; }
        QScrollArea > QWidget > QWidget { background-color: #f0f0f0; }
        
        QGroupBox {
            font-weight: bold;
            border: 1px solid #bfbfbf;
            border-radius: 4px;
            margin-top: 20px;
            background-color: white;
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            subcontrol-position: top left;
            left: 10px;
            padding: 0 5px;
            color: #333;
            background-color: transparent;
        }
        
        QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {
            padding: 5px; 
            border: 1px solid #ccc; 
            border-radius: 3px; 
            background: white; 
            color: black; 
            min-height: 25px;
        }
        QLineEdit:read-only { background-color: #e6e6e6; color: #555; }
        
        QPushButton {
            padding: 6px 15px; 
            background-color: #0078d7; 
            color: white; 
            border-radius: 4px; 
            font-weight: bold; 
            min-height: 25px;
            border: none;
        }
        QPushButton:hover { background-color: #005a9e; }
        QPushButton:disabled { background-color: #ccc; color: #888; }
        
        /* Table & Text Edit - Light Mode */
        QTextEdit, QTableWidget { 
            background-color: white; 
            color: black; 
            border: 1px solid #ccc; 
            gridline-color: #e0e0e0;
        }
        QHeaderView::section { 
            background-color: #e0e0e0; 
            padding: 4px; 
            border: 1px solid #ccc; 
            color: black;
        }
        QTableCornerButton::section {
            background-color: #e0e0e0;
            border: 1px solid #ccc;
        }
    )");
}

void MainWindow::applyDarkTheme() {
    // 1. Set Palette
    QPalette p;
    p.setColor(QPalette::Window, QColor(45, 45, 45));
    p.setColor(QPalette::WindowText, QColor(220, 220, 220));
    p.setColor(QPalette::Base, QColor(30, 30, 30));
    p.setColor(QPalette::AlternateBase, QColor(35, 35, 35));
    p.setColor(QPalette::Text, QColor(220, 220, 220));
    p.setColor(QPalette::Button, QColor(45, 45, 45));
    p.setColor(QPalette::ButtonText, QColor(220, 220, 220));
    qApp->setPalette(p);

    // 2. Set Stylesheet
    qApp->setStyleSheet(R"(
        /* Global Reset */
        QWidget {
            background-color: #2d2d2d;
            color: #e0e0e0;
            font-family: 'Segoe UI', 'Roboto', sans-serif;
            font-size: 10pt;
        }
        
        QScrollArea { border: none; background-color: #2d2d2d; }
        QScrollArea > QWidget > QWidget { background-color: #2d2d2d; }
        
        QGroupBox {
            font-weight: bold; 
            border: 1px solid #555; 
            border-radius: 4px; 
            margin-top: 20px; 
            background-color: #363636; 
        }
        QGroupBox::title {
            subcontrol-origin: margin; 
            subcontrol-position: top left; 
            left: 10px; 
            padding: 0 5px; 
            color: #4da6ff;
            background-color: transparent;
        }
        
        QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {
            padding: 5px; 
            border: 1px solid #555; 
            border-radius: 3px; 
            background: #252525; 
            color: #e0e0e0; 
            min-height: 25px;
        }
        QLineEdit:read-only { background-color: #1f1f1f; color: #888; }
        
        QPushButton {
            padding: 6px 15px; 
            background-color: #0078d7; 
            color: white; 
            border-radius: 4px; 
            font-weight: bold; 
            min-height: 25px;
            border: none;
        }
        QPushButton:hover { background-color: #0063b1; }
        QPushButton:disabled { background-color: #444; color: #777; }
        
        /* --- FIX: Table Dark Mode --- */
        QTextEdit, QTableWidget { 
            background-color: #252525; 
            color: #e0e0e0; 
            border: 1px solid #555;
            gridline-color: #444; /* Fixes light grid lines */
        }
        /* Fix Table Items */
        QTableWidget::item {
            background-color: #252525;
            color: #e0e0e0;
        }
        /* Fix Headers */
        QHeaderView::section { 
            background-color: #444; 
            color: #e0e0e0; 
            padding: 4px; 
            border: 1px solid #555; 
        }
        /* Fix the Top-Left Corner (The Light Strip) */
        QTableCornerButton::section {
            background-color: #444;
            border: 1px solid #555;
        }
        QMenuBar { background-color: #2d2d2d; color: #e0e0e0; }
        QMenuBar::item:selected { background-color: #555; }
        QMenu { background-color: #333; color: #ddd; border: 1px solid #555; }
        QMenu::item:selected { background-color: #0078d7; }
    )");
}