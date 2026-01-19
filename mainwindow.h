#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QWidget>
#include <QPushButton>
#include <QLineEdit>
#include <QLabel>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QTableWidget>
#include <QTextEdit>
#include <QProgressBar>
#include <QGroupBox>
#include <QCheckBox>
#include <QRadioButton>
#include <QButtonGroup>
#include <QTabWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QFileDialog>
#include <QMessageBox>
#include <QThread>
#include <QSettings>
#include <QTimer>
#include <QProcess>
#include <QElapsedTimer>
#include <vector>
#include <string>

class ModelGenerator;
class VisualizerX;

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    // File operations
    void browseInputFile();
    void browseOutputFile();
    void onInputFileLoaded();

    // Grid configuration
    void onPrismSizeChanged();
    void addVariableLayer();
    void removeVariableLayer();
    void updateTotalDepth();
    void onLayerTableChanged();

    // Real-time estimation
    void updateMeshEstimate();
    void onEstimatesReceived(int nx, int ny, int nz, long long vertices, long long prisms, double memoryMB);

    // Generation
    void generateModel();
    void onGenerationProgress(int value);
    void onGenerationComplete(bool success, const QString& message);

    // MPI Generation (NEW!)
    void generateModelMPI();
    void onMPIProcessOutput();
    void onMPIProcessError();
    void onMPIProcessFinished(int exitCode, QProcess::ExitStatus exitStatus);

    // Visualization
    void visualizeModel();
    void onVisualizationClosed();

    // Settings
    void saveSettings();
    void loadSettings();

    // Validation
    void validateInputs();
    void checkMeshSize();
    
    // Theme
    void toggleTheme();

private:
    void setupUI();
    void createMenuBar();
    void createInputSection();
    void createGridSection();
    void createDepthSection();
    void createOutputSection();
    void createMPISection();
    void createActionSection();
    void createLogSection();
    void createEstimatePanel();

    void logMessage(const QString& message, const QString& type = "INFO");
    bool checkFileExists(const QString& path);
    void displayDataPreview();
    void scanInputFile();  // Scan input file for extents
    void updateEstimateDisplay(int nx, int ny, int nz, long long vertices, long long prisms, double memoryMB);
    
    // Theme functions
    void applyLightTheme();
    void applyDarkTheme();
    void applyCurrentTheme();

    // MPI helper functions (NEW!)
    bool writeParametersJSON(const QString& filename);
    QString getMPIExecutablePath();
    QString getSolverExecutablePath();

    // --- Input Section ---
    QGroupBox* inputGroup;
    QLineEdit* inputFileEdit;
    QLineEdit* outputFileEdit;
    QPushButton* browseInputBtn;
    QPushButton* browseOutputBtn;
    QLabel* dataPointsLabel;
    QLabel* dataExtentLabel;
    QComboBox* eastingCombo;
    QComboBox* northingCombo;
    QComboBox* topoCombo;
    QLabel* validationLabel;

    // --- Grid Section ---
    QGroupBox* gridGroup;
    QDoubleSpinBox* dxSpinBox;
    QDoubleSpinBox* dySpinBox;
    QDoubleSpinBox* gridOffsetSpinBox;
    QLabel* prismInfoLabel;

    // --- Depth Section ---
    QGroupBox* depthGroup;
    QComboBox* cellTypeCombo;
    QDoubleSpinBox* uniformDzSpinBox;
    QDoubleSpinBox* totalDepthSpinBox;
    QTableWidget* layerTable;
    QPushButton* addLayerBtn;
    QPushButton* removeLayerBtn;
    QLabel* totalLayersLabel;
    QLabel* calculatedDepthLabel;
    QCheckBox* topoAwareCheckBox;

    // --- Mesh Estimate Panel ---
    QGroupBox* estimateGroup;
    QLabel* gridDimensionsLabel;
    QLabel* vertexCountLabel;
    QLabel* prismCountLabel;
    QLabel* totalCellsLabel;
    QLabel* memoryEstimateLabel;
    QLabel* warningLabel;

    // --- Output Format Section ---
    QGroupBox* outputGroup;
    QRadioButton* formatPMeshRadio;
    QRadioButton* formatVTKRadio;
    QRadioButton* formatBothRadio;
    QButtonGroup* formatButtonGroup;

    // --- MPI Section ---
    QGroupBox* mpiGroup;
    QCheckBox* useMPICheckbox;
    QSpinBox* mpiProcessesSpinBox;
    QLabel* mpiWarningLabel;

    // --- Actions Section ---
    QGroupBox* actionGroup;
    QPushButton* generateBtn;
    QPushButton* generateMPIBtn;  // NEW!
    QPushButton* visualizeBtn;
    QProgressBar* progressBar;

    // --- Log Section ---
    QGroupBox* logGroup;
    QTextEdit* logTextEdit;
    QPushButton* clearLogBtn;
    QPushButton* saveLogBtn;

    // --- Data Management ---
    ModelGenerator* generator;
    VisualizerX* visualizer;
    QThread* generatorThread;
    QSettings* settings;
    QTimer* estimateUpdateTimer;

    // MPI Process Management (NEW!)
    QProcess* mpiProcess;
    QElapsedTimer mpiTimer;
    int mpiProcessCount;

    QString currentInputFile;
    QString currentOutputFile;
    int totalDataPoints;
    double dataExtentX;
    double dataExtentY;

    int estimatedNx;
    int estimatedNy;
    int estimatedNz;
    long long estimatedVertices;
    long long estimatedPrisms;
    double estimatedMemoryMB;

    bool inputFileValid;
    bool meshEstimateValid;
    
    // Cached data extents (computed on file load)
    double cachedXMin = 0, cachedXMax = 0;
    double cachedYMin = 0, cachedYMax = 0;
    double cachedZMin = 0, cachedZMax = 0;
    bool extentsValid = false;
    
    // Theme state
    bool isDarkMode;
};

#endif // MAINWINDOW_H
