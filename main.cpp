#include <QApplication>
#include <QStyleFactory>
#include <QSurfaceFormat>
#include <QCoreApplication>
#include <QFont>
#include <QDir>
#include <QMessageBox>
#include <QDebug>

// VTK OpenGL widget support
#include <QVTKOpenGLNativeWidget.h>
#include <vtkAutoInit.h>

#include "mainwindow.h"

// Ensure required VTK modules are initialized
VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle)
VTK_MODULE_INIT(vtkRenderingFreeType)

/**
 * Auto-configure Intel MPI environment
 * This runs before GUI starts so users don't need to run setvars.bat
 */
bool configureIntelMPI()
{
    // Possible Intel MPI installation paths
    QStringList possiblePaths = {
        "C:/Program Files (x86)/Intel/oneAPI/mpi/2021.17",
        "C:/Program Files (x86)/Intel/oneAPI/mpi/2021.16",
        "C:/Program Files (x86)/Intel/oneAPI/mpi/latest",
        "C:/Program Files/Intel/oneAPI/mpi/2021.17",
        "C:/Program Files/Intel/oneAPI/mpi/2021.16",
        "C:/Program Files/Intel/oneAPI/mpi/latest"
    };
    
    QString intelMPIRoot;
    
    // Find Intel MPI installation
    for (const QString& path : possiblePaths) {
        if (QDir(path).exists()) {
            intelMPIRoot = path;
            qInfo() << "Found Intel MPI at:" << path;
            break;
        }
    }
    
    if (intelMPIRoot.isEmpty()) {
        qWarning() << "Intel MPI not found. MPI features will be disabled.";
        qWarning() << "To enable MPI, install Intel oneAPI MPI from:";
        qWarning() << "https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html";
        return false;
    }
    
    // Build paths using native separators for Windows
    QString binPath = QDir::toNativeSeparators(intelMPIRoot + "/bin");
    QString binReleasePath = QDir::toNativeSeparators(intelMPIRoot + "/bin/release");
    QString libfabricPath = QDir::toNativeSeparators(intelMPIRoot + "/libfabric/bin");
    
    // Get current PATH
    QString currentPath = QString::fromLocal8Bit(qgetenv("PATH"));
    
    // Prepend Intel MPI paths to PATH
    QStringList pathList;
    pathList << binPath;
    pathList << binReleasePath;
    pathList << libfabricPath;
    
    if (!currentPath.isEmpty()) {
        pathList << currentPath;
    }
    
    QString newPath = pathList.join(";");
    
    // Set environment variables for this process
    qputenv("PATH", newPath.toLocal8Bit());
    qputenv("I_MPI_ROOT", QDir::toNativeSeparators(intelMPIRoot).toLocal8Bit());
    qputenv("I_MPI_FABRICS", "shm");  // Use shared memory for single node (fastest)
    qputenv("I_MPI_PIN", "1");        // Enable process pinning for better performance
    
    qInfo() << "Intel MPI configured successfully";
    qInfo() << "I_MPI_ROOT:" << intelMPIRoot;
    qInfo() << "MPI features enabled";
    
    return true;
}

/**
 * Check if MPI solver exists
 */
bool checkMPISolver()
{
    QString appDir = QCoreApplication::applicationDirPath();
    QString solverPath = appDir + "/PolyTomoSolver.exe";
    
    if (QFile::exists(solverPath)) {
        qInfo() << "MPI Solver found:" << solverPath;
        return true;
    } else {
        qWarning() << "MPI Solver not found:" << solverPath;
        qWarning() << "MPI parallel generation will be disabled";
        return false;
    }
}

int main(int argc, char *argv[])
{
    // ----------------------------------------------------------------
    // 1. High DPI & Scaling Fixes (Must be before QApplication)
    // ----------------------------------------------------------------
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    QCoreApplication::setAttribute(Qt::AA_UseHighDpiPixmaps);
    
    #if QT_VERSION >= QT_VERSION_CHECK(5, 14, 0)
    QGuiApplication::setHighDpiScaleFactorRoundingPolicy(Qt::HighDpiScaleFactorRoundingPolicy::PassThrough);
    #endif

    // 2. VTK Surface Format
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());

    QApplication app(argc, argv);

    // 3. Application Metadata
    QCoreApplication::setOrganizationName("PolyTomo");
    QCoreApplication::setOrganizationDomain("polytomo.org");
    QCoreApplication::setApplicationName("PolyTomo MeshGen");
    QCoreApplication::setApplicationVersion("1.0");

    // 4. Style & Font
    QApplication::setStyle(QStyleFactory::create("Fusion"));
    
    // Increase font size slightly to prevent squashing on high-res screens
    QFont font = QApplication::font();
    font.setPointSize(10); 
    QApplication::setFont(font);

    // ----------------------------------------------------------------
    // 5. PRODUCTION FIX: Auto-configure Intel MPI
    // ----------------------------------------------------------------
    qInfo() << "==========================================";
    qInfo() << "PolyTomo MeshGen - Starting";
    qInfo() << "==========================================";
    
    bool mpiAvailable = configureIntelMPI();
    bool solverAvailable = checkMPISolver();
    
    if (mpiAvailable && solverAvailable) {
        qInfo() << "Status: MPI features ENABLED";
    } else {
        qInfo() << "Status: MPI features DISABLED (serial mode only)";
        if (!mpiAvailable) {
            qInfo() << "Reason: Intel MPI not found";
        }
        if (!solverAvailable) {
            qInfo() << "Reason: PolyTomoSolver.exe not found";
        }
    }
    
    qInfo() << "==========================================";

    // ----------------------------------------------------------------
    // 6. Launch Main Window
    // ----------------------------------------------------------------
    MainWindow window;
    window.show();

    return app.exec();
}
