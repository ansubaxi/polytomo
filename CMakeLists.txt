cmake_minimum_required(VERSION 3.16)
project(PolyTomoMeshGen VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Find Qt5
find_package(Qt5 REQUIRED COMPONENTS Core Widgets Concurrent)

# Find VTK
# FIX: Added FiltersGeometry (for vtkGeometryFilter) and InteractionWidgets (for axes widget)
find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    CommonColor
    CommonTransforms
    FiltersSources
    FiltersGeneral
    FiltersGeometry      # <--- ADDED: Fixes LNK2019 vtkGeometryFilter
    RenderingCore
    RenderingOpenGL2
    InteractionStyle
    InteractionWidgets   # <--- ADDED: Required for OrientationMarkerWidget
    RenderingFreeType
    RenderingAnnotation
    GUISupportQt
    IOImage
    IOXML
)

# =======================================================================
# MPI DETECTION
# =======================================================================
find_package(MPI)
if(MPI_FOUND)
    add_definitions(-DHAVE_MPI)
    message(STATUS "MPI Found. Parallel solver enabled.")
else()
    message(WARNING "MPI Not Found. Parallel solver disabled.")
endif()

# =======================================================================
# GUI EXECUTABLE
# =======================================================================
set(SOURCES
    main.cpp
    mainwindow.cpp
    modelgenerator.cpp
    visualizerx.cpp
)

set(HEADERS
    mainwindow.h
    modelgenerator.h
    visualizerx.h
)

add_executable(${PROJECT_NAME} WIN32 ${SOURCES} ${HEADERS})

include_directories(${VTK_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME}
    PRIVATE Qt5::Core Qt5::Widgets Qt5::Concurrent
    PRIVATE ${VTK_LIBRARIES}
)

if(MPI_FOUND)
    target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_CXX)
endif()

# VTK Auto-Init (Critical for factories like RenderingOpenGL2)
if(VTK_VERSION VERSION_GREATER_EQUAL "9.0.0")
    vtk_module_autoinit(
        TARGETS ${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
    )
endif()

target_compile_definitions(${PROJECT_NAME} PRIVATE USE_VTK)

# =======================================================================
# SOLVER EXECUTABLE (CLI)
# =======================================================================
if(MPI_FOUND)
    add_executable(PolyTomoSolver solver_main.cpp)
    target_link_libraries(PolyTomoSolver PRIVATE MPI::MPI_CXX)
    install(TARGETS PolyTomoSolver RUNTIME DESTINATION bin)
endif()

# =======================================================================
# DEPLOYMENT (Windows)
# =======================================================================
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W3) # Warning level 3
    
    # Qt Windeployqt helper
    if(Qt5_FOUND AND WIN32 AND TARGET Qt5::qmake AND NOT TARGET Qt5::windeployqt)
        get_target_property(_qt5_qmake_location Qt5::qmake IMPORTED_LOCATION)
        execute_process(
            COMMAND "${_qt5_qmake_location}" -query QT_INSTALL_PREFIX
            RESULT_VARIABLE return_code
            OUTPUT_VARIABLE qt5_install_prefix
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(EXISTS "${qt5_install_prefix}/bin/windeployqt.exe")
            add_executable(Qt5::windeployqt IMPORTED)
            set_target_properties(Qt5::windeployqt PROPERTIES IMPORTED_LOCATION "${qt5_install_prefix}/bin/windeployqt.exe")
        endif()
    endif()
    
    if(TARGET Qt5::windeployqt)
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND Qt5::windeployqt --no-translations "$<TARGET_FILE:${PROJECT_NAME}>"
        )
    endif()
endif()

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)