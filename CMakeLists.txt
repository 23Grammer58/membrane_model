cmake_minimum_required(VERSION 3.16)

set(AVSYM_MAJOR_VERSION 0)
set(AVSYM_MINOR_VERSION 2)
set(AVSYM_PATCH_VERSION 0)
set(AVSYM_VERSION "${AVSYM_MAJOR_VERSION}.${AVSYM_MINOR_VERSION}.${AVSYM_PATCH_VERSION}")
set(PACKAGE_NAME "AVSim")
string(TOLOWER ${PACKAGE_NAME} package_NAME)

project(${PACKAGE_NAME} VERSION ${AVSYM_VERSION} DESCRIPTION "Framework for simulation of mechanical behavior thin membranes in application to aortic valve closure" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread ")
link_libraries(stdc++fs)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" "${CMAKE_MODULE_PATH}")

option(USE_MAGNUM_GUI "Compile with Magnum, Corrade and ImGui for debug view" OFF)
option(USE_INMOST "Compile with INMOST solvers" ON)
#option(USE_MATIO "Compile with support of .mat files" OFF)
option(USE_DATADRIVEN_FORCES "Compile with Data-Driven support module" OFF)
option(USE_COLLISION "Compile with support of collision module" OFF)

option(DOWNLOAD_Ani3d "Attempts to download Ani3d if library not found" OFF)
option(DOWNLOAD_eigen3 "Attempts to download Eigen3 if library not found" OFF)
option(DOWNLOAD_casadi "Attempts to download casadi if library not found" OFF)
option(DOWNLOAD_rpoly "Attempts to download RpolyPlusPlus if library not found" ON)
option(DOWNLOAD_inmost "Attempts to download INMOST if library not found" OFF)
option(DOWNLOAD_sundials "Attempts to download SUNDIALS if library not found" OFF)
option(DOWNLOAD_CGAL "Attempts to download CGAL if library not found" OFF)
option(DOWNLOAD_GSL "Attempts to download GSL if library not found" OFF)
option(DOWNLOAD_nlohmann_json "Attempts to download nlohmann_json if library not found" OFF)
#option(DOWNLOAD_Bullet "Attempts to download Bullet (BulletCollision, LinearMath components) if library not found" OFF)
#option(DOWNLOAD_Boost "Attempts to download Boost::serialization, Boost::filesystem components if library not found" OFF)
#option(DOWNLOAD_GUI "Attempts to download Corrade, Magnum, MagnumIntegration and ImGui libraries for support of GUI if these not found" OFF)
set(LIB_DOWNLOAD_PATH "${CMAKE_BINARY_DIR}/external_packages")

include(GNUInstallDirs)
if (NOT DOWNLOAD_Ani3d)
    find_package(Ani3d REQUIRED)
else()
    find_package(Ani3d QUIET)
    if (NOT Ani3d_FOUND)
        include(cmake/DownloadAni3d.cmake)
    endif()
endif()
list(APPEND LIB_NAMES Ani3d::FRTPRM)

if (NOT DOWNLOAD_eigen3)
    find_package (Eigen3  3.3.7 REQUIRED )
else()
    find_package (Eigen3  3.3.7 QUIET )
    if (NOT Eigen3_FOUND)
        include(cmake/Downloadeigen3.cmake)
    endif()
endif()	
list(APPEND LIB_NAMES Eigen3::Eigen)

if (NOT DOWNLOAD_casadi)
    find_package(casadi REQUIRED)
else()
    find_package(casadi QUIET)
    if (NOT casadi_FOUND)
        include(cmake/Downloadcasadi.cmake)
    endif()
endif()
list(APPEND LIB_NAMES casadi)

if (NOT DOWNLOAD_GSL)
    find_package(GSL REQUIRED)  # GPLv3
else()
    find_package(GSL QUIET)
    if (NOT GSL_FOUND)
        include(cmake/Downloadgsl.cmake)
    endif()    
endif()
list(APPEND LIB_NAMES GSL::gsl)

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
if (NOT DOWNLOAD_CGAL)
    find_package (CGAL 5.3 REQUIRED)   # GPLv3
else()
    find_package (CGAL 5.3 QUIET)
    if (NOT CGAL_FOUND)
        include(cmake/Downloadcgal.cmake)
    endif()
endif() 
list(APPEND LIB_NAMES CGAL::CGAL )

if (NOT DOWNLOAD_sundials)
    find_package(SUNDIALS COMPONENTS kinsol REQUIRED)
else()
    find_package(SUNDIALS COMPONENTS kinsol QUIET)
    if (NOT SUNDIALS_FOUND)
        include(cmake/Downloadsundials.cmake)
    endif()
endif()
list(APPEND LIB_NAMES SUNDIALS::kinsol)

if (USE_INMOST)
    if(NOT DOWNLOAD_inmost)
        find_package(INMOST REQUIRED)
    else()
        find_package(INMOST QUIET)
        if (NOT INMOST_FOUND)
            include(cmake/Downloadinmost.cmake)
        endif()
    endif()
    target_compile_definitions(inmost INTERFACE -DUSE_INMOST_SOLVERS ${INMOST_DEFINITIONS})
    target_link_directories(inmost INTERFACE ${INMOST_LIBRARY_DIRS})
    target_link_libraries(inmost INTERFACE ${INMOST_LIBRARIES})
    target_include_directories(inmost INTERFACE ${INMOST_INCLUDE_DIRS})

    list(APPEND SOURCE_FILES AVSim/Solvers/LinearSolverInmost.cpp)
    list(APPEND LIB_NAMES inmost)
endif()

if (USE_COLLISION)
    find_package(Bullet REQUIRED  BulletCollision LinearMath)
    add_library(Bullet_Math INTERFACE)
    target_include_directories(Bullet_Math INTERFACE ${BULLET_INCLUDE_DIR})
    target_link_libraries(Bullet_Math INTERFACE ${BULLET_MATH_LIBRARY})

    add_library(Bullet_Collision INTERFACE)
    target_include_directories(Bullet_Collision INTERFACE ${BULLET_INCLUDE_DIR})
    target_link_libraries(Bullet_Collision INTERFACE ${BULLET_COLLISION_LIBRARY} Bullet_Math)

    list(APPEND LIB_NAMES Bullet_Collision)
    set(AVSIM_TARGETS ${AVSIM_TARGETS} Bullet_Math Bullet_Collision)
endif()

if (NOT DOWNLOAD_nlohmann_json)
    find_package(nlohmann_json 3.2.0 REQUIRED)
else()
    find_package(nlohmann_json QUIET)
    if (NOT nlohmann_json_FOUND)
        include(cmake/Downloadnlohmann_json.cmake)
    endif()
endif()    
list(APPEND LIB_NAMES nlohmann_json::nlohmann_json )

if (NOT DOWNLOAD_rpoly)
    find_package(rpoly REQUIRED)
else()
    find_package(rpoly QUIET)
    if (NOT rpoly_FOUND)
        include(cmake/DownloadRpoly.cmake)
    endif()
endif()
list(APPEND LIB_NAMES rpoly_plus_plus)

if (USE_DATADRIVEN_FORCES)
    find_package(Boost 1.70.0 COMPONENTS serialization filesystem program_options REQUIRED)
    list(APPEND LIB_NAMES Boost::serialization Boost::filesystem Boost::program_options)
endif()

if (USE_MAGNUM_GUI)
    add_definitions(-DUSE_MAGNUM_GUI)
    #set (MAGNUM_DIR "/home/alex/Desktop/prog/libs/")
    #set(CMAKE_MODULE_PATH  "${MAGNUM_DIR}/magnum-integration/modules/" ${CMAKE_MODULE_PATH})
#    set(IMGUI_DIR "/home/alex/Desktop/prog/libs/imgui")
#    set(CMAKE_MODULE_PATH  "/home/alex/Desktop/prog/libs//magnum-integration/modules/" ${CMAKE_MODULE_PATH})
    set(CMAKE_MODULE_PATH /home/alex/Desktop/Work/libs/GUI/magnum-integration/modules ${CMAKE_MODULE_PATH})
    #set(IMGUI_DIR "/home/alex/Desktop/Work/libs/GUI/magnum-integration/install/include/MagnumExternal/ImGui")
    #set(WITH_IMGUI ON CACHE BOOL "" FORCE)
    #set(WITH_ANYIMAGEIMPORTER ON CACHE BOOL "" FORCE)
    #set(WITH_STBIMAGEIMPORTER ON CACHE BOOL "" FORCE)

    find_package(Corrade REQUIRED Main)
    find_package(Magnum REQUIRED
            GL
            MeshTools
            Shaders
            Primitives
            SceneGraph
            Trade
            Sdl2Application
            )
    find_package(MagnumIntegration REQUIRED ImGui)
    list(APPEND SOURCE_FILES AVSim/GuiApplication.cpp)
    list(APPEND LIB_NAMES
            MagnumIntegration::ImGui
            Magnum::Application
            Magnum::MeshTools
            #Magnum::Trade
            Magnum::GL
            Magnum::Magnum
            Magnum::SceneGraph
            Magnum::Primitives
            Magnum::Shaders
            Corrade::Main
            )
endif()

function(prefixed_install files dir_prefix)
    #message("files = ${files}")
    foreach(file ${files})
        get_filename_component(parent_dir ${file} DIRECTORY)
        #message("file = ${file} pardir = ${parent_dir}")
        install(FILES ${file} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${dir_prefix}/${parent_dir})
    endforeach(file ${files})
endfunction()
add_subdirectory("AVSim/")
add_subdirectory("examples/")
add_subdirectory("benchmarks/")

install(TARGETS ${AVSIM_TARGETS}
        EXPORT ${package_NAME}-targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
set(INSTALL_CMAKE_FINDS ${CMAKE_INSTALL_LIBDIR}/cmake_find)      
install(FILES cmake/FindAni3d.cmake cmake/Findrpoly.cmake cmake/FindFortran.cmake DESTINATION ${INSTALL_CMAKE_FINDS})        

install(EXPORT "${package_NAME}-targets" DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${package_NAME}")
include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/${package_NAME}-config.cmake"
        INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${package_NAME}"
        PATH_VARS INSTALL_CMAKE_FINDS
        )
# generate the version file for the config file
write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/${package_NAME}ConfigVersion.cmake"
        VERSION "${PACKAGE_VERSION}"
        COMPATIBILITY SameMinorVersion
)
# install the generated configuration files
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${package_NAME}-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${package_NAME}ConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${package_NAME}"
        )       
