include(FetchContent)
include(ExternalProject)

message("-- Download CGAL")
FetchContent_Declare(cgal_get
        URL https://github.com/CGAL/cgal/releases/download/v5.3.2/CGAL-5.3.2-library.tar.xz
        UPDATE_DISCONNECTED TRUE
        PREFIX "${LIB_DOWNLOAD_PATH}"
        SOURCE_DIR "${LIB_DOWNLOAD_PATH}/CGAL"
        )
FetchContent_GetProperties(cgal_get)
if(NOT cgal_get_POPULATED)
    FetchContent_Populate(cgal_get)
endif()        

file(MAKE_DIRECTORY "${LIB_DOWNLOAD_PATH}/CGAL/install")
file(MAKE_DIRECTORY "${LIB_DOWNLOAD_PATH}/CGAL/install/${CMAKE_INSTALL_INCLUDEDIR}")
file(MAKE_DIRECTORY "${LIB_DOWNLOAD_PATH}/CGAL/install/${CMAKE_INSTALL_INCLUDEDIR}/CGAL")
file(MAKE_DIRECTORY "${LIB_DOWNLOAD_PATH}/CGAL/install/${CMAKE_INSTALL_DATAROOTDIR}")

enable_language(C)
ExternalProject_Add(cgal_get
        PREFIX "${LIB_DOWNLOAD_PATH}"
        SOURCE_DIR "${cgal_get_SOURCE_DIR}"
        BINARY_DIR "${cgal_get_SOURCE_DIR}/build"
        INSTALL_DIR "${cgal_get_SOURCE_DIR}/install"
        CMAKE_ARGS  "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
                    "-DCMAKE_C_COMPILER:PATH=${CMAKE_C_COMPILER}"
                    "-DCMAKE_INSTALL_PREFIX=${cgal_get_SOURCE_DIR}/install"
                    "-DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}"
                    "-DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}"
                    "-DCMAKE_INSTALL_DATAROOTDIR=${CMAKE_INSTALL_DATAROOTDIR}"
        )
set (CGAL_DEFINITIONS  "")
set (CGAL_INCLUDE_DIR  "${LIB_DOWNLOAD_PATH}/CGAL/install/${CMAKE_INSTALL_INCLUDEDIR}")
set (CGAL_INCLUDE_DIRS "${LIB_DOWNLOAD_PATH}/CGAL/install/${CMAKE_INSTALL_INCLUDEDIR}")
set (CGAL_ROOT_DIR     "${LIB_DOWNLOAD_PATH}/CGAL/install/") 

add_library(CGAL::CGAL INTERFACE IMPORTED GLOBAL)        
add_dependencies(CGAL::CGAL cgal_get)
target_include_directories(CGAL::CGAL INTERFACE ${CGAL_INCLUDE_DIRS})
if (CGAL_DEFINITIONS)
    target_compile_definitions(CGAL::CGAL ${CGAL_DEFINITIONS})
endif()    
set(cgal_DOWNLOADED TRUE)
install(DIRECTORY "${LIB_DOWNLOAD_PATH}/CGAL/install/" 
        DESTINATION "${CMAKE_INSTALL_PREFIX}") 

