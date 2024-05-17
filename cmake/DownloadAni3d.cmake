include(FetchContent)
include(ExternalProject)

message("-- Download Ani3d")
FetchContent_Declare(ani3d_get
        URL         https://sourceforge.net/projects/ani3d/files/ani3d/ani3d-3.1/ani3D-3.1.tar.gz
        URL_HASH    MD5=7d777454d138c2e987e11c519c68aded
        UPDATE_DISCONNECTED TRUE
        PREFIX ${LIB_DOWNLOAD_PATH}
        SOURCE_DIR ${LIB_DOWNLOAD_PATH}/Ani3d
        )
FetchContent_GetProperties(ani3d_get)
if(NOT ani3d_get_POPULATED)
    FetchContent_Populate(ani3d_get)
endif()

enable_language(C)
ExternalProject_Add(ani3d_get
        PREFIX ${LIB_DOWNLOAD_PATH}
        SOURCE_DIR ${ani3d_get_SOURCE_DIR}
        BINARY_DIR ${ani3d_get_SOURCE_DIR}/build
        INSTALL_DIR ${ani3d_get_SOURCE_DIR}
        CMAKE_ARGS  "-DCMAKE_BUILD_TYPE=Release"
                    "-DBUILD_TESTING=OFF"
                    "-DCMAKE_C_COMPILER:PATH=${CMAKE_C_COMPILER}"
                    "-DCMAKE_Fortran_FLAGS=-fallow-argument-mismatch"
        )
if( WIN32 )
    set(LIB_TYPE lib)
    set(LIB_PREF "")
else( WIN32 )
    set(LIB_TYPE a)
    set(LIB_PREF lib)
endif( WIN32 )

set(Ani3d_AFT_LIBRARY "${ani3d_get_SOURCE_DIR}/lib/${LIB_PREF}aft3D-3.0.${LIB_TYPE}")
set(Ani3d_MBA_LIBRARY "${ani3d_get_SOURCE_DIR}/lib/${LIB_PREF}mba3D-3.0.${LIB_TYPE}")
set(Ani3d_FRTPRM_LIBRARY "${ani3d_get_SOURCE_DIR}/lib/${LIB_PREF}frtprm-3.0.${LIB_TYPE}")
set(Ani3d_AFT_INCLUDE_DIR ${ani3d_get_SOURCE_DIR}/include)
set(Ani3d_AFT_INCLUDE_SRC ${ani3d_get_SOURCE_DIR}/src)
set(Ani3d_DOWNLOADED TRUE)
find_package(Fortran REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
add_library(Ani3d::AFT UNKNOWN IMPORTED GLOBAL)
set_target_properties(Ani3d::AFT PROPERTIES IMPORTED_LOCATION ${Ani3d_AFT_LIBRARY})
target_include_directories(Ani3d::AFT INTERFACE ${Ani3d_AFT_INCLUDE_DIR})
set(Ani3d_DEPENDENCY_LIBS "${LAPACK_LIBRARIES};${BLAS_LIBRARIES};${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES};m")
target_link_libraries(Ani3d::AFT INTERFACE "\$<LINK_ONLY:${Ani3d_DEPENDENCY_LIBS}>")
add_dependencies(Ani3d::AFT ani3d_get)

add_library(Ani3d::MBA UNKNOWN IMPORTED GLOBAL)
set_target_properties(Ani3d::MBA PROPERTIES IMPORTED_LOCATION ${Ani3d_MBA_LIBRARY})
target_link_libraries(Ani3d::MBA INTERFACE Ani3d::AFT)

add_library(Ani3d::FRTPRM UNKNOWN IMPORTED GLOBAL)
set_target_properties(Ani3d::FRTPRM PROPERTIES IMPORTED_LOCATION ${Ani3d_FRTPRM_LIBRARY})
target_link_libraries(Ani3d::FRTPRM INTERFACE Ani3d::AFT)
target_include_directories(Ani3d::FRTPRM INTERFACE ${Ani3d_AFT_INCLUDE_SRC})