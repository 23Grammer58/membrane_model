include(FetchContent)
include(ExternalProject)

message("-- Download GSL")
FetchContent_Declare(gsl_get
        URL https://github.com/ampl/gsl/archive/refs/tags/v2.7.0.zip
        UPDATE_DISCONNECTED TRUE
        PREFIX "${LIB_DOWNLOAD_PATH}"
        SOURCE_DIR "${LIB_DOWNLOAD_PATH}/gsl"
        )
FetchContent_GetProperties(gsl_get)
if(NOT gsl_get_POPULATED)
    FetchContent_Populate(gsl_get)
endif()

file(MAKE_DIRECTORY "${LIB_DOWNLOAD_PATH}/gsl/install")
file(MAKE_DIRECTORY "${LIB_DOWNLOAD_PATH}/gsl/install/include")
file(MAKE_DIRECTORY "${LIB_DOWNLOAD_PATH}/gsl/install/lib")
file(MAKE_DIRECTORY "${LIB_DOWNLOAD_PATH}/gsl/install/bin")

if( WIN32 )
    set(LIB_TYPE lib)
    set(LIB_PREF "")
else( WIN32 )
    set(LIB_TYPE a)
    set(LIB_PREF lib)
endif( WIN32 )
set(gsl_LIBRARY "${gsl_get_SOURCE_DIR}/install/lib/${LIB_PREF}gsl.${LIB_TYPE}")
set(gslcblas_LIBRARY "${gsl_get_SOURCE_DIR}/install/lib/${LIB_PREF}gslcblas.${LIB_TYPE}")

enable_language(C)
ExternalProject_Add(gsl_get
        PREFIX "${LIB_DOWNLOAD_PATH}"
        SOURCE_DIR "${gsl_get_SOURCE_DIR}"
        BINARY_DIR "${gsl_get_SOURCE_DIR}/build"
        INSTALL_DIR "${gsl_get_SOURCE_DIR}/install"
        CMAKE_ARGS  "-DNO_AMPL_BINDINGS=1"
                    "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
                    "-DCMAKE_C_COMPILER:PATH=${CMAKE_C_COMPILER}"
                    "-DCMAKE_INSTALL_PREFIX=${gsl_get_SOURCE_DIR}/install"
        BUILD_BYPRODUCTS "${gsl_LIBRARY}"  
                         "${gslcblas_LIBRARY}"         
        )

add_library(GSL::gslcblas STATIC IMPORTED GLOBAL)
add_dependencies(GSL::gslcblas gsl_get)
set_target_properties(GSL::gslcblas PROPERTIES IMPORTED_LOCATION ${gslcblas_LIBRARY}) 
set_target_properties(GSL::gslcblas PROPERTIES INTERFACE_LINK_LIBRARIES "m")
target_include_directories(GSL::gslcblas INTERFACE "${LIB_DOWNLOAD_PATH}/gsl/install/include")  

add_library(GSL::gsl STATIC IMPORTED GLOBAL)
add_dependencies(GSL::gsl gsl_get)
set_target_properties(GSL::gsl PROPERTIES IMPORTED_LOCATION ${gsl_LIBRARY}) 
set_target_properties(GSL::gsl PROPERTIES INTERFACE_LINK_LIBRARIES "m;GSL::gslcblas")
target_include_directories(GSL::gsl INTERFACE "${LIB_DOWNLOAD_PATH}/gsl/install/include")

set(gsl_DOWNLOADED TRUE)
install(DIRECTORY "${LIB_DOWNLOAD_PATH}/gsl/install/" 
        DESTINATION "${CMAKE_INSTALL_PREFIX}") 