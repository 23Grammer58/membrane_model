include(FetchContent)
include(ExternalProject)

message("-- Download RpolyPlusPlus")
FetchContent_Declare(rpoly_get
        GIT_REPOSITORY https://github.com/sweeneychris/RpolyPlusPlus.git
        GIT_TAG        origin/master
        UPDATE_DISCONNECTED TRUE
        PREFIX ${LIB_DOWNLOAD_PATH}
        SOURCE_DIR ${LIB_DOWNLOAD_PATH}/RpolyPlusPlus
        )
FetchContent_GetProperties(rpoly_get)
if(NOT rpoly_get_POPULATED)
    FetchContent_Populate(rpoly_get)
endif()

file(MAKE_DIRECTORY ${LIB_DOWNLOAD_PATH}/RpolyPlusPlus/install)
file(MAKE_DIRECTORY ${LIB_DOWNLOAD_PATH}/RpolyPlusPlus/install/include)
file(MAKE_DIRECTORY ${LIB_DOWNLOAD_PATH}/RpolyPlusPlus/install/lib)
file(GLOB RPOLY_HEADERS "${rpoly_get_SOURCE_DIR}/src/*.h")

enable_language(C)
set(rpoly_args  )
if(eigen3_DOWNLOADED)
        set(rpoly_args ${rpoly_args} "-DEIGEN_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR}")
elseif(Eigen3_DIR) 
        set(rpoly_args ${rpoly_args} "-DEigen_DIR=${Eigen3_DIR}")
endif()                
ExternalProject_Add(rpoly_get
        PREFIX ${LIB_DOWNLOAD_PATH}
        SOURCE_DIR ${rpoly_get_SOURCE_DIR}
        BINARY_DIR ${rpoly_get_SOURCE_DIR}/build
        INSTALL_DIR ${rpoly_get_SOURCE_DIR}/install
        CMAKE_ARGS  
                "-DCMAKE_BUILD_TYPE=Release"
                "-DBUILD_TESTING=OFF"
                "-DCMAKE_C_COMPILER:PATH=${CMAKE_C_COMPILER}"
                "-DBUILD_SHARED_LIBS=ON"
                "-DCMAKE_CXX_FLAGS=-fpermissive"
                ${rpoly_args}
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy ${RPOLY_HEADERS} ${rpoly_get_SOURCE_DIR}/install/include
                        && ${CMAKE_COMMAND} -E copy_directory ${rpoly_get_SOURCE_DIR}/build/lib ${rpoly_get_SOURCE_DIR}/install/lib
        )
unset(rpoly_args)        
if( WIN32 )
    set(LIB_TYPE dylib)
    set(LIB_PREF "")
else( WIN32 )
    set(LIB_TYPE so)
    set(LIB_PREF lib)
endif( WIN32 )

set(rpoly_LIBRARY "${rpoly_get_SOURCE_DIR}/install/lib/${LIB_PREF}rpoly_plus_plus.${LIB_TYPE}")
set(rpoly_INCLUDE_DIR "${rpoly_get_SOURCE_DIR}/install/include")
add_library(rpoly_plus_plus UNKNOWN IMPORTED GLOBAL)
add_dependencies(rpoly_plus_plus rpoly_get)
set_target_properties(rpoly_plus_plus PROPERTIES IMPORTED_LOCATION ${rpoly_LIBRARY})
target_include_directories(rpoly_plus_plus INTERFACE ${rpoly_INCLUDE_DIR})
set(rpoly_DOWNLOADED TRUE)

install(DIRECTORY ${rpoly_get_SOURCE_DIR}/src
        DESTINATION ${rpoly_get_SOURCE_DIR}/install/include
        FILES_MATCHING PATTERN "*.h" )
install(DIRECTORY ${rpoly_get_SOURCE_DIR}/build/lib
        DESTINATION ${rpoly_get_SOURCE_DIR}/install/lib
        FILES_MATCHING PATTERN "*" )