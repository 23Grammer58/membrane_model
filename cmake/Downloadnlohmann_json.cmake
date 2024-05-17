include(FetchContent)
include(ExternalProject)

message("-- Download nlohmann_json")
FetchContent_Declare(
        njson_get
        URL https://github.com/nlohmann/json/releases/download/v3.2.0/include.zip
        UPDATE_DISCONNECTED TRUE
        PREFIX "${LIB_DOWNLOAD_PATH}"
        SOURCE_DIR "${LIB_DOWNLOAD_PATH}/nlohmann_json"
)

FetchContent_GetProperties(njson_get)
if(NOT njson_get_POPULATED)
    FetchContent_Populate(njson_get)
endif()

include(FindPackageHandleStandardArgs)
set(NLOHMANN_INCLUDE_DIR ${LIB_DOWNLOAD_PATH}/nlohmann_json)
find_package_handle_standard_args(nlohmann_json DEFAULT_MSG NLOHMANN_INCLUDE_DIR)

if (nlohmann_json_FOUND)
    add_library(nlohmann_json::nlohmann_json INTERFACE IMPORTED GLOBAL)        
    add_dependencies(nlohmann_json::nlohmann_json njson_get)
    target_include_directories(nlohmann_json::nlohmann_json INTERFACE ${NLOHMANN_INCLUDE_DIR})
    set(nlohmann_json_FOUND CACHED true)
endif()

set(nlohmann_json_DOWNLOADED TRUE)
install(DIRECTORY "${LIB_DOWNLOAD_PATH}/nlohmann_json" 
        DESTINATION "${CMAKE_INSTALL_PREFIX}") 

