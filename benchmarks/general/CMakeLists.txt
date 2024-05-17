add_executable(model main.cpp model.h model.cpp examples.h examples.cpp)
target_link_libraries(model PRIVATE DomainCollection Core_Forces Core_IO Core_ForceAppliers Core_Renderers)
target_include_directories(model PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})