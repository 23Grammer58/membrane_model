project(hausdorff)

add_executable(hausdorff main.cpp)

target_link_libraries(hausdorff PRIVATE Core_Core)

install(TARGETS hausdorff RUNTIME DESTINATION bin)