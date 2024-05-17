if (TARGET AorticValve)
    add_executable(IdealSuturedLeaf runBench.cpp VanLoonBench.cpp VanLoonBench.h)
    target_link_libraries(IdealSuturedLeaf PRIVATE LeafCollection AorticValve Core)
endif()