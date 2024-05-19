if (USE_COLLISION AND USE_BOOST_LIB)
    add_executable(Sewing SewOnCilinder.cpp)# VirtualSuturer.cpp VirtualSuturer.h)
    target_link_libraries(Sewing PRIVATE LeafCollection AorticValve Core Boost::program_options)

    add_library(VirtualSuturer VirtualSuturer.cpp VirtualSuturer.h TestVirtualSuturer.cpp TestVirtualSuturer.h)
    target_link_libraries(VirtualSuturer PUBLIC LeafCollection AorticValve Core Boost::program_options)
    target_include_directories(VirtualSuturer INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/AVSim/benchmarks/include>)
    install(FILES VirtualSuturer.h TestVirtualSuturer.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/AVSim/benchmarks/include)

    add_executable(Suturing  ForOptimizerOneLeaf.cpp ForOptimizerOneLeaf.h)
    target_link_libraries(Suturing PRIVATE VirtualSuturer)
endif()
