if (USE_COLLISION)
    set(AV_SOURCES AVSimulator.cpp Configuration.cpp Messurer.cpp Sewer/MinEnergyDeformator.cpp VirtSuture/ContactSurface.cpp)
    set(AV_HEADERS AVSimulator.h Configuration.h Messurer.h Sewer/MinEnergyDeformator.h VirtSuture/ContactSurface.h)
    
    add_library(AorticValve ${AV_SOURCES} ${AV_HEADERS})
    target_link_libraries(AorticValve PRIVATE LeafCollection Core GSL::gsl nlohmann_json::nlohmann_json rpoly_plus_plus GSL::gsl)
    target_include_directories(AorticValve INTERFACE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    install(FILES AVSimulator.h Configuration.h Messurer.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/AVSim/AorticSimulator/)
    install(FILES Sewer/MinEnergyDeformator.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/AVSim/AorticSimulator/Sewer)
    install(FILES VirtSuture/ContactSurface.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/AVSim/AorticSimulator/VirtSuture/)

    set(AS_TARGETS AorticValve PARENT_SCOPE)
endif()