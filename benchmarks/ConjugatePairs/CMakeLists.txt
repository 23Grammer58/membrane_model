if (USE_DATADRIVEN_FORCES)
    add_executable(ConjugatePairs CircleTest.cpp CircleTest.h)
    target_link_libraries(ConjugatePairs PRIVATE DomainCollection Core Boost::program_options)
endif()