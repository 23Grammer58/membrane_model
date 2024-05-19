if (USE_BOOST_LIB)
    add_executable(ConjugatePairs CircleTest.cpp CircleTest.h)
    target_link_libraries(ConjugatePairs PRIVATE DomainCollection Core Boost::program_options)
endif()