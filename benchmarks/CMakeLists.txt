set(BENCH_TARGETS)
macro(set_target target_name)
    if (TARGET ${target_name})
        set(BENCH_TARGETS ${BENCH_TARGETS} ${target_name})
    endif()
endmacro()    

add_subdirectory(ContactBenchs)
set_target(ContactBenchs)  

add_subdirectory(ElasticBenchs)
set_target(ElasticBenchs) 

add_subdirectory(LinVsNonlin)
set_target(LinVsNonlin) 

add_subdirectory(Sewing)
set_target(Sewing) 
set_target(VirtualSuturer) 
set_target(Suturing) 

add_subdirectory(IdealSuturedLeaf)
set_target(IdealSuturedLeaf) 

add_subdirectory(ConjugatePairs)
set_target(ConjugatePairs) 

add_subdirectory(general)
set_target(model) 

set(AVSIM_TARGETS ${BENCH_TARGETS} ${AVSIM_TARGETS} PARENT_SCOPE)
string(REGEX REPLACE " " ";" BENCH_TARGETS "${BENCH_TARGETS}")

if (UNIX AND USE_MAGNUM_GUI)
    #this flag helps to link transitive dependencies of Magnum libraries
    #see https://sm-stackoverflow.azurefd.net/questions/52018092
    foreach(target ${BENCH_TARGETS})
        target_link_options(${target} PRIVATE "LINKER:--disable-new-dtags")
    endforeach()
endif ()



