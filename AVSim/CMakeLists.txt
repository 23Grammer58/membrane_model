add_subdirectory(Solvers)
add_subdirectory(Core)
add_subdirectory(DomainCollection)
add_subdirectory(LeafTemplates)
add_subdirectory(AorticSimulator)

set(AVSIM_TARGETS ${SOLVERS_TARGETS} ${CORE_TARGETS} ${DC_TARGETS} ${LT_TARGETS} ${AS_TARGETS} ${AVSIM_TARGETS} PARENT_SCOPE)