set(LINSOL_SOURCES LinearSolverInterface.cpp LinearSolverEigen.cpp LinearSolverInmost.cpp LinearSolver.cpp)
set(LINSOL_HEADERS LinearSolverInterface.h LinearSolverEigen.h LinearSolverInmost.h LinearSolver.h)

add_library(LinSol  ${LINSOL_SOURCES} ${LINSOL_HEADERS})
target_link_libraries(LinSol PUBLIC Eigen3::Eigen)
if (USE_INMOST)
    target_link_libraries(LinSol PUBLIC inmost)
endif()

target_include_directories(LinSol INTERFACE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

set(NONLINSOL_SOURCES   NonLinearSolverInterface.cpp NonLinearSolverKinsol.cpp NonLinearSolverPetsc.cpp)
set(NONLINSOL_HEADERS   NonLinearSolverInterface.h NonLinearSolverKinsol.h NonLinearSolverPetsc.h)

add_library(NonLinSol ${NONLINSOL_SOURCES} ${NONLINSOL_HEADERS})
target_link_libraries(NonLinSol PUBLIC SUNDIALS::kinsol LinSol)
target_include_directories(NonLinSol INTERFACE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(FILES ${LINSOL_HEADERS} ${NONLINSOL_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/AVSim/Solvers/)

set(SOLVERS_TARGETS LinSol NonLinSol PARENT_SCOPE)
