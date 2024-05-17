set(CORE_HEADERS    AVMesh.h AVMeshConnectivityHelpers.h TriangularMeshHelpers.h
                    Object3D.h ForceInterface.h ForceAppliersInterface.h CollisionManagerInterface.h CollisionManagerInterface.inl RendererInterface.h
                    World.h NSWorldWrapper.h NonLinearSolverCustom.h SignedDistanceField.h IO/Obj3dVTKSaver.h)
set(CORE_SOURCES    TriangularMeshHelpers.cpp Object3D.cpp ForceInterface.cpp World.cpp NSWorldWrapper.cpp NonLinearSolverCustom.cpp CollisionManagerInterface.cpp IO/Obj3dVTKSaver.cpp)
prefixed_install("${CORE_HEADERS}" "AVSim/Core")

add_subdirectory(Renderers)
list(TRANSFORM RENDERERS_SOURCES PREPEND Renderers/)
list(TRANSFORM RENDERERS_HEADERS PREPEND Renderers/)
prefixed_install("${RENDERERS_HEADERS}" "AVSim/Core")

add_subdirectory(ForceAppliers)
list(TRANSFORM FAPL_SOURCES PREPEND ForceAppliers/)
list(TRANSFORM FAPL_HEADERS PREPEND ForceAppliers/)
prefixed_install("${FAPL_HEADERS}" "AVSim/Core")

add_subdirectory(Collision)
list(TRANSFORM COLLISSION_SOURCES PREPEND Collision/)
list(TRANSFORM COLLISSION_HEADERS PREPEND Collision/)

add_subdirectory(Forces)
list(TRANSFORM FORCE_SOURCES PREPEND Forces/)
list(TRANSFORM FORCE_HEADERS PREPEND Forces/)
list(TRANSFORM DATA_DRIVEN_SOURCES PREPEND Forces/)
list(TRANSFORM DATA_DRIVEN_HEADERS PREPEND Forces/)
list(TRANSFORM CONTACT_FORCE_SOURCES PREPEND Forces/)
list(TRANSFORM CONTACT_FORCE_HEADERS PREPEND Forces/)
prefixed_install("${FORCE_HEADERS}" "AVSim/Core")

add_subdirectory(MeshGen)
list(TRANSFORM MGEN_SOURCES PREPEND MeshGen/)
list(TRANSFORM MGEN_HEADERS PREPEND MeshGen/)
prefixed_install("${MGEN_HEADERS}" "AVSim/Core")

add_subdirectory(IO)
list(TRANSFORM IO_SOURCES PREPEND IO/)
list(TRANSFORM IO_HEADERS PREPEND IO/)
prefixed_install("${IO_HEADERS}" "AVSim/Core")

add_library(Core_Core ${CORE_SOURCES} ${CORE_HEADERS})
target_link_libraries(Core_Core NonLinSol CGAL::CGAL Ani3d::FRTPRM)
target_compile_definitions(Core_Core PRIVATE -DUSE_ANI3D)

add_library(Core_Renderers ${RENDERERS_SOURCES} ${RENDERERS_HEADERS})
target_link_libraries(Core_Renderers PUBLIC Core_Core)
if (USE_MAGNUM_GUI)
    target_link_libraries(Core_Renderers PUBLIC Magnum::Trade MagnumIntegration::ImGui Magnum::Primitives Magnum::Shaders Magnum::MeshTools Magnum::Application Magnum::GL Magnum::SceneGraph Magnum::Magnum Corrade::Main)
endif()

add_library(Core_ForceAppliers ${FAPL_SOURCES} ${FAPL_HEADERS})
target_link_libraries(Core_ForceAppliers PUBLIC Core_Core)

add_library(Core_MeshGen ${MGEN_SOURCES} ${MGEN_HEADERS})
target_link_libraries(Core_MeshGen PUBLIC Core_Core Ani3d::FRTPRM)
target_compile_definitions(Core_MeshGen PRIVATE -DUSE_ANI3D)

add_library(Core_IO ${IO_SOURCES} ${IO_HEADERS})
target_link_libraries(Core_IO PUBLIC Core_Core casadi)
if (USE_INMOST)
    target_link_libraries(Core_IO PRIVATE inmost)
endif()    

add_library(Core_Forces ${FORCE_SOURCES} ${FORCE_HEADERS})
target_link_libraries(Core_Forces PUBLIC Core_Core casadi)

add_library(Core INTERFACE)
target_link_libraries(Core INTERFACE Core_Forces Core_IO Core_MeshGen Core_ForceAppliers Core_Renderers Core_Core)
set(CORE_TARGETS Core Core_Forces Core_IO Core_MeshGen Core_ForceAppliers Core_Renderers Core_Core)

if (USE_COLLISION)
    prefixed_install("${COLLISSION_HEADERS}" "AVSim/Core")
    add_library(Core_Collision ${COLLISSION_SOURCES} ${COLLISSION_HEADERS})
    target_link_libraries(Core_Collision PUBLIC Core_Core Bullet_Math Bullet_Collision)

    prefixed_install("${CONTACT_FORCE_HEADERS}" "AVSim/Core")
    add_library(Core_ContactForce ${CONTACT_FORCE_SOURCES} ${CONTACT_FORCE_HEADERS})
    target_link_libraries(Core_ContactForce PUBLIC Core_Collision)

    target_link_libraries(Core INTERFACE Core_ContactForce Core_Collision)
    set(CORE_TARGETS Core_ContactForce Core_Collision ${CORE_TARGETS})
endif()

if (USE_DATADRIVEN_FORCES)
    prefixed_install("${DATA_DRIVEN_HEADERS}" "AVSim/Core")
    add_library(Core_DataDriven ${DATA_DRIVEN_SOURCES} ${DATA_DRIVEN_HEADERS})
    target_link_libraries(Core_DataDriven PUBLIC Core_Forces Core_ForceAppliers Boost::serialization Boost::filesystem GSL::gsl)

    target_link_libraries(Core INTERFACE Core_DataDriven)
    set(CORE_TARGETS Core_DataDriven ${CORE_TARGETS}) 
endif()

target_include_directories(Core INTERFACE  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set(CORE_TARGETS ${CORE_TARGETS} PARENT_SCOPE)                