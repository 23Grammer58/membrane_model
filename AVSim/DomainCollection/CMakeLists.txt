add_library(DomainCollection
                            AniMeshCollection.h
                            EightSpherePart.cpp
                            Annulus.cpp
                            CilinderPart.cpp
                            Circle.cpp
                            CutAnnulus.cpp
                            FragmentalRectangle.cpp
                            HalfCircle.cpp
                            HalfSphere.cpp
                            HoleyRectangle.cpp
                            IdealHexagon.cpp
                            RectangleWithOxLine.cpp
                            RectangleWithOxOyLine.cpp
                            CilEllCoaptTemplate.cpp
                            EllipsoidPart.cpp
                            CrossPiece.cpp
                            )

target_link_libraries(DomainCollection PRIVATE Core_MeshGen)
target_include_directories(DomainCollection INTERFACE  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(FILES AniMeshCollection.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/AVSim/DomainCollection/)

set(DC_TARGETS DomainCollection PARENT_SCOPE)