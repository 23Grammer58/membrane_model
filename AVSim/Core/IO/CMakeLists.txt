set(IO_SOURCES MetaTriFormat.cpp VTKFormat.cpp MetaTriMesh.cpp HyperElasticFormat.cpp)
set(IO_HEADERS MetaTriFormat.h VTKFormat.h MetaTriMesh.h MetaTriMesh.inl HyperElasticFormat.h HyperElasticFormat.inl)
if (USE_INMOST)
    set(IO_SOURCES ${IO_SOURCES} MeshlabPickPoint.cpp)
    set(IO_HEADERS ${IO_HEADERS} MeshlabPickPoint.h)
endif()
set(IO_SOURCES ${IO_SOURCES} PARENT_SCOPE)
set(IO_HEADERS ${IO_HEADERS} PARENT_SCOPE)