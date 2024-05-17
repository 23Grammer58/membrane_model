#
#
# - Try to find mesh-generator tools of library Ani3d
#
# Usage:
#    Control the search through ANI3DDIR or setting environment variable
#    ANI3DDIR or through variable Ani3d_ROOT
#
#    Set targets Ani3d::AFT Ani3d::MBA Ani3d::FRTPRM
#
# #############################################################################

include(FindPackageHandleStandardArgs)

find_library(  Ani3d_AFT_LIBRARY
        NAMES aft3D aft3D-3.0
        PATHS $ENV{ANI3DDIR} ${ANI3DDIR} $ENV{ANI3DDIR}/lib ${ANI3DDIR}/lib ${LIB_INSTALL_DIR})
find_library(  Ani3d_MBA_LIBRARY
        NAMES mba3D mba3D-3.0
        PATHS $ENV{ANI3DDIR} ${ANI3DDIR} $ENV{ANI3DDIR}/lib ${ANI3DDIR}/lib ${LIB_INSTALL_DIR})
find_library(  Ani3d_FRTPRM_LIBRARY
        NAMES frtprm frtprm-3.0
        PATHS $ENV{ANI3DDIR} ${ANI3DDIR} $ENV{ANI3DDIR}/lib ${ANI3DDIR}/lib ${LIB_INSTALL_DIR})
find_path(Ani3d_AFT_INCLUDE_DIR
        NAMES libaft.h
        PATHS $ENV{ANI3DDIR} ${ANI3DDIR} $ENV{ANI3DDIR}/include ${ANI3DDIR}/include)
if (Ani3d_AFT_INCLUDE_DIR)
    set(Ani3d_AFT_INCLUDE_SRC ${Ani3d_AFT_INCLUDE_DIR}/../src)
endif()
find_package(Fortran)
find_package(BLAS)
find_package(LAPACK)

find_package_handle_standard_args(Ani3d DEFAULT_MSG
        Ani3d_AFT_LIBRARY Ani3d_MBA_LIBRARY Ani3d_FRTPRM_LIBRARY Ani3d_AFT_INCLUDE_DIR Ani3d_AFT_INCLUDE_SRC)

if (Ani3d_FOUND)
    mark_as_advanced(Ani3d_AFT_LIBRARY
                     Ani3d_MBA_LIBRARY
                     Ani3d_FRTPRM_LIBRARY
                     Ani3d_AFT_INCLUDE_DIR
                     Ani3d_AFT_INCLUDE_SRC)
endif()

if (Ani3d_FOUND)
    add_library(Ani3d::AFT UNKNOWN IMPORTED)
    set_target_properties(Ani3d::AFT PROPERTIES IMPORTED_LOCATION ${Ani3d_AFT_LIBRARY})
    target_include_directories(Ani3d::AFT INTERFACE ${Ani3d_AFT_INCLUDE_DIR})
    set(Ani3d_DEPENDENCY_LIBS "${LAPACK_LIBRARIES};${BLAS_LIBRARIES};${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES};m")
    target_link_libraries(Ani3d::AFT INTERFACE "\$<LINK_ONLY:${Ani3d_DEPENDENCY_LIBS}>")

    add_library(Ani3d::MBA UNKNOWN IMPORTED)
    set_target_properties(Ani3d::MBA PROPERTIES IMPORTED_LOCATION ${Ani3d_MBA_LIBRARY})
    target_link_libraries(Ani3d::MBA INTERFACE Ani3d::AFT)

    add_library(Ani3d::FRTPRM UNKNOWN IMPORTED)
    set_target_properties(Ani3d::FRTPRM PROPERTIES IMPORTED_LOCATION ${Ani3d_FRTPRM_LIBRARY})
    target_link_libraries(Ani3d::FRTPRM INTERFACE Ani3d::AFT)
    target_include_directories(Ani3d::FRTPRM INTERFACE ${Ani3d_AFT_INCLUDE_SRC})
endif()
