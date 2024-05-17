#
#
# - Try to find RpolyPlusPlus library
#
# Usage:
#    Control the search through RPOLYDIR or setting environment variable
#    RPOLYDIR or through variable RPOLY_ROOT
#
#    Set target rpoly_plus_plus
#
# #############################################################################

find_library(  rpoly_LIBRARY
        NAMES rpoly_plus_plus
        PATHS $ENV{RPOLYDIR} ${RPOLYDIR} $ENV{RPOLYDIR}/lib ${RPOLYDIR}/lib ${LIB_INSTALL_DIR})

find_path(rpoly_INCLUDE_DIR
        NAMES find_polynomial_roots_jenkins_traub.h
        PATHS $ENV{RPOLYDIR} ${RPOLYDIR} $ENV{RPOLYDIR}/include ${RPOLYDIR}/include $ENV{RPOLYDIR}../src ${RPOLYDIR}../src)

find_package_handle_standard_args(rpoly DEFAULT_MSG rpoly_LIBRARY rpoly_INCLUDE_DIR)
if (rpoly_FOUND)
    mark_as_advanced(rpoly_LIBRARY rpoly_INCLUDE_DIR)

    add_library(rpoly_plus_plus UNKNOWN IMPORTED)
    set_target_properties(rpoly_plus_plus PROPERTIES IMPORTED_LOCATION ${rpoly_LIBRARY})
    target_include_directories(rpoly_plus_plus INTERFACE ${rpoly_INCLUDE_DIR})
endif()