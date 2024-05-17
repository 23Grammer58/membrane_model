#
#
# - Try to find matio library
#
# Usage:
#    Control the search through MATIODIR or setting environment variable
#    MATIODIR or through variable MATIO_ROOT
#
#    Set target MATIO
#
# #############################################################################

find_library(  MATIO_LIBRARY
        NAMES matio
        PATHS $ENV{MATIODIR} ${MATIODIR} $ENV{MATIODIR}/lib ${MATIODIR}/lib ${LIB_INSTALL_DIR})

find_path(MATIO_INCLUDE_DIR
        NAMES matio.h
        PATHS $ENV{MATIODIR} ${MATIODIR} $ENV{MATIODIR}/include ${MATIODIR}/include)

find_package_handle_standard_args(MATIO DEFAULT_MSG MATIO_LIBRARY MATIO_INCLUDE_DIR)
if (MATIO_FOUND)
    mark_as_advanced(MATIO_LIBRARY MATIO_INCLUDE_DIR)

    add_library(MATIO UNKNOWN IMPORTED)
    set_target_properties(MATIO PROPERTIES IMPORTED_LOCATION ${MATIO_LIBRARY})
    target_include_directories(MATIO INTERFACE ${MATIO_INCLUDE_DIR})
endif()
