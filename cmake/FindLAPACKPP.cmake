# find the LAPACK++ include path
find_path(LAPACKPP_INCLUDE_DIR lapack.hh
  NAMES lapack.hh
  HINTS ${LAPACKPP_DIR} $ENV{LAPACKPP_DIR} ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with LAPACKPP header"
)

# find the main LAPACK++ library
find_library(LAPACKPP_LIBRARIES
  NAMES lapackpp
  HINTS ${LAPACKPP_DIR} $ENV{LAPACKPP_DIR} ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
  DOC "The LAPACK++ library."
)

find_package_handle_standard_args(LAPACKPP
  REQUIRED_VARS
  LAPACKPP_LIBRARIES
  LAPACKPP_INCLUDE_DIR
  VERSION_VAR
  LAPACKPP_VERSION
)

# Create target for LAPACK++
if(LAPACKPP_FOUND)

  if(NOT TARGET XSDK::LAPACKPP)
    add_library(XSDK::LAPACKPP INTERFACE IMPORTED)
  endif()

  message(STATUS "Created XSDK::LAPACKPP target")
  message(STATUS "   INTERFACE_INCLUDE_DIRECTORIES: ${LAPACKPP_INCLUDE_DIR}")
  message(STATUS "   INTERFACE_LINK_LIBRARIES: ${LAPACKPP_LIBRARIES}")

  set_target_properties(XSDK::LAPACKPP PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${LAPACKPP_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${LAPACKPP_LIBRARIES}")

endif()