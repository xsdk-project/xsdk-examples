# find the PLASMA include path
find_path(PLASMA_INCLUDE_DIR plasma.h
  NAMES plasma.h
  HINTS ${PLASMA_DIR} $ENV{PLASMA_DIR} ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with PLASMA header"
)

# find the main PLASMA library
find_library(PLASMA_LIBRARY
  NAMES plasma
  HINTS ${PLASMA_DIR} $ENV{PLASMA_DIR} ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
  DOC "The PLASMA library."
)

find_library(PLASMA_CORE_BLAS_LIBRARY
  NAMES plasma_core_blas
  HINTS ${PLASMA_DIR} $ENV{PLASMA_DIR} ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
  DOC "The PLASMA core blas library."
)
  
set(PLASMA_LIBRARIES "${PLASMA_LIBRARY};${PLASMA_CORE_BLAS_LIBRARY}")

find_package_handle_standard_args(PLASMA
  REQUIRED_VARS
    PLASMA_LIBRARIES
    PLASMA_INCLUDE_DIR
  VERSION_VAR
    PLASMA_VERSION
)

# Create target for PLASMA
if(PLASMA_FOUND)
  if(NOT TARGET XSDK::PLASMA)
    add_library(XSDK::PLASMA INTERFACE IMPORTED)
  endif()
  
  message(STATUS "Created XSDK::PLASMA target")
  message(STATUS "   INTERFACE_INCLUDE_DIRECTORIES: ${PLASMA_INCLUDE_DIR}")
  message(STATUS "   INTERFACE_LINK_LIBRARIES: ${PLASMA_LIBRARIES}")

  set_target_properties(XSDK::PLASMA PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${PLASMA_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${PLASMA_LIBRARIES}")
endif()
