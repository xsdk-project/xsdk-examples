# find the SLATE include path
find_path(SLATE_INCLUDE_DIR slate.hh
  NAMES slate/slate.hh
  HINTS ${SLATE_DIR} $ENV{SLATE_DIR} ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with SLATE header"
)

# find the main SLATE library
find_library(SLATE_LIBRARY_DIR
  NAMES slate
  HINTS ${SLATE_DIR} $ENV{SLATE_DIR} ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
  DOC "The SLATE library."
)

set(SLATE_LIBRARIES "slate;lapackpp;blaspp")

find_package_handle_standard_args(SLATE
  REQUIRED_VARS
    SLATE_LIBRARY_DIR
    SLATE_LIBRARIES
    SLATE_INCLUDE_DIR
  VERSION_VAR
    SLATE_VERSION
)

# Create target for SLATE
if(SLATE_FOUND)

  if(NOT TARGET XSDK::SLATE)
    add_library(XSDK::SLATE UNKNOWN IMPORTED)
  endif()

  set_target_properties(XSDK::SLATE PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SLATE_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${SLATE_LIBRARIES}"
    IMPORTED_LOCATION "${SLATE_LIBRARY_DIR}")

endif()
