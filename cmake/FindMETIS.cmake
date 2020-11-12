
find_path(METIS_INCLUDE_DIRS metis.h
    HINTS ${METIS_DIR} $ENV{METIS_DIR}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
    DOC "The directory with the Metis header file.")

find_library(METIS_LIBRARY metis
    HINTS ${METIS_DIR} $ENV{METIS_DIR}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
    DOC "The Metis library.")

find_library(ZLIB_LIBRARY z
    HINTS ${ZLIB_LIBRARY_DIR} $ENV{ZLIB_LIBRARY_DIR}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
    DOC "The zlib library.")

# set package variables including METIS_FOUND
find_package_handle_standard_args(METIS
  REQUIRED_VARS
    METIS_LIBRARY
    METIS_INCLUDE_DIRS
  )

# Create target for METIS
if(METIS_FOUND)

  if(NOT TARGET METIS)
    add_library(METIS UNKNOWN IMPORTED)
  endif()

  set_target_properties(METIS PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${ZLIB_LIBRARY}"
    IMPORTED_LOCATION "${METIS_LIBRARY}")

endif()
