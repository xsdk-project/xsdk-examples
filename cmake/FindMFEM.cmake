# MFEM requires metis.
find_package(METIS REQUIRED)

# Find MFEM
# Find the MFEM header files
find_path(MFEM_INCLUDE_DIRS mfem.hpp
  HINTS ${MFEM_DIR} $ENV{MFEM_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with MFEM header.")

# Find the MFEM library
find_library(MFEM_LIBRARY mfem
  HINTS ${MFEM_DIR} $ENV{MFEM_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The MFEM library.")

# set package variables including MFEM_FOUND
find_package_handle_standard_args(MFEM
  REQUIRED_VARS
    MFEM_LIBRARY
    MFEM_INCLUDE_DIRS
  )

# Create target for MFEM
if(MFEM_FOUND)

  if(NOT TARGET XSDK::MFEM)
    add_library(XSDK::MFEM UNKNOWN IMPORTED)
  endif()

  set_target_properties(XSDK::MFEM PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${MFEM_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${METIS_LIBRARY}"
    IMPORTED_LOCATION "${MFEM_LIBRARY}")

endif()
