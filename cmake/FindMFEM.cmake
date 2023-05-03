# Find MFEM
# Find the MFEM header files
find_path(MFEM_INCLUDE_DIRS mfem.hpp
  HINTS ${MFEM_DIR} $ENV{MFEM_DIR} ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with MFEM header.")

# Find the MFEM library
find_library(MFEM_LIBRARY mfem
  HINTS ${MFEM_DIR} $ENV{MFEM_DIR} ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The MFEM library.")

# Find config.mk
find_file(MFEM_CONFIG_MK config.mk
  HINTS ${MFEM_DIR} $ENV{MFEM_DIR} ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES share/mfem config
  REQUIRED
  NO_DEFAULT_PATH
  DOC "MFEM's config.mk file.")

# From config.mk, read the values of MFEM_TPLFLAGS and MFEM_EXT_LIBS
file(STRINGS "${MFEM_CONFIG_MK}" MFEM_TPLFLAGS_LINE
  REGEX "^MFEM_TPLFLAGS * = .*$")
string(REGEX REPLACE "^MFEM_TPLFLAGS * = *" ""
  MFEM_TPLFLAGS "${MFEM_TPLFLAGS_LINE}")
file(STRINGS "${MFEM_CONFIG_MK}" MFEM_EXT_LIBS_LINE
  REGEX "^MFEM_EXT_LIBS * = .*$")
string(REGEX REPLACE "^MFEM_EXT_LIBS * = *" ""
  MFEM_EXT_LIBS "${MFEM_EXT_LIBS_LINE}")
# message(STATUS "MFEM_TPLFLAGS: ${MFEM_TPLFLAGS}")
# message(STATUS "MFEM_EXT_LIBS: ${MFEM_EXT_LIBS}")

# set package variables including MFEM_FOUND
find_package_handle_standard_args(MFEM
  REQUIRED_VARS
    MFEM_LIBRARY
    MFEM_INCLUDE_DIRS
  )

# Create target for MFEM
if(MFEM_FOUND)

  if(NOT TARGET XSDK::MFEM)
    add_library(XSDK::MFEM INTERFACE IMPORTED)
  endif()

  set_target_properties(XSDK::MFEM PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${MFEM_INCLUDE_DIRS}"
    INTERFACE_COMPILE_OPTIONS ${MFEM_TPLFLAGS}
    INTERFACE_LINK_LIBRARIES "${MFEM_LIBRARY};${MFEM_EXT_LIBS}")

  # Assuming MPI build of MFEM:
  target_link_libraries(XSDK::MFEM INTERFACE MPI::MPI_C)

   if(ENABLE_CUDA)
     target_link_libraries(XSDK::MFEM INTERFACE CUDA::cudart CUDA::cusparse)
   endif()
   if(ENABLE_HIP)
     target_link_libraries(XSDK::MFEM INTERFACE hip::amdhip64 roc::hipsparse)
   endif()

endif()
