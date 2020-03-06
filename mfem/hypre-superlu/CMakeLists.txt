cmake_minimum_required(VERSION 3.8)
project( convdiff
         DESCRIPTION "Convective diffusion using MFEM + SUPERLU + HYPRE"
         LANGUAGES C)

option(MFEM_DIR "absolute/path/to/mfem@develop/install" [.])
option(SUPERLU_DIR "absolute/path/to/superlu-dist@develop/install" [.])
option(HYPRE_DIR "absolute/path/to/hypre@2.18.2+superlu-dist+shared/install" [.])

find_package(MPI REQUIRED)

# Find MFEM
# Find the header files
find_path(MFEM_INCLUDE_DIRS mfem.h
  HINTS ${MFEM_DIR} $ENV{MFEM_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with MFEM header.")
find_path(MFEM_INCLUDE_DIRS mfem.h)

# Find the library
find_library(_LIBRARY MFEM
  HINTS ${MFEM_DIR} $ENV{MFEM_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The MFEM library.")
find_library(MFEM_LIBRARY MFEM)

# Find SUPERLU
# Find the header files
find_path(SUPERLU_INCLUDE_DIRS superlu.h
  HINTS ${SUPERLU_DIR} $ENV{SUPERLU_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with SUPERLU header.")
find_path(SUPERLU_INCLUDE_DIRS superlu.h)

# Find the library
find_library(_LIBRARY SUPERLU
  HINTS ${SUPERLU_DIR} $ENV{SUPERLU_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The SUPERLU library.")
find_library(SUPERLU_LIBRARY SUPERLU)


# Find HYPRE
# Find the header files
find_path(HYPRE_INCLUDE_DIRS HYPRE.h
  HINTS ${HYPRE_DIR} $ENV{HYPRE_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with HYPRE header.")
find_path(HYPRE_INCLUDE_DIRS HYPRE.h)

# Find the library
find_library(HYPRE_LIBRARY HYPRE
  HINTS ${HYPRE_DIR} $ENV{HYPRE_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The HYPRE library.")
find_library(HYPRE_LIBRARY HYPRE)

add_executable(convdiff convdiff.cpp)
target_compile_options(convdiff PRIVATE -g -Wall)
target_link_libraries( convdiff PRIVATE
  ${MFEM_LIBRARY} ${SUPERLU_LIBRARY} ${HYPRE_LIBRARY} ${MPI_C_LIBRARIES})
target_include_directories(convdiff PRIVATE
  ${MFEM_INCLUDE_DIRS} ${SUPERLU_INCLUDE_DIRS} ${HYPRE_INCLUDE_DIRS} ${MPI_C_INCLUDE_PATH})
