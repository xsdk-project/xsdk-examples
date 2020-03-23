cmake_minimum_required(VERSION 3.8)
project( convdiff
         DESCRIPTION "Convective diffusion using MFEM + SUPERLU + HYPRE"
         LANGUAGES C)

option(SUPERLUDIST_DIR "absolute/path/to/superlu-dist@develop/install" [.])

find_package(MPI REQUIRED)

option(ENABLE_MFEM "Enable MFEM" ON)

find_library(ZLIB_LIBRARY z
    HINTS ${ZLIB_LIBRARY_DIR} $ENV{METIS_LIBRARY_DIR}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
    DOC "The zlib library.")
find_library(ZLIB_LIBRARY z)

find_library(METIS_LIBRARY metis
    HINTS ${METIS_DIR} $ENV{METIS_DIR}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
    DOC "The Metis library.")
find_library(METIS_LIBRARY metis)

# Find MFEM
# Find the MFEM header files
find_path(MFEM_INCLUDE_DIRS mfem.hpp
  HINTS ${MFEM_DIR} $ENV{MFEM_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with MFEM header.")
find_path(MFEM_INCLUDE_DIRS mfem.hpp)

# Find the MFEM library
find_library(MFEM_LIBRARY mfem
  HINTS ${MFEM_DIR} $ENV{MFEM_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The MFEM library.")
find_library(MFEM_LIBRARY MFEM)

# Find SUPERLUDIST
# Find the header files
find_path(SUPERLUDIST_INCLUDE_DIRS superlu_defs.h
    HINTS ${SUPERLUDIST_DIR} $ENV{SUPERLUDIST_DIR}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
    DOC "Directory with SUPERLUDIST header.")
find_path(SUPERLUDIST_INCLUDE_DIRS superlu_defs.h)

# Find the library
find_library(SUPERLUDIST_LIBRARY superlu_dist
    HINTS ${SUPERLUDIST_DIR} $ENV{SUPERLUDIST_DIR}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
    DOC "The SUPERLUDIST library.")
find_library(SUPERLUDIST_LIBRARY SUPERLU)


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
target_link_libraries( convdiff PRIVATE
    ${MFEM_LIBRARY} ${SUPERLU_LIBRARY} ${HYPRE_LIBRARY} ${MPI_C_LIBRARIES})
target_include_directories(convdiff PRIVATE
    ${MFEM_INCLUDE_DIRS} ${SUPERLUDIST_INCLUDE_DIR} ${HYPRE_INCLUDE_DIRS} ${MPI_C_INCLUDE_PATH})
