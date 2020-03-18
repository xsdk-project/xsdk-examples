cmake_minimum_required(VERSION 3.8)
project( hypre-superludist
         DESCRIPTION "HYPRE + SuperLU_DIST Example"
         LANGUAGES C)

option(HYPRE_DIR "absolute/path/to/hypre@2.18.2+superlu-dist+shared/install" [.])
option(DSUPERLU_DIR "absolute/path/to/superlu-dist@6.1.1/install" [.])

find_package(MPI REQUIRED)

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

# Link libraries for Unix systems
if (UNIX)
  list (APPEND HYPRE_LIBRARY m)
endif (UNIX)

# Find SuperLU_DIST
find_path(SUPERLUDIST_INCLUDE_DIR superlu_defs.h
  HINTS ${SUPERLUDIST_DIR} $ENV{SUPERLUDIST_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with SuperLU_DIST defs header.")
find_path(DSUPERLU_INCLUDE_DIRS superlu_defs.h)

# Find the library
find_library(SUPERLUDIST_LIBRARY superlu_dist
  HINTS ${SUPERLUDIST_DIR} $ENV{SUPERLUDIST_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The SuperLU_DIST library.")
find_library(SUPERLUDIST_LIBRARY superlu_dist)

add_executable(ij_laplacian ij_laplacian.c)
target_compile_options(ij_laplacian PRIVATE -g -Wall)
target_link_libraries(ij_laplacian PRIVATE
  ${HYPRE_LIBRARY} ${SUPERLUDIST_LIBRARY} ${MPI_C_LIBRARIES})
target_include_directories(ij_laplacian PRIVATE
  ${HYPRE_INCLUDE_DIRS} ${SUPERLUDIST_INCLUDE_DIRS} ${MPI_C_INCLUDE_PATH})

# Copy glvis helper script and directory
file(COPY vis DESTINATION .)
