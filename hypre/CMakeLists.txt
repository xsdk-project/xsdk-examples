cmake_minimum_required(VERSION 3.8)
project( hypre-superludist
         DESCRIPTION "HYPRE + SuperLU_DIST Example"
         LANGUAGES C)

set(HYPRE_DIR "/path/to/hypre@2.18.2+superlu-dist+shared")
set(DSUPERLU_DIR "/path/to/superlu-dist@6.1.1")

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
find_path(DSUPERLU_INCLUDE_DIRS superlu_defs.h
  HINTS ${DSUPERLU_DIR} $ENV{DSUPERLU_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with SuperLU_DIST defs header.")
find_path(DSUPERLU_INCLUDE_DIRS superlu_defs.h)

# Find the library
find_library(DSUPERLU_LIBRARY superlu_dist
  HINTS ${DSUPERLU_DIR} $ENV{DSUPERLU_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The SuperLU_DIST library.")
find_library(DSUPERLU_LIBRARY superlu_dist)

add_executable(ij_laplacian ij_laplacian.c)
target_compile_options(ij_laplacian PRIVATE -g -Wall)
target_link_libraries(ij_laplacian PRIVATE
  ${HYPRE_LIBRARY} ${DSUPERLU_LIBRARY} ${MPI_C_LIBRARIES})
target_include_directories(ij_laplacian PRIVATE
  ${HYPRE_INCLUDE_DIRS} ${DSUPERLU_INCLUDE_DIRS} ${MPI_C_INCLUDE_PATH})

# Copy glvis helper script and directory
file(COPY vis DESTINATION .)
