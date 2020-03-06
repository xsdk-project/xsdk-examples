cmake_minimum_required(VERSION 3.8)
project( obstacle
         DESCRIPTION "Obstacle problem using MFEM + PETSC/Tao"
         LANGUAGES C)

option(MFEM_DIR "absolute/path/to/mfem@develop/install" [.])
option(PETSC_DIR "absolute/path/to/petsc@develop/install" [.])

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

# Find PETSC
# Find the header files
find_path(PETSC_INCLUDE_DIRS petsc.h
  HINTS ${PETSC_DIR} $ENV{PETSC_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with PETSC header.")
find_path(PETSC_INCLUDE_DIRS petsc.h)

# Find the library
find_library(_LIBRARY PETSC
  HINTS ${PETSC_DIR} $ENV{PETSC_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The PETSC library.")
find_library(PETSC_LIBRARY PETSC)

add_executable(obstacle obstacle.cpp)
target_compile_options(obstacle PRIVATE -g -Wall)
target_link_libraries( obstacle PRIVATE
  ${MFEM_LIBRARY} ${PETSC_LIBRARY} ${HYPRE_LIBRARY} ${MPI_C_LIBRARIES})
target_include_directories(obstacle PRIVATE
  ${MFEM_INCLUDE_DIRS} ${PETSC_INCLUDE_DIRS} ${HYPRE_INCLUDE_DIRS} ${MPI_C_INCLUDE_PATH})
