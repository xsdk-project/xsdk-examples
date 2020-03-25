cmake_minimum_required(VERSION 3.8)
project( obstacle
         DESCRIPTION "Obstacle problem using MFEM + PETSC/Tao"
         LANGUAGES CXX)

find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

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
find_library(MFEM_LIBRARY mfem)

find_package(PETSC REQUIRED)

find_package(SUNDIALS REQUIRED)

# Find PETSC
# Find the header files
find_path(PETSC_INCLUDE_DIRS petsc.h
  HINTS ${PETSC_DIR} $ENV{PETSC_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with PETSC header.")
find_path(PETSC_INCLUDE_DIRS petsc.h)

# Find the library
find_library(PETSC_LIBRARY petsc
  HINTS ${PETSC_DIR} $ENV{PETSC_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The PETSC library.")
find_library(PETSC_LIBRARY petsc)

add_executable(obstacle obstacle.cpp)
target_link_libraries( obstacle PRIVATE
    ${MFEM_LIBRARY} SUNDIALS::nvecserial SUNDIALS::nvecparhyp PETSC::ALL ${HYPRE_LIBRARY} ${METIS_LIBRARY} ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})
target_include_directories(obstacle PRIVATE
    ${MFEM_INCLUDE_DIRS} ${PETSC_DIR}/include ${HYPRE_INCLUDE_DIRS} ${MPI_C_INCLUDE_PATH})

install(TARGETS obstacle RUNTIME DESTINATION bin)
