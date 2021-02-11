cmake_minimum_required(VERSION 3.8)
project(mfem-gko LANGUAGES CXX)

find_package(OpenMP)

find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

option(ENABLE_MFEM "Enable MFEM" ON)

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

find_library(GINKGO_LIBRARY ginkgo
    HINTS ${GINKGO_LIBRARY_DIR} ${GINKGO_DIR})

find_library(GINKGO_OMP_LIBRARY ginkgo_omp
    HINTS ${GINKGO_LIBRARY_DIR} ${GINKGO_DIR})

find_library(GINKGO_CUDA_LIBRARY ginkgo_cuda
    HINTS ${GINKGO_LIBRARY_DIR} ${GINKGO_DIR})

add_executable(mfem_ex1_gko mfem_ex1_gko.cpp)
target_sources(mfem_ex1_gko PRIVATE mfem_ex1_gko.cpp mfem_wrapper.cpp)

target_link_libraries(mfem_ex1_gko PRIVATE ${MFEM_LIBRARY} XSDK::SUNDIALS XSDK::HYPRE ${GINKGO_LIBRARY} ${GINKGO_OMP_LIBRARY} ${GINKGO_CUDA_LIBRARY} PETSC::ALL ${METIS_LIBRARY} ZLIB::ZLIB)

target_include_directories(mfem_ex1_gko PRIVATE ${MFEM_INCLUDE_DIR} ${GINKGO_INCLUDE_DIR})
