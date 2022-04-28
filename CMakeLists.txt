cmake_minimum_required(VERSION 3.18)
project(
  xsdk-examples
  DESCRIPTION "xSDK Examples"
  LANGUAGES CXX C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(FindPackageHandleStandardArgs)

# build options
option(ENABLE_AMREX "Enable AMReX" ON)
set(AMREX_DIR
  "${AMREX_DIR}"
  CACHE PATH "Path to AMReX installation directory")

option(ENABLE_CUDA "Enable CUDA" OFF)

option(ENABLE_HYPRE "Enable hypre" ON)
set(HYPRE_DIR
    "${HYPRE_DIR}"
    CACHE PATH "Path to hypre installation directory")

option(ENABLE_MFEM "Enable MFEM" ON)
set(MFEM_DIR
    "${MFEM_DIR}"
    CACHE PATH "Path to MFEM installation directory")

option(ENABLE_MAGMA "Enable MAGMA" OFF)
set(MAGMA_DIR
    "${MAGMA_DIR}"
    CACHE PATH "Path to MAGMA installation directory")

option(ENABLE_PETSC "Enable PETSc" ON)
set(PETSc_DIR
    "${PETSc_DIR}"
    CACHE PATH "Path to PETSc installation directory")

option(ENABLE_SUNDIALS "Enable SUNDIALS" ON)
set(SUNDIALS_DIR
    "${SUNDIALS_DIR}"
    CACHE PATH "Path to SUNDIALS installation directory")

option(ENABLE_SUPERLU "Enable SuperLU" ON)
set(SUPERLUDIST_DIR
    "${SUPERLUDIST_DIR}"
    CACHE PATH "Path to SuperLU_DIST installation directory")

option(ENABLE_STRUMPACK "Enable STRUMPACK" OFF)
set(STRUMPACK_DIR
    "${STRUMPACK_DIR}"
    CACHE PATH "Path to STRUMPACK installation directory")

option(ENABLE_TRILINOS "Enable TRILINOS" OFF)
set(TRILINOS_DIR
    "${Trilinos_DIR}"
    CACHE PATH "Path to Trilinos installation directory")

set(METIS_DIR
    "${METIS_DIR}"
    CACHE PATH "Path to Metis installation directory")

# compiler options
if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

# check for AMReX
if(ENABLE_AMREX)
  find_package(AMReX REQUIRED)
endif()

# setup CUDA
if(ENABLE_CUDA)
  enable_language(CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 70)
  endif()
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  find_package(CUDAToolkit REQUIRED)
endif()

# check for hypre
if(ENABLE_HYPRE)
  find_package(HYPRE REQUIRED)
endif()

# check for MFEM
if(ENABLE_MFEM)
  find_package(ZLIB REQUIRED)
  find_package(MFEM REQUIRED)
  find_package(Ginkgo REQUIRED)
endif()

# check for MAGMA
if(ENABLE_MAGMA)
  find_package(MAGMA REQUIRED)
endif()

# check for PETSC
if(ENABLE_PETSC)
  find_package(PETSc REQUIRED)
endif()

# check for SUNDIALS
if(ENABLE_SUNDIALS)
  find_package(SUNDIALS REQUIRED)
endif()

# check for SuperLU DIST
if(ENABLE_SUPERLU)
  find_package(SUPERLUDIST REQUIRED)
endif()

# check for STRUMPACK
if(ENABLE_STRUMPACK)
  find_package(STRUMPACK REQUIRED)
endif()

# check for math
find_library(MATH_LIBRARY NAMES m)

# check for metis
find_package(METIS)

# check for MPI
find_package(MPI REQUIRED)

# example subdirectories
if(ENABLE_AMREX)
  add_subdirectory(amrex)
endif()
if(ENABLE_HYPRE)
  add_subdirectory(hypre)
endif()
if(ENABLE_MFEM)
  add_subdirectory(mfem)
endif()
if(ENABLE_PETSC)
  add_subdirectory(petsc)
endif()
if(ENABLE_SUNDIALS)
  add_subdirectory(sundials)
endif()
if(ENABLE_TRILINOS)
  add_subdirectory(trilinos)
endif()
if(ENABLE_STRUMPACK)
  add_subdirectory(strumpack)
endif()
