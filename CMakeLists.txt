cmake_minimum_required(VERSION 3.12)
project(xsdk-examples
        DESCRIPTION "xSDK Examples"
        LANGUAGES CXX C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(FindPackageHandleStandardArgs)

# build options
option(ENABLE_CUDA "Enable CUDA" OFF)
option(CMAKE_CUDA_ARCHITECTURES "CUDA architecture(s) to target" "70")

option(ENABLE_HYPRE "Enable hypre" ON)
set(HYPRE_DIR "${HYPRE_DIR}" CACHE PATH "Path to hypre installation directory")

option(ENABLE_MFEM "Enable MFEM" ON)
set(MFEM_DIR "${MFEM_DIR}" CACHE PATH "Path to MFEM installation directory")

option(ENABLE_PETSC "Enable PETSc" ON)
set(PETSC_DIR "${PETSC_DIR}" CACHE PATH "Path to PETSc installation directory")

option(ENABLE_SUNDIALS "Enable SUNDIALS" ON)
set(SUNDIALS_DIR "${SUNDIALS_DIR}" CACHE PATH "Path to SUNDIALS installation directory")

option(ENABLE_SUPERLU "Enable SuperLU" ON)
set(SUPERLUDIST_DIR "${SUPERLUDIST__DIR}" CACHE PATH "Path to SuperLU_DIST installation directory")

option(ENABLE_TRILINOS "Enable TRILINOS" ON)
set(TRILINOS_DIR "${Trilinos_DIR}" CACHE PATH "Path to Trilinos installation directory")

set(METIS_DIR "${METIS_DIR}" CACHE PATH "Path to Metis installation directory")

if(ENABLE_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
endif()

# check for hypre
if(ENABLE_HYPRE)
  find_package(HYPRE REQUIRED)
endif()

# check for PETSc
if(ENABLE_MFEM)
  find_package(ZLIB REQUIRED)
  find_package(MFEM REQUIRED)
  find_package(Ginkgo REQUIRED)
endif()

# check for PETSC
if(ENABLE_PETSC)
  find_package(PETSC REQUIRED)
endif()

# check for SUNDIALS
if(ENABLE_SUNDIALS)
  find_package(SUNDIALS REQUIRED)
endif()

# check for SuperLU DIST
if(ENABLE_SUPERLU)
  find_package(SUPERLUDIST REQUIRED)
endif()

# check for math
find_library(MATH_LIBRARY NAMES m)

# check for metis
find_package(METIS)

# check for MPI
find_package(MPI REQUIRED)

# example subdirectories
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
