# ---------------------------------------------------------------
# Programmer: Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------

cmake_minimum_required(VERSION 3.10)

# Prohibit in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source build prohibited.")
endif()

# set the project name and default languages
project(SUNDIALS-XSDK-EXAMPLES CXX C)

# location of additional CMake modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)

# build options
set(SUNDIALS_DIR "${SUNDIALS_DIR}" CACHE PATH "Path to SUNDIALS installation directory")
option(ENABLE_PETSC "Enable PETSc" ON)
set(PETSC_DIR "${PETSC_DIR}" CACHE PATH "Path to PETSc installation directory")
option(ENABLE_SUPERLU "Enable SuperLU" ON)
set(SUPERLUDIST_LIBRARY_DIR "${SUPERLUDIST_LIBRARY_DIR}" CACHE PATH "Path to SuperLU_DIST libraries")
set(SUPERLUDIST_INCLUDE_DIR "${SUPERLUDIST_INCLUDE_DIR}" CACHE PATH "Path to SuperLU_DIST includes")

# check for MPI
find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# check for SuperLU_DIST
if(ENABLE_SUPERLU)
  find_package(SUPERLUDIST REQUIRED)
endif()

# check for PETSc
if(ENABLE_PETSC)
  find_package(PETSC REQUIRED)
endif()

# check for SUNDIALS
find_package(SUNDIALS REQUIRED)

find_library(MATH_LIBRARY NAMES m)

# add the executables
if(ENABLE_PETSC)
  add_executable(cv_petsc_ex7 cv_petsc_ex7.c)
  target_link_libraries(cv_petsc_ex7 SUNDIALS::cvode SUNDIALS::nvecpetsc
                                     SUNDIALS::sunnonlinsolpetscsnes
                                     PETSC::ALL ${SUPERLU_INCLUDE_DIRS})
    if(MATH_LIBRARY)
        target_link_libraries(cv_petsc_ex7 ${MATH_LIBRARY})
    endif()
endif()
if(ENABLE_SUPERLU)
  add_executable(ark_brusselator1D_FEM_sludist ark_brusselator1D_FEM_sludist.cpp)
  target_link_libraries(ark_brusselator1D_FEM_sludist SUNDIALS::arkode SUNDIALS::nvecserial
                                                      SUNDIALS::sunmatrixslunrloc
                                                      SUNDIALS::sunlinsolsuperludist
                                                      SuperLU_DIST::SuperLU_DIST)
    if(MATH_LIBRARY)
        target_link_libraries(ark_brusselator1D_FEM_sludist ${MATH_LIBRARY})
    endif()
endif()
