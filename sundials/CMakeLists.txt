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

cmake_minimum_required(VERSION 3.12)

# set the project name and default languages
project(sundials-xsdk-examples CXX C)

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER})

# add the executables
if(ENABLE_PETSC)
  add_executable(cv_petsc_ex7 cv_petsc_ex7.c)
  target_link_libraries(cv_petsc_ex7 SUNDIALS::cvode SUNDIALS::sunnonlinsolpetscsnes PETSC::ALL)
  if(MATH_LIBRARY)
      target_link_libraries(cv_petsc_ex7 ${MATH_LIBRARY})
  endif()
  install(TARGETS cv_petsc_ex7 RUNTIME DESTINATION bin)
endif()
if(ENABLE_SUPERLU)
  add_executable(ark_brusselator1D_FEM_sludist ark_brusselator1D_FEM_sludist.cpp)
  target_link_libraries(ark_brusselator1D_FEM_sludist PRIVATE SUNDIALS::arkode SUNDIALS::sunlinsolsuperludist XSDK::SUPERLU)
  if(MATH_LIBRARY)
    target_link_libraries(ark_brusselator1D_FEM_sludist PRIVATE ${MATH_LIBRARY})
  endif()
  install(TARGETS ark_brusselator1D_FEM_sludist RUNTIME DESTINATION bin)
endif()
if(ENABLE_MAGMA)
  set_source_files_properties(cvRoberts_blockdiag_magma.cpp PROPERTIES LANGUAGE CUDA)
  add_executable(cvRoberts_blockdiag_magma cvRoberts_blockdiag_magma.cpp)
  target_link_libraries(cvRoberts_blockdiag_magma PRIVATE SUNDIALS::cvode SUNDIALS::sunlinsolmagma XSDK::MAGMA)
  target_compile_definitions(cvRoberts_blockdiag_magma PRIVATE SUNDIALS_MAGMA_BACKENDS_CUDA)
  if(MATH_LIBRARY)
    target_link_libraries(cvRoberts_blockdiag_magma PRIVATE ${MATH_LIBRARY})
  endif()
  install(TARGETS cvRoberts_blockdiag_magma RUNTIME DESTINATION bin)
endif()

