# ---------------------------------------------------------------
# Programmer: Cody J. Balos, David J. Gardner  @ LLNL
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

# determine SUNDIALS components needed
if(NOT SUNDIALS_FIND_COMPONENTS)
  set(SUNDIALS_FIND_COMPONENTS
    "arkode"
    "cvodes"
    "nvecserial"
    "nvecparallel"
    "nvecmpiplusx")
endif()

if(ENABLE_CUDA)
  list(APPEND SUNDIALS_FIND_COMPONENTS "nveccuda")
endif()
if(ENABLE_MFEM)
  list(APPEND SUNDIALS_FIND_COMPONENTS "kinsol")
endif()
if(ENABLE_MAGMA)
  list(APPEND SUNDIALS_FIND_COMPONENTS "sunmatrixmagmadense")
  list(APPEND SUNDIALS_FIND_COMPONENTS "sunlinsolmagmadense")
endif()
if(ENABLE_PETSC)
  list(APPEND SUNDIALS_FIND_COMPONENTS "nvecpetsc")
  list(APPEND SUNDIALS_FIND_COMPONENTS "sunnonlinsolpetscsnes")
endif()
if(ENABLE_SUPERLU)
  list(APPEND SUNDIALS_FIND_COMPONENTS "sunmatrixslunrloc")
  list(APPEND SUNDIALS_FIND_COMPONENTS "sunlinsolsuperludist")
endif()

find_package(SUNDIALS REQUIRED COMPONENTS ${SUNDIALS_FIND_COMPONENTS}
    HINTS ${SUNDIALS_DIR} $ENV{SUNDIALS_DIR} ${CMAKE_PREFIX_PATH}
    NO_DEFAULT_PATH)

if(NOT TARGET XSDK::SUNDIALS)
  add_library(XSDK_SUNDIALS INTERFACE)
  foreach(_component ${SUNDIALS_FIND_COMPONENTS})
      target_link_libraries(XSDK_SUNDIALS INTERFACE SUNDIALS::${_component})
  endforeach()
  add_library(XSDK::SUNDIALS ALIAS XSDK_SUNDIALS)
endif()
