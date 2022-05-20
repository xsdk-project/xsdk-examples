# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# SuperLUDIST find module that creates an imported target for
# SuperLU_DIST. The target is XSDK::SUPERLU.
#
# This module also defines variables, but it is best to use
# the defined target to ensure includes and compile/link
# options are correctly passed to consumers.
#
#   SUPERLUDIST_FOUND        - system has SuperLU_DIST library
#   SUPERLUDIST_INCLUDE_DIR  - Directory with the SuperLU_DIST header
#   SUPERLUDIST_LIBRARY      - the SuperLU_DIST library
# ---------------------------------------------------------------

### find include dir
find_path(SUPERLUDIST_INCLUDE_DIR superlu_defs.h
  HINTS ${SUPERLU_DIR} $ENV{SUPERLU_DIR} ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with the SuperLU DIST header")

### find library
find_library(SUPERLUDIST_LIBRARY superlu_dist
  HINTS ${SUPERLU_DIR} $ENV{SUPERLU_DIR} ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
  DOC "The SuperLU DIST library")

# set package variables including SUPERLUDIST_FOUND
find_package_handle_standard_args(SUPERLUDIST
  REQUIRED_VARS
    SUPERLUDIST_LIBRARY
    SUPERLUDIST_INCLUDE_DIR
  )

# Create target for SuperLU_DIST
if(SUPERLUDIST_FOUND AND (NOT TARGET XSDK::SUPERLU))
  add_library(XSDK::SUPERLU UNKNOWN IMPORTED)
  set_target_properties(XSDK::SUPERLU PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SUPERLUDIST_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${SUPERLUDIST_LIBRARIES}"
    IMPORTED_LOCATION "${SUPERLUDIST_LIBRARY}")
endif()
