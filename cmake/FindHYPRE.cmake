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
# HYPRE find module that creates an imported target for HYPRE.
# The target is XSDK::HYPRE.
#
# This module also defines variables, but it is best to use
# the defined target to ensure includes and compile/link
# options are correctly passed to consumers.
#
#   HYPRE_FOUND       - system has HYPRE library
#   HYPRE_LIBRARY     - the HYPRE library
#   HYPRE_INCLUDE_DIR - the HYPRE include path
#   HYPRE_LIBRARIES   - all of the libraries needed for HYPRE
# ---------------------------------------------------------------

### find include dir
find_path(HYPRE_INCLUDE_DIR
  NAMES HYPRE.h hypre.h
  HINTS ${HYPRE_DIR} $ENV{HYPRE_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with hypre header.")

### find library
find_library(HYPRE_LIBRARY
  NAMES HYPRE hypre
  HINTS ${HYPRE_DIR} $ENV{HYPRE_DIR}
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
  DOC "The hypre library.")

list(FIND HYPRE_LIBRARIES ${HYPRE_LIBRARY} _idx)
if (_idx EQUAL -1)
  set(HYPRE_LIBRARIES "${HYPRE_LIBRARY};${HYPRE_LIBRARIES}" CACHE STRING "" FORCE)
endif ()

# set package variables including HYPRE_FOUND
find_package_handle_standard_args(HYPRE
  REQUIRED_VARS
    HYPRE_LIBRARY
    HYPRE_LIBRARIES
    HYPRE_INCLUDE_DIR
  )

# Create target for HYPRE
if(HYPRE_FOUND)

  if(NOT TARGET XSDK::HYPRE)
    add_library(XSDK::HYPRE UNKNOWN IMPORTED)
  endif()

  set_target_properties(XSDK::HYPRE PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${HYPRE_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${HYPRE_LIBRARIES}"
    IMPORTED_LOCATION "${HYPRE_LIBRARY}")

endif()
