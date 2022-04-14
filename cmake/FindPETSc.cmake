# FindPETSc
# ---------
#
# Locates the PETSc library using pkg-config
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
# 
# This module defines the followwing IMPORTED target:
#
#  PETSc::PETSc        - the PETSc library
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#
#  PETSc_FOUND          - if false, do not try to link to PETSc
#  PETSc_LIBRARIES      - a list of the full paths to all libraries
#  PETSc_INCLUDE_DIRS   - a list of all include directories
#  PETSc_VERSION        - the full version of PETSc MAJOR.MINOR.PATCH
#  PETSc_VERSION_MAJOR  - the MAJOR part of PETSc_VERSION
#  PETSc_VERSION_MINOR  - the MINOR part of PETSc_VERSION
#  PETSc_VERSION_PATCH  - the PATCH part of PETSc_VERSION
# 
# Variables for locating PETSc
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Additional CMake variables for locating PETSc
#  PETSc_DIR            - the path to the root directory of PETSc 
#  PETSc_ARCH           - the PETSc architecture
#  PETSc_NO_ENV         - instructs the module not to use the environment variables 'PETSC_DIR' and 'PETSC_ENV' to find PETSc
#
# Environment Variables for locating PETSc
#  PETSC_DIR            - the path to the root directory of PETSc, part of the PETSc installation process
#  PETSC_ARCH           - the PETSc architecture, part of the PETSc installation process
#
#
# The orignal author is Frédéric Simonis @fsimonis.
# This file is derived from https://github.com/precice/precice/blob/develop/cmake/modules/FindPETSc.cmake,
# as such it falls under the original license, copied below:
#
#                    GNU LESSER GENERAL PUBLIC LICENSE
#                        Version 3, 29 June 2007
#
#  Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
#  Everyone is permitted to copy and distribute verbatim copies
#  of this license document, but changing it is not allowed.
#
#
#   This version of the GNU Lesser General Public License incorporates
# the terms and conditions of version 3 of the GNU General Public
# License, supplemented by the additional permissions listed below.
#
#   0. Additional Definitions.
#
#   As used herein, "this License" refers to version 3 of the GNU Lesser
# General Public License, and the "GNU GPL" refers to version 3 of the GNU
# General Public License.
#
#   "The Library" refers to a covered work governed by this License,
# other than an Application or a Combined Work as defined below.
#
#   An "Application" is any work that makes use of an interface provided
# by the Library, but which is not otherwise based on the Library.
# Defining a subclass of a class defined by the Library is deemed a mode
# of using an interface provided by the Library.
#
#   A "Combined Work" is a work produced by combining or linking an
# Application with the Library.  The particular version of the Library
# with which the Combined Work was made is also called the "Linked
# Version".
#
#   The "Minimal Corresponding Source" for a Combined Work means the
# Corresponding Source for the Combined Work, excluding any source code
# for portions of the Combined Work that, considered in isolation, are
# based on the Application, and not on the Linked Version.
#
#   The "Corresponding Application Code" for a Combined Work means the
# object code and/or source code for the Application, including any data
# and utility programs needed for reproducing the Combined Work from the
# Application, but excluding the System Libraries of the Combined Work.
#
#   1. Exception to Section 3 of the GNU GPL.
#
#   You may convey a covered work under sections 3 and 4 of this License
# without being bound by section 3 of the GNU GPL.
#
#   2. Conveying Modified Versions.
#
#   If you modify a copy of the Library, and, in your modifications, a
# facility refers to a function or data to be supplied by an Application
# that uses the facility (other than as an argument passed when the
# facility is invoked), then you may convey a copy of the modified
# version:
#
#    a) under this License, provided that you make a good faith effort to
#    ensure that, in the event an Application does not supply the
#    function or data, the facility still operates, and performs
#    whatever part of its purpose remains meaningful, or
#
#    b) under the GNU GPL, with none of the additional permissions of
#    this License applicable to that copy.
#
#   3. Object Code Incorporating Material from Library Header Files.
#
#   The object code form of an Application may incorporate material from
# a header file that is part of the Library.  You may convey such object
# code under terms of your choice, provided that, if the incorporated
# material is not limited to numerical parameters, data structure
# layouts and accessors, or small macros, inline functions and templates
# (ten or fewer lines in length), you do both of the following:
#
#    a) Give prominent notice with each copy of the object code that the
#    Library is used in it and that the Library and its use are
#    covered by this License.
#
#    b) Accompany the object code with a copy of the GNU GPL and this license
#    document.
#
#   4. Combined Works.
#
#   You may convey a Combined Work under terms of your choice that,
# taken together, effectively do not restrict modification of the
# portions of the Library contained in the Combined Work and reverse
# engineering for debugging such modifications, if you also do each of
# the following:
#
#    a) Give prominent notice with each copy of the Combined Work that
#    the Library is used in it and that the Library and its use are
#    covered by this License.
#
#    b) Accompany the Combined Work with a copy of the GNU GPL and this license
#    document.
#
#    c) For a Combined Work that displays copyright notices during
#    execution, include the copyright notice for the Library among
#    these notices, as well as a reference directing the user to the
#    copies of the GNU GPL and this license document.
#
#    d) Do one of the following:
#
#        0) Convey the Minimal Corresponding Source under the terms of this
#        License, and the Corresponding Application Code in a form
#        suitable for, and under terms that permit, the user to
#        recombine or relink the Application with a modified version of
#        the Linked Version to produce a modified Combined Work, in the
#        manner specified by section 6 of the GNU GPL for conveying
#        Corresponding Source.
#
#        1) Use a suitable shared library mechanism for linking with the
#        Library.  A suitable mechanism is one that (a) uses at run time
#        a copy of the Library already present on the user's computer
#        system, and (b) will operate properly with a modified version
#        of the Library that is interface-compatible with the Linked
#        Version.
#
#    e) Provide Installation Information, but only if you would otherwise
#    be required to provide such information under section 6 of the
#    GNU GPL, and only to the extent that such information is
#    necessary to install and execute a modified version of the
#    Combined Work produced by recombining or relinking the
#    Application with a modified version of the Linked Version. (If
#    you use option 4d0, the Installation Information must accompany
#    the Minimal Corresponding Source and Corresponding Application
#    Code. If you use option 4d1, you must provide the Installation
#    Information in the manner specified by section 6 of the GNU GPL
#    for conveying Corresponding Source.)
#
#   5. Combined Libraries.
#
#   You may place library facilities that are a work based on the
# Library side by side in a single library together with other library
# facilities that are not Applications and are not covered by this
# License, and convey such a combined library under terms of your
# choice, if you do both of the following:
#
#    a) Accompany the combined library with a copy of the same work based
#    on the Library, uncombined with any other library facilities,
#    conveyed under the terms of this License.
#
#    b) Give prominent notice with the combined library that part of it
#    is a work based on the Library, and explaining where to find the
#    accompanying uncombined form of the same work.
#
#   6. Revised Versions of the GNU Lesser General Public License.
#
#   The Free Software Foundation may publish revised and/or new versions
# of the GNU Lesser General Public License from time to time. Such new
# versions will be similar in spirit to the present version, but may
# differ in detail to address new problems or concerns.
#
#   Each version is given a distinguishing version number. If the
# Library as you received it specifies that a certain numbered version
# of the GNU Lesser General Public License "or any later version"
# applies to it, you have the option of following the terms and
# conditions either of that published version or of any later version
# published by the Free Software Foundation. If the Library as you
# received it does not specify a version number of the GNU Lesser
# General Public License, you may choose any version of the GNU Lesser
# General Public License ever published by the Free Software Foundation.
#
#   If the Library as you received it specifies that a proxy can decide
# whether future versions of the GNU Lesser General Public License shall
# apply, that proxy's public statement of acceptance of any version is
# permanent authorization for you to choose that version for the
# Library.
#


cmake_policy(VERSION 3.10)


# Macro to print the search context used by pkg-config
macro(_petsc_print_pkg_env)
  if(NOT PETSc_FIND_QUIETLY)
    set(_env_mess "pkg-config will search the following paths:")
    if(DEFINED ENV{PKG_CONFIG_PATH})
      set(_env_mess "${_env_mess}\nPKG_CONFIG_PATH")
      string(REPLACE ":" ";" _env_pkg_list "$ENV{PKG_CONFIG_PATH}")
      foreach(p IN LISTS _env_pkg_list)
        set(_env_mess "${_env_mess}\n   ${p}")
      endforeach()
      unset(_env_pkg_list)
    endif()
    if(CMAKE_PREFIX_PATH)
      set(_env_mess "${_env_mess}\nCMAKE_PREFIX_PATH")
      foreach(p IN LISTS CMAKE_PREFIX_PATH)
        set(_env_mess "${_env_mess}\n   ${p}")
      endforeach()
    endif()
    if(CMAKE_FRAMEWORK_PATH)
      set(_env_mess "${_env_mess}\nCMAKE_FRAMEWORK_PATH")
      foreach(p IN LISTS CMAKE_FRAMEWORK_PATH)
        set(_env_mess "${_env_mess}\n   ${p}")
      endforeach()
    endif()
    if(CMAKE_APPBUNDLE_PATH)
      set(_env_mess "${_env_mess}\nCMAKE_APPBUNDLE_PATH")
      foreach(p IN LISTS CMAKE_APPBUNDLE_PATH)
        set(_env_mess "${_env_mess}\n   ${p}")
      endforeach()
    endif()
    message(STATUS "${_env_mess}")
    unset(_env_mess)
  endif()
endmacro()


# Message macro which respects the QUIET arguemnt of the package
macro(_message)
  if(NOT PETSc_FIND_QUIETLY)
    message(${ARGV})
  endif()
endmacro()

set(_petsc_quiet_arg "")
if(PETSc_FIND_QUIETLY)
  set(_petsc_quiet_arg "QUIET")
endif()
find_package(PkgConfig ${_petsc_quiet_arg})

if(PKG_CONFIG_FOUND)
  # Detect additional pefix paths
  set(_petsc_detected_prefixes "")
  if(DEFINED PETSc_DIR)
    list(APPEND _petsc_detected_prefixes "${PETSc_DIR}")
    if(DEFINED PETSc_ARCH)
      list(APPEND _petsc_detected_prefixes "${PETSc_DIR}/${PETSc_ARCH}")
    endif()
  endif()

  if(DEFINED ENV{PETSC_DIR} AND NOT PETSc_NO_ENV)
    list(APPEND _petsc_detected_prefixes "$ENV{PETSC_DIR}")
    if(DEFINED ENV{PETSC_ARCH})
      list(APPEND _petsc_detected_prefixes "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}")
    endif()
  endif()
  list(REMOVE_DUPLICATES _petsc_detected_prefixes)

  set(_petsc_prefixes "")
  set(_petsc_skipped_prefixes "")
  _message(STATUS "Detecting additional PETSc prefixes")
  foreach(prefix IN LISTS _petsc_detected_prefixes )
    if(EXISTS "${prefix}/lib/pkgconfig")
      _message(STATUS "Detected ${prefix}")
      list(APPEND _petsc_prefixes "${prefix}")
    else()          
      list(APPEND _petsc_skipped_prefixes "${prefix}")
    endif()
  endforeach()
  if(_petsc_skipped_prefixes)
    _message(STATUS "Skipped the following invalid prefixes: ${_petsc_skipped_prefixes}")
  endif()
  unset(_petsc_detected_prefixes)

  # Remember the previous state of CMAKE_PREFIX_PATH
  set(_petsc_prefix_unset True)
  if(DEFINED CMAKE_PREFIX_PATH)
    set(_petsc_prefix_unset False)
    set(_petsc_prefix_old ${CMAKE_PREFIX_PATH})
  endif()

  list(APPEND CMAKE_PREFIX_PATH ${_petsc_prefixes})
  _petsc_print_pkg_env()

  # Build the pkg-config version spec
  set(_pkg_version_spec "")
  if(DEFINED PETSc_FIND_VERSION)
    if(PETSc_FIND_VERSION_EXACT)
      set(_pkg_version_spec "=${PETSc_FIND_VERSION}")
    else()
      set(_pkg_version_spec ">=${PETSc_FIND_VERSION}")
    endif()
  endif()

  # Set PKG_CONFIG_ALLOW_SYSTEM_CFLAGS
  set(_petsc_prev_allow_system_cflags $ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS})
  set(ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS} 1)

  # Use pkg-config to find PETSc
  set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH "YES")
  pkg_check_modules(PC_PETSc ${_petsc_quiet_arg} "PETSc${_pkg_version_spec}")

  # Restore/Reset PKG_CONFIG_USE_CMAKE_PREFIX_PATH
  set(ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS} ${_petsc_prev_allow_system_cflags})

  unset(PKG_CONFIG_USE_CMAKE_PREFIX_PATH)
  unset(_pkg_version_spec)

  # Restore the previous state of CMAKE_PREFIX_PATH
  if(_petsc_prefix_unset)
    unset(CMAKE_PREFIX_PATH)
  else()
    set(CMAKE_PREFIX_PATH "${_petsc_prefix_old}")
  endif()

  # Set straight forward result variables
  set(PETSc_FOUND ${PC_PETSc_FOUND})
  set(PETSc_INCLUDE_DIRS ${PC_PETSc_INCLUDE_DIRS})

  # libm is always required
  set(_petsc_libs "m")
  set(_petsc_missing_libs "")

  # Find main PETSc libraries
  foreach(_next_lib IN LISTS PC_PETSc_LIBRARIES)
    find_library(_petsc_lib_${_next_lib} NAMES ${_next_lib} HINTS ${PC_PETSc_LIBRARY_DIRS})
    if(_petsc_lib_${_next_lib})
      list(APPEND _petsc_libs "${_petsc_lib_${_next_lib}}")
    else()
      list(APPEND _petsc_missing_libs "${_next_lib}")
    endif()
  endforeach()

  # Link against MPI if it is used.
  # This adds all required link directories.
  foreach(_next_lib IN LISTS PC_PETSc_STATIC_LIBRARIES)
    if(_next_lib STREQUAL "mpi")
      find_package(MPI ${_petsc_quiet_arg})
      if(MPI_FOUND)
        # Prefer to use the CXX dependencies if enabled otherwise use the C
        if(DEFINED CMAKE_CXX_COMPILER)
          list(APPEND _petsc_libs "MPI::MPI_CXX")
        else()
          enable_language(C)
          list(APPEND _petsc_libs "MPI::MPI_C")
        endif()
        break()
      else()
        list(APPEND _petsc_missing_libs "MPI")
      endif()
    endif()
  endforeach()

  # Check if everything was detected
  if(_petsc_missing_libs AND NOT PETSc_FIND_QUIETLY)
    message("The following libraries were not detected: ${_petsc_missing_libs}")
  elseif(NOT _petsc_missing_libs)
    # Set the visible variable. This will let the module to succeed.
    set(PETSc_LIBRARIES "${_petsc_libs}")
  endif()
  unset(_petsc_libs)
  unset(_petsc_missing_libs)

  # Extract version parts from the version information
  if(PC_PETSc_VERSION)
    set(_petsc_versions "")
    string(REGEX MATCHALL "[0-9]+" _petsc_versions ${PC_PETSc_VERSION})
    list(GET _petsc_versions 0 _petsc_version_major)
    list(GET _petsc_versions 1 _petsc_version_minor)
    list(GET _petsc_versions 2 _petsc_version_patch)

    set(PETSc_VERSION ${PC_PETSc_VERSION} CACHE STRING "Full version of PETSc")
    set(PETSc_VERSION_MAJOR ${_petsc_version_major} CACHE INTERNAL "Major version of PETSc")
    set(PETSc_VERSION_MINOR ${_petsc_version_minor} CACHE INTERNAL "Minor version of PETSc")
    set(PETSc_VERSION_PATCH ${_petsc_version_patch} CACHE INTERNAL "Patch version of PETSc")

    unset(_petsc_versions)
    unset(_petsc_version_major)
    unset(_petsc_version_minor)
    unset(_petsc_version_patch)
  endif()
endif()
unset(_petsc_quiet_arg)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PETSc
  REQUIRED_VARS PETSc_FOUND PETSc_INCLUDE_DIRS PETSc_LIBRARIES
  VERSION_VAR PETSc_VERSION
  )

if(NOT TARGET XSDK::PETSc)
  add_library(XSDK::PETSc INTERFACE IMPORTED)
  set_target_properties(XSDK::PETSc PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${PETSc_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${PETSc_LIBRARIES}"
    )
endif()

mark_as_advanced(PETSc_INCLUDE_DIRS PETSc_LIBRARIES PETSc_VERSION_MAJOR PETSc_VERSION_MINOR PETSc_VERSION_PATCH VERSION_VAR PETSc_VERSION)
