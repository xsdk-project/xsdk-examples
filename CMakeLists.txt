cmake_minimum_required(VERSION 3.21)
project(
    xsdk-examples
    DESCRIPTION "xSDK Examples"
    LANGUAGES CXX C
)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(CTest)
include(FindPackageHandleStandardArgs)
include(XsdkAddTest)

# build options
option(ENABLE_CUDA "Enable CUDA" OFF)
option(ENABLE_HIP "Enable HIP" OFF)

option(ENABLE_AMREX "Enable AMReX" ON)
set(AMREX_DIR
    "${AMREX_DIR}"
    CACHE PATH "Path to AMReX installation directory"
)

option(ENABLE_GINKGO "Enable Ginkgo" ON)
set(Ginkgo_DIR
    "${Ginkgo_DIR}"
    CACHE PATH "Path to Ginkgo installation directory"
)

option(ENABLE_HIOP "Enable HiOp" ON)
set(HIOP_DIR
    "${HIOP_DIR}"
    CACHE PATH "Path to HiOp installation directory"
)

option(ENABLE_HYPRE "Enable hypre" ON)
set(HYPRE_DIR
    "${HYPRE_DIR}"
    CACHE PATH "Path to hypre installation directory"
)

option(ENABLE_MAGMA "Enable MAGMA" OFF)
set(MAGMA_DIR
    "${MAGMA_DIR}"
    CACHE PATH "Path to MAGMA installation directory"
)

option(ENABLE_MFEM "Enable MFEM" ON)
set(MFEM_DIR
    "${MFEM_DIR}"
    CACHE PATH "Path to MFEM installation directory"
)

option(ENABLE_PETSC "Enable PETSc" ON)
set(PETSc_DIR
    "${PETSc_DIR}"
    CACHE PATH "Path to PETSc installation directory"
)

option(ENABLE_PLASMA "Enable PLASMA" ON)
set(PLASMA_DIR
    "${PLASMA_DIR}"
    CACHE PATH "Path to PLASMA installation directory"
)

option(ENABLE_PUMI "Enable PUMI" ON)
set(PUMI_DIR
    "${PUMI_DIR}"
    CACHE PATH "Path to PUMI installation directory"
)

option(ENABLE_SUNDIALS "Enable SUNDIALS" ON)
set(SUNDIALS_DIR
    "${SUNDIALS_DIR}"
    CACHE PATH "Path to SUNDIALS installation directory"
)

option(ENABLE_SUPERLU "Enable SuperLU" ON)
set(SUPERLUDIST_DIR
    "${SUPERLUDIST_DIR}"
    CACHE PATH "Path to SuperLU_DIST installation directory"
)

option(ENABLE_STRUMPACK "Enable STRUMPACK" OFF)
set(STRUMPACK_DIR
    "${STRUMPACK_DIR}"
    CACHE PATH "Path to STRUMPACK installation directory"
)

option(ENABLE_TRILINOS "Enable TRILINOS" OFF)
set(TRILINOS_DIR
    "${Trilinos_DIR}"
    CACHE PATH "Path to Trilinos installation directory"
)

set(METIS_DIR
    "${METIS_DIR}"
    CACHE PATH "Path to Metis installation directory"
)

# Don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# When building, don't use the install RPATH already (but later on when
# installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_FULL_LIBDIR}")

# Add the automatically determined parts of the RPATH which point to directories
# outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# The RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
endif()

# check for MPI
find_package(MPI REQUIRED)

# check for OpenMP
find_package(OpenMP)

# compiler options
if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

# setup CUDA
if(ENABLE_CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    find_package(CUDAToolkit REQUIRED)
endif()

if(ENABLE_HIP)
    enable_language(HIP)
    find_package(hip REQUIRED)
    find_package(hipsparse REQUIRED)
    find_package(hiprand REQUIRED)
    find_package(rocrand REQUIRED)
    find_package(rocprim REQUIRED)
    find_package(rocsparse REQUIRED)
endif()

# check for AMReX
if(ENABLE_AMREX)
    find_package(AMReX REQUIRED)
endif()

# check for HiOp
if(ENABLE_HIOP)
    find_package(HiOp REQUIRED)
endif()

# check for hypre
if(ENABLE_HYPRE)
    find_package(HYPRE REQUIRED)
endif()

# check for MAGMA
if(ENABLE_MAGMA)
    find_package(MAGMA REQUIRED)
endif()

# check for MFEM
if(ENABLE_MFEM)
    find_package(ZLIB REQUIRED)
    find_package(MFEM REQUIRED)
    if(ENABLE_GINKGO)
        find_package(Ginkgo REQUIRED)
    endif()
endif()

# check for PETSC
if(ENABLE_PETSC)
    find_package(PETSc REQUIRED)
endif()

# check for PLASMA
if(ENABLE_PLASMA)
    find_package(PLASMA REQUIRED)
    find_package(BLASPP REQUIRED)
    find_package(LAPACKPP REQUIRED)
    find_package(SLATE REQUIRED)
endif()

# check for PUMI
if(ENABLE_PUMI)
    # not really needed: MFEM links with PUMI if it is built with it
    # find_package(SCOREC REQUIRED)
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
if(ENABLE_PLASMA)
    add_subdirectory(plasma)
endif()
if(ENABLE_STRUMPACK)
    add_subdirectory(strumpack)
endif()
if(ENABLE_SUNDIALS)
    add_subdirectory(sundials)
endif()
if(ENABLE_TRILINOS)
    add_subdirectory(trilinos)
endif()
