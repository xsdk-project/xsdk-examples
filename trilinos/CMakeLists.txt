cmake_minimum_required(VERSION 3.8)
project( amesos-superlu
                  DESCRIPTION "Trilinos Examples"
                  LANGUAGES CXX)

set(CMAKE_CXX_COMPILER ${MPI_DIR}/bin/mpicxx)

set(MPI_INCLUDE_DIRS ${MPI_DIR}/include)

set(SUPERLU_DIR ${SPACK_BUILD_DIR}/superlu-dist-6.1.1-kodx6yljltoqu2wz4mhptt7rk57wqenv)
set(SUPERLU_INCLUDE_DIRS ${SUPERLU_DIR}/include)
set(SUPERLU_LIBRARIES ${SUPERLU_DIR}/lib)

set(TRILINOS_INCLUDE_DIRS ${TRILINOS_DIR}/include)
set(TRILINOS_LIBRARIES ${TRILINOS_DIR}/lib)

add_executable(SimpleSolve_WithParameters SimpleSolve_WithParameters.cpp)
target_compile_options(SimpleSolve_WithParameters PRIVATE -g -Wall -fopenmp)
set(CMAKE_PREFIX_PATH ${TRILINOS_LIBRARIES};${SUPERLU_INCLUDE_DIRS})
find_library(AMESOS2_LIB NAMES amesos2 PATHS ${TRILINOS_LIBRARIES})
find_library(TPETRA_LIB  NAMES tpetra PATHS ${TRILINOS_LIBRARIES})
find_library(TEUCHOSCOMM_LIB NAMES teuchoscomm PATHS ${TRILINOS_LIBRARIES})
find_library(TEUCHOSCORE_LIB NAMES teuchoscore PATHS ${TRILINOS_LIBRARIES})
find_library(TEUCHOSLIST_LIB NAMES teuchosparameterlist PATHS ${TRILINOS_LIBRARIES})
find_library(KOKKOSCORE_LIB NAMES kokkoscore PATHS ${TRILINOS_LIBRARIES})
find_library(SUPERLU_DIST_LIB NAMES superlu_dist PATHS ${SUPERLU_LIBRARIES})
target_link_libraries( SimpleSolve_WithParameters
                       ${AMESOS2_LIB} ${TPETRA_LIB} ${TEUCHOSCOMM_LIB} ${TEUCHOSCORE_LIB} ${TEUCHOSLIST_LIB} ${KOKKOSCORE_LIB} ${SUPERLU_DIST_LIB} -fopenmp -lm -ldl)

target_include_directories(SimpleSolve_WithParameters PRIVATE  ${TRILINOS_INCLUDE_DIRS} ${SUPERLU_INCLUDE_DIRS} )
