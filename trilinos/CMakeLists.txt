cmake_minimum_required(VERSION 3.8)
project( amesos-superlu
                  DESCRIPTION "Trilinos Examples"
                  LANGUAGES CXX)

# location of additional CMake modules
#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)
set(SPACK_BUILD_DIR "/home/knteran/test/spack/opt/spack/linux-rhel7-sandybridge/gcc-7.2.0")
#set(TRILINOS_DIR "/path/to/trilinos")
#set(SUPERLU_DIR "/path/to/superlu")
#set(MPI_DIR "path/to/mpi")

set(MPI_DIR ${SPACK_BUILD_DIR}/openmpi-3.1.4-sr7eobrwauqv7ija4pyxpto4mgapjghv/)
set(CMAKE_CXX_COMPILER ${MPI_DIR}/bin/mpicxx)

set(MPI_INCLUDE_DIRS ${MPI_DIR}/include)
set(SUPERLU_DIR ${SPACK_BUILD_DIR}/superlu-dist-6.1.1-kodx6yljltoqu2wz4mhptt7rk57wqenv)
set(SUPERLU_INCLUDE_DIRS ${SUPERLU_DIR}/include)
set(SUPERLU_LIBRARIES ${SUPERLU_DIR}/lib)


set(TRILINOS_DIR ${SPACK_BUILD_DIR}/trilinos-12.18-pf3jgq5nlgxgtjpbx6a7v5icqdotj7iu)
set(TRILINOS_INCLUDE_DIRS ${TRILINOS_DIR}/include)
set(TRILINOS_LIBRARIES ${TRILINOS_DIR}/lib)

#find_package(SUPERLUDIST REQUIRED)
#find_package(TRILINOS  REQUIRED)

add_executable(SimpleSolve_WithParameters SimpleSolve_WithParameters.cpp)
target_compile_options(SimpleSolv_WithParameterse PRIVATE -g -Wall -fopenmp)
set(CMAKE_PREFIX_PATH ${TRILINOS_LIBRARIES};${SUPERLU_INCLUDE_DIRS})
find_library(AMESOS2_LIB NAMES amesos2 PATHS ${TRILINOS_LIBRARIES})
find_library(TPETRA_LIB  NAMES tpetra PATHS ${TRILINOS_LIBRARIES})
find_library(TEUCHOSCOMM_LIB NAMES teuchoscomm PATHS ${TRILINOS_LIBRARIES})
find_library(TEUCHOSCORE_LIB NAMES teuchoscore PATHS ${TRILINOS_LIBRARIES})
find_library(TEUCHOSLIST_LIB NAMES teuchosparameterlist PATHS ${TRILINOS_LIBRARIES})
find_library(SUPERLU_DIST_LIB NAMES superlu_dist PATHS ${SUPERLU_LIBRARIES})

                       ${AMESOS2_LIB} ${TPETRA_LIB} ${TEUCHOSCOMM_LIB} ${TEUCHOSCORE_LIB}  ${SUPERLU_DIST_LIB} -fopenmp -lm -ldl)
target_link_libraries( SimpleSolve_WithParameters
                       ${AMESOS2_LIB} ${TPETRA_LIB} ${TEUCHOSCOMM_LIB}  ${TEUCHOSCORE_LIB}  ${TEUCHOSLIST_LIB} ${SUPERLU_DIST_LIB} -fopenmp -lm -ldl)

target_include_directories(SimpleSolve_WithParameters PRIVATE  ${TRILINOS_INCLUDE_DIRS} ${SUPERLU_INCLUDE_DIRS} )
