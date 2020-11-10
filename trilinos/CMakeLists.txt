cmake_minimum_required(VERSION 3.12)
project( amesos-superlu
                  DESCRIPTION "Trilinos Examples"
                  LANGUAGES CXX)

set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

set(TRILINOS_INCLUDE_DIRS ${TRILINOS_DIR}/include)
set(TRILINOS_LIBRARIES ${TRILINOS_DIR}/lib)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif()

add_executable(SimpleSolve_WithParameters SimpleSolve_WithParameters.cpp)
target_compile_options(SimpleSolve_WithParameters PRIVATE ${OPENMP_C_FLAGS})
target_link_libraries(SimpleSolve_WithParameters PRIVATE XSDK::SUPERLU METIS)

set(CMAKE_PREFIX_PATH ${TRILINOS_LIBRARIES})
find_library(AMESOS2_LIB NAMES amesos2 PATHS ${TRILINOS_LIBRARIES})
find_library(TPETRA_LIB  NAMES tpetra PATHS ${TRILINOS_LIBRARIES})
find_library(TEUCHOSCOMM_LIB NAMES teuchoscomm PATHS ${TRILINOS_LIBRARIES})
find_library(TEUCHOSCORE_LIB NAMES teuchoscore PATHS ${TRILINOS_LIBRARIES})
find_library(TEUCHOSLIST_LIB NAMES teuchosparameterlist PATHS ${TRILINOS_LIBRARIES})
find_library(KOKKOSCORE_LIB NAMES kokkoscore PATHS ${TRILINOS_LIBRARIES})

target_link_libraries(SimpleSolve_WithParameters PRIVATE
  ${AMESOS2_LIB} ${TPETRA_LIB} ${TEUCHOSCOMM_LIB} ${TEUCHOSCORE_LIB} ${TEUCHOSLIST_LIB} ${KOKKOSCORE_LIB})

target_include_directories(SimpleSolve_WithParameters PRIVATE ${TRILINOS_INCLUDE_DIRS})

if(MATH_LIBRARY)
  target_link_libraries(SimpleSolve_WithParameters PRIVATE ${MATH_LIBRARY})
endif()

find_library(DYNLNK_LIBRARY NAMES dl)
if(DYNLNK_LIBRARY)
  target_link_libraries(SimpleSolve_WithParameters PRIVATE ${DYNLNK_LIBRARY})
endif()

install(TARGETS SimpleSolve_WithParameters RUNTIME DESTINATION bin)
