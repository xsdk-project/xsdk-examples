cmake_minimum_required(VERSION 3.0)
project(petsc-example DESCRIPTION "PETSc Example" LANGUAGES C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)

set(PETSC_DIR "${PETSC_DIR}" CACHE PATH "Path to PETSc installation directory")

find_package(MPI REQUIRED)

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

find_package(PETSC REQUIRED)

if (UNIX)
    set(MATH_LIB m)
else (UNIX)
    set(MATH_LIB)
endif (UNIX)

add_executable(ex19 ex19.c)
target_link_libraries(ex19 PRIVATE PETSC::ALL ${MATH_LIB})
target_include_directories(ex19 PRIVATE ${PETSC_DIR}/include)
