cmake_minimum_required(VERSION 3.0)
project(
  petsc-example
  DESCRIPTION "PETSc Example"
  LANGUAGES C)

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

add_executable(ex19 ex19.c)
target_link_libraries(ex19 PRIVATE XSDK::PETSc)
target_include_directories(ex19 PRIVATE ${PETSC_DIR}/include)

if(MATH_LIBRARY)
  target_link_libraries(ex19 PRIVATE ${MATH_LIBRARY})
endif()

install(TARGETS ex19 RUNTIME DESTINATION bin)
