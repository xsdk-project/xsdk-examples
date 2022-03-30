cmake_minimum_required(VERSION 3.12)
project(
  hypre-superludist
  DESCRIPTION "HYPRE + SuperLU_DIST Example"
  LANGUAGES C)

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

add_executable(ij_laplacian ij_laplacian.c)

target_link_libraries(ij_laplacian PRIVATE XSDK::HYPRE XSDK::SUPERLU)
if(MATH_LIBRARY)
  target_link_libraries(ij_laplacian PRIVATE ${MATH_LIBRARY})
endif()

install(TARGETS ij_laplacian RUNTIME DESTINATION bin)

# Copy glvis helper script and directory
file(COPY vis DESTINATION .)

install(
  DIRECTORY vis
  DESTINATION share/xsdk-examples/hypre
  PATTERN
    vis/*
    PERMISSIONS
      OWNER_EXECUTE
      OWNER_WRITE
      OWNER_READ
      GROUP_EXECUTE
      GROUP_READ
      WORLD_READ
      WORLD_EXECUTE)
