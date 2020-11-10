cmake_minimum_required(VERSION 3.8)
project( convdiff
         DESCRIPTION "Convective diffusion using MFEM + SUPERLU + HYPRE"
         LANGUAGES CXX)

set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

add_executable(convdiff convdiff.cpp)
target_link_libraries(convdiff PRIVATE PETSC::ALL SUNDIALS::cvodes XSDK::HYPRE XSDK::SUPERLU XSDK::MFEM)

install(TARGETS convdiff RUNTIME DESTINATION bin)
