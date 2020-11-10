cmake_minimum_required(VERSION 3.12)
project(transient-heat
        DESCRIPTION "Transient Heat Conduction using MFEM + SUNDIALS + HYPRE"
        LANGUAGES C)

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

add_executable(transient-heat transient-heat.cpp)
target_link_libraries(transient-heat PRIVATE XSDK::MFEM XSDK::SUNDIALS XSDK::HYPRE PETSC::ALL)

install(TARGETS transient-heat RUNTIME DESTINATION bin)
