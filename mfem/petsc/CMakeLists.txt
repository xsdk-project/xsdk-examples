cmake_minimum_required(VERSION 3.8)
project( obstacle
         DESCRIPTION "Obstacle problem using MFEM + PETSC/Tao"
         LANGUAGES CXX)

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

add_executable(obstacle obstacle.cpp)
target_link_libraries(obstacle PRIVATE XSDK::MFEM XSDK::HYPRE XSDK::SUNDIALS PETSC::ALL ZLIB::ZLIB)

install(TARGETS obstacle RUNTIME DESTINATION bin)
