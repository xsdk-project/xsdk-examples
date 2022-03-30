cmake_minimum_required(VERSION 3.8)
project(
  obstacle
  DESCRIPTION "Obstacle problem using MFEM + PETSC/Tao"
  LANGUAGES CXX)

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

add_executable(obstacle obstacle.cpp)
target_link_libraries(obstacle PRIVATE XSDK::MFEM XSDK::HYPRE XSDK::PETSc
                                       XSDK::SUNDIALS ZLIB::ZLIB)
if(ENABLE_CUDA)
  # find_package(Thrust REQUIRED CONFIG) thrust_create_target(Thrust)
  # target_link_libraries(obstacle PRIVATE Thrust)
endif()

install(TARGETS obstacle RUNTIME DESTINATION bin)
