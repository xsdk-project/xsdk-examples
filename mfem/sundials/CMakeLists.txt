cmake_minimum_required(VERSION 3.18)
project(
  transient-heat
  DESCRIPTION "Transient Heat Conduction using MFEM + SUNDIALS + HYPRE"
  LANGUAGES C CXX CUDA)

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})

add_executable(transient-heat transient-heat.cpp)
target_link_libraries(transient-heat PRIVATE XSDK::SUNDIALS XSDK::PETSc
                                             XSDK::HYPRE XSDK::MFEM ZLIB::ZLIB)
install(TARGETS transient-heat RUNTIME DESTINATION bin)

if(ENABLE_CUDA)
  add_executable(advection advection.cpp)
  target_link_libraries(advection PRIVATE XSDK::PETSc XSDK::SUNDIALS
                                          XSDK::HYPRE XSDK::MFEM ZLIB::ZLIB)
  install(TARGETS advection RUNTIME DESTINATION bin)
endif()
