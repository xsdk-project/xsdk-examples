cmake_minimum_required(VERSION 3.8)
project(mfem-gko LANGUAGES CXX)

find_package(OpenMP)

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

option(ENABLE_MFEM "Enable MFEM" ON)

add_executable(mfem_ex1_gko mfem_ex1_gko.cpp)
target_sources(mfem_ex1_gko PRIVATE mfem_ex1_gko.cpp mfem_wrapper.cpp)

target_link_libraries(mfem_ex1_gko PRIVATE XSDK::MFEM XSDK::SUNDIALS XSDK::HYPRE Ginkgo::ginkgo PETSC::ALL ${METIS_LIBRARY} ZLIB::ZLIB)
