cmake_minimum_required(VERSION 3.17)
project(
  strumpack-butterflypack
  DESCRIPTION "STRUMPACK + ButterflyPACK Example"
  LANGUAGES CXX)

# set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

add_executable(sparse sparse.cpp)

target_link_libraries(sparse PRIVATE STRUMPACK::strumpack)
if(MATH_LIBRARY)
  target_link_libraries(sparse PRIVATE ${MATH_LIBRARY})
endif()

install(TARGETS sparse RUNTIME DESTINATION bin)
