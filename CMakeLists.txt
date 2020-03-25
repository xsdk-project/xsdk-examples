cmake_minimum_required(VERSION 3.10)
project( xsdk-examples
         DESCRIPTION "xSDK Examples"
         LANGUAGES CXX C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

add_subdirectory(hypre)
add_subdirectory(mfem)
add_subdirectory(petsc)
add_subdirectory(sundials)
add_subdirectory(trilinos)
