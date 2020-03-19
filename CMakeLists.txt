cmake_minimum_required(VERSION 3.10)
project( xsdk-examples
         DESCRIPTION "xSDK Examples"
         LANGUAGES CXX C)

add_subdirectory(hypre)
add_subdirectory(petsc)
add_subdirectory(sundials)
add_subdirectory(trilinos)
