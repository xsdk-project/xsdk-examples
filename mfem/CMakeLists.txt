cmake_minimum_required(VERSION 3.12)
project(mfem-examples
        DESCRIPTION "MFEM Examples"
        LANGUAGES CXX)

add_subdirectory(hypre-superlu)
add_subdirectory(ginkgo)
add_subdirectory(petsc)
add_subdirectory(sundials)
