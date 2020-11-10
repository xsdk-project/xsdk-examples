cmake_minimum_required(VERSION 3.10)
project(mfem-examples
        DESCRIPTION "MFEM Examples"
        LANGUAGES CXX)

#add_subdirectory(ginkgo)
add_subdirectory(hypre-superlu)
add_subdirectory(petsc)
add_subdirectory(sundials)
