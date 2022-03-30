cmake_minimum_required(VERSION 3.12)
project(mfem-examples
        DESCRIPTION "MFEM Examples"
        LANGUAGES CXX)

if(ENABLE_GINKGO)
  add_subdirectory(ginkgo)
endif()
if(ENABLE_PETSC AND ENABLE_SUPERLU AND ENABLE_HYPRE)
  add_subdirectory(hypre-superlu)
endif()
if(ENABLE_PETSC)
  add_subdirectory(petsc)
endif()
if(ENABLE_PETSC AND ENABLE_SUNDIALS)
  add_subdirectory(sundials)
endif()

install(FILES data/star.mesh data/beam-hex.mesh DESTINATION data)
