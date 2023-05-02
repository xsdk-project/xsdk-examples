# SUNDIALS Copyright Start Copyright (c) 2002-2019, Lawrence Livermore National Security and
# Southern Methodist University. All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause SUNDIALS Copyright End
# ---------------------------------------------------------------

# add the executables
if(ENABLE_PETSC)
    add_executable(cv_petsc_ex7 cv_petsc_ex7.c)
    target_link_libraries(
        cv_petsc_ex7 PRIVATE XSDK::SUNDIALS XSDK::PETSc MPI::MPI_C
    )
    if(MATH_LIBRARY)
        target_link_libraries(cv_petsc_ex7 PRIVATE ${MATH_LIBRARY})
    endif()
    if (TARGET OpenMP::OpenMP_C)
      target_link_libraries(cv_petsc_ex7 PRIVATE OpenMP::OpenMP_C)
    endif()
    xsdk_add_test(NAME SUNDIALS-cv_petsc_ex7_1 COMMAND $<TARGET_FILE:cv_petsc_ex7> MPI_NPROCS 4)
    xsdk_add_test(
        NAME
        SUNDIALS-cv_petsc_ex7_2
        COMMAND
        $<TARGET_FILE:cv_petsc_ex7>
        --snes_type
        anderson
        MPI_NPROCS
        4
    )
    install(TARGETS cv_petsc_ex7 RUNTIME DESTINATION bin)
endif()
if(ENABLE_SUPERLU)
    set(tgt ark_brusselator1D_FEM_sludist)
    add_executable(${tgt} ark_brusselator1D_FEM_sludist.cpp)
    target_link_libraries(${tgt} PRIVATE XSDK::SUNDIALS XSDK::SUPERLU MPI::MPI_C)
    if(MATH_LIBRARY)
        target_link_libraries(ark_brusselator1D_FEM_sludist PRIVATE ${MATH_LIBRARY})
    endif()
    xsdk_add_test(NAME SUNDIALS-${tgt} COMMAND $<TARGET_FILE:${tgt}> MPI_NPROCS 1)
    install(TARGETS ${tgt} RUNTIME DESTINATION bin)
endif()
if(ENABLE_MAGMA)
    if(ENABLE_CUDA)
      set_source_files_properties(cv_bruss_batched_magma.cpp PROPERTIES LANGUAGE CUDA)
    elseif(ENABLE_HIP)
      set_source_files_properties(cv_bruss_batched_magma.cpp PROPERTIES LANGUAGE HIP)
    endif()
    add_executable(cv_bruss_batched_magma cv_bruss_batched_magma.cpp)
    target_link_libraries(cv_bruss_batched_magma PRIVATE XSDK::SUNDIALS XSDK::MAGMA)
    if(MATH_LIBRARY)
        target_link_libraries(cv_bruss_batched_magma PRIVATE ${MATH_LIBRARY})
    endif()
    xsdk_add_test(
        NAME SUNDIALS-cv_bruss_batched_magma COMMAND $<TARGET_FILE:cv_bruss_batched_magma>
    )
    install(TARGETS cv_bruss_batched_magma RUNTIME DESTINATION bin)
endif()
