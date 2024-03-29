add_executable(ex19 ex19.c)
target_link_libraries(ex19 PRIVATE XSDK::PETSc MPI::MPI_C)
if (TARGET OpenMP::OpenMP_C)
  target_link_libraries(ex19 PRIVATE OpenMP::OpenMP_C)
endif()

# The makefile to run examples requires rm command
find_program(RM rm)

# If PETSc was found without setting PETSc_DIR, e.g. using CMAKE_PREFIX_PATH
if(NOT PETSc_DIR)
    set(PETSc_DIR ${PC_PETSc_PREFIX})
endif()
# Handle cases like Flux where MPIEXEC_EXECUTABLE can be "flux;run"
string(REPLACE ";" " " MPIEXEC_EXECUTABLE "${MPIEXEC_EXECUTABLE}")

# We use the makefile to launch the tests
configure_file(makefile ${CMAKE_CURRENT_BINARY_DIR}/makefile-run)
file(COPY output DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/output)

if(MATH_LIBRARY)
    target_link_libraries(ex19 PRIVATE ${MATH_LIBRARY})
endif()

xsdk_add_test(
    NAME
    PETSc-ex19_1
    COMMAND
    make
    -f
    makefile-run
    runex19
    ENVIRONMENT
    PETSC_DIR=${PETSc_DIR}
)
if(ENABLE_CUDA)
    xsdk_add_test(
        NAME
        PETSc-ex19_cuda
        COMMAND
        make
        -f
        makefile-run
        runex19_cuda
        ENVIRONMENT
        PETSC_DIR=${PETSc_DIR}
    )
endif()
if(ENABLE_HIP)
    xsdk_add_test(
        NAME
        PETSc-ex19_hip
        COMMAND
        make
        -f
        makefile-run
        runex19_hip
        ENVIRONMENT
        PETSC_DIR=${PETSc_DIR}
    )
endif()
if(ENABLE_HYPRE)
    if(ENABLE_CUDA)
        xsdk_add_test(
            NAME
            PETSc-ex19_hypre_cuda
            COMMAND
            make
            -f
            makefile-run
            runex19_hypre_cuda
            ENVIRONMENT
            PETSC_DIR=${PETSc_DIR}
        )
    else()
        xsdk_add_test(
            NAME
            PETSc-ex19_hypre
            COMMAND
            make
            -f
            makefile-run
            runex19_hypre
            ENVIRONMENT
            PETSC_DIR=${PETSc_DIR}
        )
    endif()
endif()
if(ENABLE_SUPERLU)
    xsdk_add_test(
        NAME
        PETSc-ex19_superlu_dist
        COMMAND
        make
        -f
        makefile-run
        runex19_superlu_dist
        ENVIRONMENT
        PETSC_DIR=${PETSc_DIR}
    )
endif()

install(FILES output/ex19_1.testout output/ex19_cuda_1.out output/ex19_hypre.out
              output/ex19_hip.out output/ex19_superlu.out DESTINATION output
)
install(TARGETS ex19 RUNTIME DESTINATION bin)
