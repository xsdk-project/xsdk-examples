if (ENABLE_CUDA)
    set_source_files_properties(transient-heat.cpp PROPERTIES LANGUAGE CUDA)
endif()
if (ENABLE_HIP)
    set_source_files_properties(transient-heat.cpp PROPERTIES LANGUAGE HIP)
endif()
add_executable(transient-heat transient-heat.cpp)
target_link_libraries(transient-heat PRIVATE ${MFEM_ALL_DEPS})
if(NOT ENABLE_CUDA AND NOT ENABLE_HIP)
    # This example fails (as of 2023/05/08) when HYPRE is built with CUDA or
    # HIP, so it is disabled in that case for now.
    xsdk_add_test(
        NAME
        MFEM-transient-heat
        MPI_NPROCS
        4
        COMMAND
        $<TARGET_FILE:transient-heat>
        --mesh
        ../data/star.mesh
        --kappa
        0.5
        --alpha
        0.01
        --ode-solver
        8
        --no-visualization
    )
endif()
install(TARGETS transient-heat RUNTIME DESTINATION bin)

if (ENABLE_CUDA)
    set_source_files_properties(advection.cpp PROPERTIES LANGUAGE CUDA)
endif()
if (ENABLE_HIP)
    set_source_files_properties(advection.cpp PROPERTIES LANGUAGE HIP)
endif()
add_executable(advection advection.cpp)
target_link_libraries(advection PRIVATE ${MFEM_ALL_DEPS})
if(NOT ENABLE_CUDA AND NOT ENABLE_HIP)
    # This example fails (as of 2023/05/08) when CUDA or HIP is enabled, so we
    # disable it here for now.
    xsdk_add_test(
        NAME
        MFEM-advection--cpu
        MPI_NPROCS
        4
        COMMAND
        $<TARGET_FILE:advection>
        --device
        cpu
        --partial-assembly
        --ode-solver
        8
        --no-visualization
    )
endif()
if(ENABLE_CUDA)
    xsdk_add_test(
        NAME
        MFEM-advection--gpu
        MPI_NPROCS
        4
        COMMAND
        $<TARGET_FILE:advection>
        --device
        cuda
        --partial-assembly
        --ode-solver
        8
        --no-visualization
    )
endif()
if (ENABLE_HIP)
    # MFEM+SUNDIALS+HIP requires MFEM PR #3596 which as of 2023/05/08 is under
    # review.
endif()
install(TARGETS advection RUNTIME DESTINATION bin)
