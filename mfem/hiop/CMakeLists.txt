if (ENABLE_HIP)
    set_source_files_properties(adv.cpp PROPERTIES LANGUAGE HIP)
endif()
add_executable(adv adv.cpp)
target_link_libraries(adv PRIVATE ${MFEM_ALL_DEPS})
if (NOT ENABLE_CUDA AND NOT ENABLE_HIP)
    # This test does not seem to work with HYPRE built with CUDA or HIP.
    xsdk_add_test(
        NAME
        MFEM-HIOP-adv
        MPI_NPROCS
        4
        COMMAND
        $<TARGET_FILE:adv>
        -no-vis
    )
endif()
install(TARGETS adv RUNTIME DESTINATION bin)
