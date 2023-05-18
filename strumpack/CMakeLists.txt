add_executable(sparse sparse.cpp)

target_link_libraries(sparse PRIVATE STRUMPACK::strumpack)
if(MATH_LIBRARY)
    target_link_libraries(sparse PRIVATE ${MATH_LIBRARY})
endif()

if (NOT ENABLE_CUDA)
    # This example randomly fails with strumpack@7.0.1 and slate@2022.07.00,
    # i.e. xsdk@0.8.0, with CUDA enabled, so it is disabled for now.
    # See https://github.com/xsdk-project/xsdk-examples/issues/47
    xsdk_add_test(
        NAME
        STRUMPACK-sparse
        MPI_NPROCS
        4
        COMMAND
        $<TARGET_FILE:sparse>
        24
        --sp_compression
        hodlr
        --hodlr_butterfly_levels
        10
        ENVIRONMENT
        OMP_NUM_THREADS=4
    )
endif()

install(TARGETS sparse RUNTIME DESTINATION bin)
