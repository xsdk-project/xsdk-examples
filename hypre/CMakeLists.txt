if(ENABLE_SUPERLU AND (NOT ENABLE_CUDA))
    add_executable(ij_laplacian ij_laplacian.c)

    target_link_libraries(ij_laplacian PRIVATE XSDK::HYPRE XSDK::SUPERLU MPI::MPI_C)

    if(MATH_LIBRARY)
        target_link_libraries(ij_laplacian PRIVATE ${MATH_LIBRARY})
    endif()

    xsdk_add_test(
        NAME
        HYPRE-ij_laplacian
        MPI_NPROCS
        4
        COMMAND
        $<TARGET_FILE:ij_laplacian>
        -dslu_th
        50
    )

    install(TARGETS ij_laplacian RUNTIME DESTINATION bin)

    # Copy glvis helper script and directory
    file(COPY vis DESTINATION .)
    install(
        DIRECTORY vis
        DESTINATION share/xsdk-examples/hypre
        PATTERN
            vis/*
            PERMISSIONS
                OWNER_EXECUTE
                OWNER_WRITE
                OWNER_READ
                GROUP_EXECUTE
                GROUP_READ
                WORLD_READ
                WORLD_EXECUTE
    )
else()
    message(
        STATUS
            "SKIPPED HYPRE-ij_laplacian example. Requires ENABLE_SUPERLU=ON and ENABLE_CUDA=OFF but got ${ENABLE_SUPERLU} and ${ENABLE_CUDA}"
    )
endif()
