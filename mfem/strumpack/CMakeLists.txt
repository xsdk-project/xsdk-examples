if (ENABLE_CUDA)
    set_source_files_properties(diffusion-eigen.cpp PROPERTIES LANGUAGE CUDA)
endif()
if (ENABLE_HIP)
    set_source_files_properties(diffusion-eigen.cpp PROPERTIES LANGUAGE HIP)
endif()
add_executable(diffusion-eigen diffusion-eigen.cpp)
target_link_libraries(diffusion-eigen PRIVATE ${MFEM_ALL_DEPS})

xsdk_add_test(
    NAME
    MFEM-diffusion-eigen--strumpack
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:diffusion-eigen>
    -m
    ../data/star.mesh
    -ls
    strumpack
)
if(ENABLE_SUPERLU)
    xsdk_add_test(
        NAME
        MFEM-diffusion-eigen--superlu
        MPI_NPROCS
        4
        COMMAND
        $<TARGET_FILE:diffusion-eigen>
        -m
        ../data/star.mesh
        -ls
        superlu
    )
endif()
xsdk_add_test(
    NAME
    MFEM-diffusion-eigen--hypre-boomeramg
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:diffusion-eigen>
    -m
    ../data/star.mesh
    -ls
    hypre
)

install(TARGETS diffusion-eigen RUNTIME DESTINATION bin)
