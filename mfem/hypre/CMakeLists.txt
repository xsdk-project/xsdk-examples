set(example magnetic-diffusion)

add_executable(${example} ${example}.cpp)
target_link_libraries(${example} PRIVATE ${MFEM_ALL_DEPS})

xsdk_add_test(
    NAME
    MFEM-magnetic-diffusion--cpu
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:magnetic-diffusion>
    -m
    ../data/star.mesh
    -no-vis
)
# TODO: add the other CPU sample runs from README.md
if(ENABLE_CUDA)
    xsdk_add_test(
        NAME
        MFEM-magnetic-diffusion--gpu
        MPI_NPROCS
        4
        COMMAND
        $<TARGET_FILE:magnetic-diffusion>
        -m
        ../data/star.mesh
        -pa
        -d
        cuda
        -no-vis
    )
    # TODO: add the other GPU sample runs from README.md
endif()
if(ENABLE_HIP)
    xsdk_add_test(
        NAME
        MFEM-magnetic-diffusion--gpu
        MPI_NPROCS
        4
        COMMAND
        $<TARGET_FILE:magnetic-diffusion>
        -m
        ../data/star.mesh
        -pa
        -d
        hip
        -no-vis
    )
    # TODO: add the other GPU sample runs from README.md
endif()

install(TARGETS ${example} RUNTIME DESTINATION bin)
