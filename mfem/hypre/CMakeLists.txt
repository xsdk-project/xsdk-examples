set(example magnetic-diffusion)

add_executable(${example} ${example}.cpp)
target_link_libraries(
    ${example}
    PRIVATE XSDK::MFEM
            XSDK::HYPRE
            XSDK::SUNDIALS
            ZLIB::ZLIB
            MPI::MPI_CXX
)

xsdk_add_test(
    NAME
    MFEM-magnetic-diffusion--cpu
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:magnetic-diffusion>
    -m
    ../data/star.mesh
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
    )
    # TODO: add the other GPU sample runs from README.md
endif()

install(TARGETS ${example} RUNTIME DESTINATION bin)
