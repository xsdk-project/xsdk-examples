set(example magnetic-diffusion)

add_executable(${example} ${example}.cpp)
target_link_libraries(
    ${example}
    PRIVATE XSDK::MFEM
            XSDK::PETSc
            XSDK::HYPRE
            XSDK::SUPERLU
            XSDK::SUNDIALS
            ZLIB::ZLIB
            MPI::MPI_CXX
)

if(ENABLE_CUDA)
    xsdk_add_test(
        NAME
        MFEM-magnetic-diffusion-gpu
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
elseif(NOT ENABLE_CUDA)
    xsdk_add_test(
        NAME
        MFEM-magnetic-diffusion-cpu
        MPI_NPROCS
        4
        COMMAND
        $<TARGET_FILE:magnetic-diffusion>
        -m
        ../data/star.mesh
    )
endif()

install(TARGETS ${example} RUNTIME DESTINATION bin)
