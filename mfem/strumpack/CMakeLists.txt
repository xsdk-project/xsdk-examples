add_executable(diffusion-eigen diffusion-eigen.cpp)

target_link_libraries(
    diffusion-eigen
    PRIVATE XSDK::MFEM
            XSDK::PETSc
            XSDK::HYPRE
            XSDK::SUPERLU
            XSDK::SUNDIALS
            STRUMPACK::strumpack
            ZLIB::ZLIB
)

xsdk_add_test(
    NAME
    MFEM-diffusion-eigen-no-strumpack
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:diffusion-eigen>
    -m
    ../data/star.mesh
    --no-strumpack
)
xsdk_add_test(
    NAME
    MFEM-diffusion-eigen-strumpack
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:diffusion-eigen>
    -m
    ../data/star.mesh
    --strumpack
)
xsdk_add_test(
    NAME
    MFEM-diffusion-eigen-superlu
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:diffusion-eigen>
    -m
    ../data/star.mesh
    --superlu
)

install(TARGETS diffusion-eigen RUNTIME DESTINATION bin)
