add_executable(transient-heat transient-heat.cpp)
target_link_libraries(
    transient-heat PRIVATE XSDK::SUNDIALS XSDK::PETSc XSDK::HYPRE XSDK::MFEM ZLIB::ZLIB
                           MPI::MPI_CXX
)
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
)
install(TARGETS transient-heat RUNTIME DESTINATION bin)

add_executable(advection advection.cpp)
target_link_libraries(
    advection PRIVATE XSDK::PETSc XSDK::SUNDIALS XSDK::HYPRE XSDK::MFEM ZLIB::ZLIB MPI::MPI_CXX
)
if(ENABLE_CUDA)
    xsdk_add_test(
        NAME
        MFEM-advection
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
else()
    xsdk_add_test(
        NAME
        MFEM-advection
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
install(TARGETS advection RUNTIME DESTINATION bin)
