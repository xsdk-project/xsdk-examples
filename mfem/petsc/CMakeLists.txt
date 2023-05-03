add_executable(obstacle obstacle.cpp)

target_link_libraries(obstacle PRIVATE ${MFEM_ALL_DEPS})

# This is a serial example linked with MPI, so launch it with 'mpirun' or
# equivalent to avoid warnings on some platforms.
xsdk_add_test(
    NAME
    MFEM-obstacle
    MPI_NPROCS
    1
    COMMAND
    $<TARGET_FILE:obstacle>
    --order
    2
)

install(TARGETS obstacle RUNTIME DESTINATION bin)
