add_executable(obstacle obstacle.cpp)

target_link_libraries(
    obstacle PRIVATE XSDK::MFEM XSDK::HYPRE XSDK::PETSc XSDK::SUNDIALS ZLIB::ZLIB MPI::MPI_CXX
)

xsdk_add_test(
    NAME
    MFEM-obstacle
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:obstacle>
    --order
    2
)

install(TARGETS obstacle RUNTIME DESTINATION bin)
