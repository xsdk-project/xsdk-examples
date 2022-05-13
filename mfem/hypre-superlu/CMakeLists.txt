add_executable(convdiff convdiff.cpp)
target_link_libraries(
    convdiff
    PRIVATE XSDK::MFEM
            XSDK::PETSc
            XSDK::HYPRE
            XSDK::SUPERLU
            XSDK::SUNDIALS
            ZLIB::ZLIB
            MPI::MPI_CXX
)

xsdk_add_test(
    NAME
    MFEM-convdiff-no-superlu
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:convdiff>
    --no-superlu
    --velocity
    100.0
)
xsdk_add_test(
    NAME
    MFEM-convdiff-superlu
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:convdiff>
    --superlu
    --velocity
    100.0
)

install(TARGETS convdiff RUNTIME DESTINATION bin)
