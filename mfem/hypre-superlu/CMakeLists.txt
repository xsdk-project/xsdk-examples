add_executable(convdiff convdiff.cpp)
target_link_libraries(convdiff PRIVATE XSDK::MFEM XSDK::PETSc XSDK::HYPRE
    XSDK::SUPERLU XSDK::SUNDIALS ZLIB::ZLIB MPI::MPI_CXX)

install(TARGETS convdiff RUNTIME DESTINATION bin)
