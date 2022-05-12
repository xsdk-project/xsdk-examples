add_executable(obstacle obstacle.cpp)
target_link_libraries(obstacle PRIVATE XSDK::MFEM XSDK::HYPRE XSDK::PETSc
    XSDK::SUNDIALS ZLIB::ZLIB MPI::MPI_CXX)
install(TARGETS obstacle RUNTIME DESTINATION bin)
