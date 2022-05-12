set(example magnetic-diffusion)

add_executable(${example} ${example}.cpp)
target_link_libraries(${example} PRIVATE
    XSDK::MFEM XSDK::PETSc XSDK::HYPRE XSDK::SUPERLU XSDK::SUNDIALS ZLIB::ZLIB MPI::MPI_CXX)

install(TARGETS ${example} RUNTIME DESTINATION bin)
