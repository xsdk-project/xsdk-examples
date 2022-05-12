add_executable(transient-heat transient-heat.cpp)
target_link_libraries(
    transient-heat PRIVATE XSDK::SUNDIALS XSDK::PETSc XSDK::HYPRE XSDK::MFEM ZLIB::ZLIB
                           MPI::MPI_CXX
)
install(TARGETS transient-heat RUNTIME DESTINATION bin)

add_executable(advection advection.cpp)
target_link_libraries(
    advection PRIVATE XSDK::PETSc XSDK::SUNDIALS XSDK::HYPRE XSDK::MFEM ZLIB::ZLIB MPI::MPI_CXX
)
install(TARGETS advection RUNTIME DESTINATION bin)
