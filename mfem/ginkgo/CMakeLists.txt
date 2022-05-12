add_executable(mfem_ex1_gko mfem_ex1_gko.cpp)
target_sources(mfem_ex1_gko PRIVATE mfem_wrapper.cpp)

target_link_libraries(
    mfem_ex1_gko
    PRIVATE XSDK::MFEM
            XSDK::SUNDIALS
            XSDK::HYPRE
            XSDK::PETSc
            Ginkgo::ginkgo
            ZLIB::ZLIB
            MPI::MPI_CXX
)

install(TARGETS mfem_ex1_gko RUNTIME DESTINATION bin)
