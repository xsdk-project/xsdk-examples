find_package(Heffte 2.2 PATHS "${HEFFTE_DIR}" REQUIRED)

add_executable(heffte_example heffte_example_gpu.cpp)
target_link_libraries(heffte_example Heffte::Heffte)

xsdk_add_test(NAME heFFTe-MAGMA MPI_NPROCS 2 COMMAND heffte_example)
