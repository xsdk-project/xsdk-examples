add_executable(heffte_example heffte_example_gpu.cpp)
target_link_libraries(heffte_example Heffte::Heffte)

xsdk_add_test(NAME heFFTe-MAGMA MPI_NPROCS 2 COMMAND $<TARGET_FILE:heffte_example>)
