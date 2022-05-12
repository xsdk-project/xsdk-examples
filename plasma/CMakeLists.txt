add_executable(ex1solve ex1solve.cpp)
if(NOT ENABLE_CUDA)
    target_compile_definitions(ex1solve PRIVATE SLATE_NO_CUDA)
endif()
target_compile_definitions(ex1solve PRIVATE SLATE_NO_HIP)

target_link_libraries(
    ex1solve PRIVATE XSDK::PLASMA XSDK::SLATE XSDK::LAPACKPP XSDK::BLASPP MPI::MPI_CXX
)
if(MATH_LIBRARY)
    target_link_libraries(ex1solve PRIVATE ${MATH_LIBRARY})
endif()
if(ENABLE_CUDA)
    target_link_libraries(ex1solve PRIVATE CUDA::cudart)
endif()

install(TARGETS ex1solve RUNTIME DESTINATION bin)
