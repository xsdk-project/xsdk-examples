# enforce OpenMP requirement; OpenMP search results should be cached at
# this point
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

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

xsdk_add_test(
    NAME
    PLASMA-ex1solve
    COMMAND
    env
    OMP_NUM_THREADS=1
    $<TARGET_FILE:ex1solve>
    --n=1000
)

install(TARGETS ex1solve RUNTIME DESTINATION bin)
