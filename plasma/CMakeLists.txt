cmake_minimum_required(VERSION 3.8)
project(
  plasmapp
  DESCRIPTION "Example of PLASMA with SLATE, LAPACK++, and BLAS++"
  LANGUAGES CXX
)

message(STATUS "p ${PLASMA_DIR}/include s ${SLATE_DIR}/include l ${LAPACKPP_DIR}/include b ${BLASPP_DIR}/include")
add_executable(ex1solve ex1solve.cpp)
target_include_directories(ex1solve PRIVATE ${PLASMA_DIR}/include ${SLATE_DIR}/include ${LAPACKPP_DIR}/include ${BLASPP_DIR}/include)
if (NOT ENABLE_CUDA)
    target_compile_definitions(ex1solve PRIVATE SLATE_NO_CUDA)
endif()
target_compile_definitions(ex1solve PRIVATE SLATE_NO_HIP)

target_link_libraries(ex1solve PRIVATE XSDK::PLASMA XSDK::SLATE XSDK::LAPACKPP XSDK::BLASPP)
if(MATH_LIBRARY)
  target_link_libraries(ex1solve PRIVATE ${MATH_LIBRARY})
endif()

install(TARGETS ex1solve RUNTIME DESTINATION bin)
