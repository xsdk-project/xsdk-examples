cmake_minimum_required(VERSION 3.8)
project(mfem-gko LANGUAGES CXX)

enable_language(CUDA)

set(MFEM_DIR "/path/to/mfem")
find_package(MFEM 4.0 REQUIRED
  HINTS ${MFEM_DIR})

set(Ginkgo_DIR "/path/to/ginkgo/")
find_package(Ginkgo REQUIRED)

add_executable(${PROJECT_NAME} mfem_ex1_gko.cpp)
target_sources(${PROJECT_NAME} PRIVATE mfem_ex1_gko.cpp mfem_wrapper.cpp)
target_compile_options(${PROJECT_NAME} PRIVATE -g -Wall -fopenmp)

target_link_libraries(${PROJECT_NAME}
         PRIVATE ${MFEM_LIBRARIES} ${GINKGO_LIB_DIR} -fopenmp)

target_include_directories(${PROJECT_NAME}
         PRIVATE ${GINKGO_INC_DIR} ${GINKGO_LIB_DIR} 
                ${MFEM_INCLUDE_DIRS} ${MFEM_LIB_DIR})
