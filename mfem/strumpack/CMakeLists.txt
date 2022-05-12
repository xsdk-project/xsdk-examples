add_executable(diffusion-eigen diffusion-eigen.cpp)
target_link_libraries(diffusion-eigen PRIVATE
  XSDK::MFEM XSDK::PETSc XSDK::HYPRE XSDK::SUPERLU XSDK::SUNDIALS
  STRUMPACK::strumpack ZLIB::ZLIB)
install(TARGETS diffusion-eigen RUNTIME DESTINATION bin)
