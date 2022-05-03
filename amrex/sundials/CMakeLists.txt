add_executable(amrex_sundials_advection_diffusion
  amrex_sundials_advection_diffusion.cpp)
if(ENABLE_CUDA)
  set_source_files_properties(amrex_sundials_advection_diffusion.cpp
    PROPERTIES LANGUAGE CUDA)
  set_target_properties(amrex_sundials_advection_diffusion
    PROPERTIES CUDA_SEPARABLE_COMPILATION ON) # This adds -dc
endif()
target_link_libraries(amrex_sundials_advection_diffusion
  XSDK::AMReX XSDK::SUNDIALS)
install(TARGETS amrex_sundials_advection_diffusion RUNTIME DESTINATION bin)
