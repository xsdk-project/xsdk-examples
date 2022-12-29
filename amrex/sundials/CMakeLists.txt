if(ENABLE_CUDA)
    set_source_files_properties(amrex_sundials_advection_diffusion.cpp PROPERTIES LANGUAGE CUDA)
    set_target_properties(
        amrex_sundials_advection_diffusion PROPERTIES CUDA_SEPARABLE_COMPILATION ON
    ) # This adds -dc
endif()

add_executable(amrex_sundials_advection_diffusion amrex_sundials_advection_diffusion.cpp)

target_link_libraries(amrex_sundials_advection_diffusion PRIVATE XSDK::AMReX XSDK::SUNDIALS MPI::MPI_CXX)

xsdk_add_test(
    NAME AMREX-amrex_sundials_advection_diffusion COMMAND
    $<TARGET_FILE:amrex_sundials_advection_diffusion>
)

install(TARGETS amrex_sundials_advection_diffusion RUNTIME DESTINATION bin)
