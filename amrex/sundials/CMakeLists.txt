if(ENABLE_CUDA)
    set(CMAKE_CUDA_HOST_LINK_LAUNCHER ${CMAKE_CXX_COMPILER})
    set_source_files_properties(amrex_sundials_advection_diffusion.cpp PROPERTIES LANGUAGE CUDA)
endif()
# Due to some linking issues when HIP is enabled, we skip this example in that
# case; for details, see the following GitHub issue and PR:
#     https://github.com/xsdk-project/xsdk-examples/issues/50
#     https://github.com/xsdk-project/xsdk-examples/pull/52
if (ENABLE_HIP)
    return()
    # set_source_files_properties(amrex_sundials_advection_diffusion.cpp PROPERTIES LANGUAGE HIP)
endif()

add_executable(amrex_sundials_advection_diffusion amrex_sundials_advection_diffusion.cpp)

target_link_libraries(amrex_sundials_advection_diffusion PRIVATE XSDK::AMReX XSDK::SUNDIALS MPI::MPI_CXX)
if(ENABLE_CUDA)
    set_target_properties(amrex_sundials_advection_diffusion
        PROPERTIES CUDA_SEPARABLE_COMPILATION ON
    ) # This adds -dc
endif()

xsdk_add_test(
    NAME AMREX-amrex_sundials_advection_diffusion
    MPI_NPROCS 1
    COMMAND $<TARGET_FILE:amrex_sundials_advection_diffusion>
)

install(TARGETS amrex_sundials_advection_diffusion RUNTIME DESTINATION bin)
