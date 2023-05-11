if(ENABLE_CUDA)
    set(CMAKE_CUDA_HOST_LINK_LAUNCHER ${CMAKE_CXX_COMPILER})
    set_source_files_properties(amrex_sundials_advection_diffusion.cpp PROPERTIES LANGUAGE CUDA)
endif()
if (ENABLE_HIP)
    set_source_files_properties(amrex_sundials_advection_diffusion.cpp PROPERTIES LANGUAGE HIP)
endif()

add_executable(amrex_sundials_advection_diffusion amrex_sundials_advection_diffusion.cpp)

target_link_libraries(amrex_sundials_advection_diffusion PRIVATE XSDK::AMReX XSDK::SUNDIALS MPI::MPI_CXX)
if(ENABLE_CUDA)
    set_target_properties(amrex_sundials_advection_diffusion
        PROPERTIES CUDA_SEPARABLE_COMPILATION ON
    ) # This adds -dc
endif()

xsdk_add_test(
    NAME AMREX-amrex_sundials_advection_diffusion COMMAND
    $<TARGET_FILE:amrex_sundials_advection_diffusion>
)

install(TARGETS amrex_sundials_advection_diffusion RUNTIME DESTINATION bin)
