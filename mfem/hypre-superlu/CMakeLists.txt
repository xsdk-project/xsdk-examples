if (ENABLE_CUDA)
    set_source_files_properties(convdiff.cpp PROPERTIES LANGUAGE CUDA)
endif()
if (ENABLE_HIP)
    set_source_files_properties(convdiff.cpp PROPERTIES LANGUAGE HIP)
endif()
add_executable(convdiff convdiff.cpp)
target_link_libraries(convdiff PRIVATE ${MFEM_ALL_DEPS})

xsdk_add_test(
    NAME
    MFEM-convdiff--hypre-boomeramg
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:convdiff>
    --no-superlu
    --velocity
    100.0
)
xsdk_add_test(
    NAME
    MFEM-convdiff--superlu
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:convdiff>
    --superlu
    --velocity
    100.0
)

install(TARGETS convdiff RUNTIME DESTINATION bin)
