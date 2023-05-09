add_executable(adapt adapt.cpp)
target_link_libraries(adapt PRIVATE ${MFEM_ALL_DEPS})
xsdk_add_test(
    NAME
    MFEM-PUMI-adapt
    MPI_NPROCS
    2
    COMMAND
    $<TARGET_FILE:adapt>
    -no-vis
)
install(TARGETS adapt RUNTIME DESTINATION bin)
