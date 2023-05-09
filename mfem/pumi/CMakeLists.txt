add_executable(adapt adapt.cpp)
target_link_libraries(adapt PRIVATE ${MFEM_ALL_DEPS})
# Temporarily disable the test until we have input PUMI mesh files
if (FALSE)
    xsdk_add_test(
        NAME
        MFEM-PUMI-adapt
        MPI_NPROCS
        8
        COMMAND
        $<TARGET_FILE:adapt>
        -no-vis
    )
endif()
install(TARGETS adapt RUNTIME DESTINATION bin)
