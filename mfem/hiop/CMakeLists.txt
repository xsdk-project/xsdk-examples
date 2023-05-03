add_executable(adv adv.cpp)
target_link_libraries(adv PRIVATE ${MFEM_ALL_DEPS})
xsdk_add_test(
    NAME
    MFEM-HIOP-adv
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:adv>
    -no-vis
)
install(TARGETS adv RUNTIME DESTINATION bin)
