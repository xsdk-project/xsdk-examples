add_executable(sparse sparse.cpp)

target_link_libraries(sparse PRIVATE STRUMPACK::strumpack)
if(MATH_LIBRARY)
    target_link_libraries(sparse PRIVATE ${MATH_LIBRARY})
endif()

xsdk_add_test(
    NAME
    STRUMPACK-sparse
    MPI_NPROCS
    4
    COMMAND
    $<TARGET_FILE:sparse>
    24
    --sp_compression
    hodlr
    --hodlr_butterfly_levels
    10
    ENVIRONMENT
    OMP_NUM_THREADS=4
)

install(TARGETS sparse RUNTIME DESTINATION bin)
