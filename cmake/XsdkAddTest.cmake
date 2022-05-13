macro(xsdk_add_test)

    set(options)

    # MPI_NPROCS = number of mpi tasks to use in parallel tests
    set(oneValueArgs "NAME" "MPI_NPROCS")

    set(multiValueArgs "COMMAND" "ENVIRONMENT")

    # parse inputs and create variables SUNDIALS_ADD_TEST_<keyword>
    cmake_parse_arguments(xsdk_add_test "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(xsdk_add_test_MPI_NPROCS)
        add_test(NAME ${xsdk_add_test_NAME}
                 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${xsdk_add_test_MPI_NPROCS}
                         ${MPIEXEC_PREFLAGS} ${xsdk_add_test_COMMAND} ${MPIEXEC_POSTFLAGS}
        )
    else()
        add_test(NAME ${xsdk_add_test_NAME} COMMAND ${xsdk_add_test_COMMAND})
    endif()

    if(xsdk_add_test_ENVIRONMENT)
        set_tests_properties(
            ${xsdk_add_test_NAME} PROPERTIES ENVIRONMENT ${xsdk_add_test_ENVIRONMENT}
        )
    endif()

endmacro(xsdk_add_test)
