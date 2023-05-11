
# At least in some cases, deal.II needs to use the MPI wrapper:
if (MPI_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

# build Laplace problem, which requires both deal.II and preCICE
add_executable(laplace_problem laplace_problem.cc)
DEAL_II_SETUP_TARGET(laplace_problem)
target_link_libraries(laplace_problem precice::precice)

# build the boundary condition, which only requires preCICE
add_executable(boundary_condition boundary_condition.cc)
target_link_libraries(boundary_condition precice::precice)

# Copy the shell script 'runtargets.sh' and config file 'precice-config.xml'
# from the source tree to the build tree so that the example can be run there.
add_custom_command(
    OUTPUT runtargets.sh precice-config.xml
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/runtargets.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/precice-config.xml
        ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Copying the deal.II + preCICE run script and config ..."
)
add_custom_target(copy_precise_data ALL
    DEPENDS runtargets.sh precice-config.xml
)

xsdk_add_test(
    NAME 
    "deal.II-preCICE-laplace_problem"
    COMMAND
    runtargets.sh
)
