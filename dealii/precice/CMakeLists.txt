PROJECT("dealii-precice")

CMAKE_MINIMUM_REQUIRED(VERSION 3.3.0)

FIND_PACKAGE(deal.II 9.2.0 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )

FIND_PACKAGE(precice REQUIRED
  HINTS ${PRECICE_DIR}
  )

DEAL_II_INITIALIZE_CACHED_VARIABLES()

# build Laplace problem, which requires both deal.II and preCICE
ADD_EXECUTABLE(laplace_problem laplace_problem.cc)
DEAL_II_SETUP_TARGET(laplace_problem)
TARGET_LINK_LIBRARIES(laplace_problem precice::precice)

# build the boundary condition, which only requires preCICE
ADD_EXECUTABLE(boundary_condition boundary_condition.cc)
TARGET_LINK_LIBRARIES(boundary_condition precice::precice)

include(CTest)
add_test(
    NAME 
    ${CMAKE_PROJECT_NAME} 
    COMMAND
    runtargets.sh
)

#xsdk_add_test(
#    NAME
#    ${CMAKE_PROJECT_NAME}
#    MPI_NPROCS
#    1
#    COMMAND
#    runtargets.sh
#)
