##
#  CMake script for the sundials tutorial program:
##

# Set the name of the project and target:
set(TARGET "sundials")

# Declare all source files the target consists of. Here, this is only
# the one step-X.cc file, but as you expand your project you may wish
# to add other source files as well. If your project becomes much larger,
# you may want to either replace the following statement by something like
#  file(GLOB_RECURSE TARGET_SRC  "source/*.cc")
#  file(GLOB_RECURSE TARGET_INC  "include/*.h")
#  set(TARGET_SRC ${TARGET_SRC}  ${TARGET_INC})
# or switch altogether to the large project CMakeLists.txt file discussed
# in the "CMake in user projects" page accessible from the "User info"
# page of the documentation.
set(TARGET_SRC
  ${TARGET}.cc
  )

# Usually, you will not need to modify anything beyond this point...

cmake_minimum_required(VERSION 3.13.4)

if(NOT DEAL_II_WITH_SUNDIALS) # keep in one line
  message(FATAL_ERROR "
Error! This tutorial requires a deal.II library that was configured with the following options:
    DEAL_II_WITH_SUNDIALS = ON
However, the deal.II library found at ${DEAL_II_PATH} was configured with these options:
    DEAL_II_WITH_SUNDIALS = ${DEAL_II_WITH_SUNDIALS}
This conflicts with the requirements."
    )
endif()

deal_ii_initialize_cached_variables()
project(${TARGET})
add_executable(${TARGET} ${TARGET_SRC})
deal_ii_setup_target(${TARGET})

xsdk_add_test(
    NAME
    dealii-${TARGET}
    MPI_NPROCS
    1
    COMMAND
    $<TARGET_FILE:${TARGET}>
)
