# Ginkgo is included at the upper level when MFEM is enabled
if(ENABLE_GINKGO)
    add_subdirectory(ginkgo)
endif()
if(ENABLE_HYPRE)
    add_subdirectory(hypre)
endif()
if(ENABLE_SUPERLU AND ENABLE_HYPRE)
    add_subdirectory(hypre-superlu)
endif()
if(ENABLE_PETSC)
    add_subdirectory(petsc)
endif()
if(ENABLE_SUNDIALS AND ENABLE_HYPRE)
    add_subdirectory(sundials)
endif()
if(ENABLE_STRUMPACK AND ENABLE_HYPRE)
    add_subdirectory(strumpack)
endif()

# Copy the data directory from the source tree to the build tree so that the examples can be run
# from in the build tree.
add_custom_command(
    OUTPUT data_is_copied
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data data
    COMMAND ${CMAKE_COMMAND} -E touch data_is_copied
    COMMENT "Copying the MFEM data directory ..."
)
add_custom_target(copy_mfem_data ALL DEPENDS data_is_copied)

install(FILES data/star.mesh data/beam-hex.mesh data/periodic-hexagon.mesh DESTINATION data)
