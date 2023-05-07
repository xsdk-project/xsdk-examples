# All MFEM dependencies
set(MFEM_ALL_DEPS XSDK::MFEM)

# Ginkgo is included at the upper level when MFEM is enabled
if(ENABLE_GINKGO)
    add_subdirectory(ginkgo)
endif()
if(ENABLE_HIOP)
    add_subdirectory(hiop)
endif()
if(ENABLE_HYPRE)
    add_subdirectory(hypre)
endif()
if(ENABLE_HYPRE AND ENABLE_STRUMPACK)
    add_subdirectory(strumpack)
endif()
if(ENABLE_HYPRE AND ENABLE_SUNDIALS)
    add_subdirectory(sundials)
endif()
if(ENABLE_HYPRE AND ENABLE_SUPERLU)
    add_subdirectory(hypre-superlu)
endif()
if(ENABLE_PETSC)
    add_subdirectory(petsc)
endif()
if(ENABLE_PUMI)
    add_subdirectory(pumi)
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

install(FILES
    data/amr-quad.mesh
    data/beam-hex.mesh
    data/disc-nurbs.mesh
    data/periodic-cube.mesh
    data/periodic-hexagon.mesh
    data/periodic-segment.mesh
    data/periodic-square.mesh
    data/star.mesh
    DESTINATION data)
# TODO: install PUMI files
