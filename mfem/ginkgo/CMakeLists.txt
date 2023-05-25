if (ENABLE_CUDA)
    set_source_files_properties(mfem_ex22_gko.cpp PROPERTIES LANGUAGE CUDA)
endif()
if (ENABLE_HIP)
    set_source_files_properties(mfem_ex22_gko.cpp PROPERTIES LANGUAGE HIP)
endif()
add_executable(mfem_ex22_gko mfem_ex22_gko.cpp)
target_link_libraries(mfem_ex22_gko PRIVATE ${MFEM_ALL_DEPS})

xsdk_add_test(NAME MFEM-mfem_ex22_gko COMMAND $<TARGET_FILE:mfem_ex22_gko> --no-partial-assembly)
if(ENABLE_CUDA)
    xsdk_add_test(
        NAME MFEM-mfem_ex22_gko--gpu
        COMMAND $<TARGET_FILE:mfem_ex22_gko>
        --partial-assembly --device cuda)
endif()
if(ENABLE_HIP)
    xsdk_add_test(
        NAME MFEM-mfem_ex22_gko--gpu
        COMMAND $<TARGET_FILE:mfem_ex22_gko>
        --partial-assembly --device hip)
endif()

install(TARGETS mfem_ex22_gko RUNTIME DESTINATION bin)
