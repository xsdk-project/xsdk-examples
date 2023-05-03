add_executable(mfem_ex1_gko mfem_ex1_gko.cpp)
target_sources(mfem_ex1_gko PRIVATE mfem_wrapper.cpp)
target_link_libraries(mfem_ex1_gko PRIVATE ${MFEM_ALL_DEPS})

xsdk_add_test(NAME MFEM-mfem_ex1_gko COMMAND $<TARGET_FILE:mfem_ex1_gko> --no-partial-assembly)
if(ENABLE_CUDA)
    xsdk_add_test(
        NAME MFEM-mfem_ex1_gko--gpu
        COMMAND $<TARGET_FILE:mfem_ex1_gko>
        --partial-assembly --device cuda)
endif()

install(TARGETS mfem_ex1_gko RUNTIME DESTINATION bin)
