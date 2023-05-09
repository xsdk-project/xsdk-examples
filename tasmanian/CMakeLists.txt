add_executable(tasmanian_magma example_unstructured_grid.cpp)
target_link_libraries(tasmanian_magma Tasmanian::Tasmanian)

xsdk_add_test(NAME tasmanian COMMAND tasmanian_magma)
