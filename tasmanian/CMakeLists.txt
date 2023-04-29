find_package(Tasmanian 7.5 PATHS "${TASMANIAN_DIR}" REQUIRED)

add_executable(tasmanian_magma example_unstructured_grid.cpp)
target_link_libraries(tasmanian_magma Tasmanian::Tasmanian)

xsdk_add_test(tasmanian tasmanian_magma)
