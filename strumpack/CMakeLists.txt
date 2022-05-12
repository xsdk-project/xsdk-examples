add_executable(sparse sparse.cpp)

target_link_libraries(sparse PRIVATE STRUMPACK::strumpack)
if(MATH_LIBRARY)
  target_link_libraries(sparse PRIVATE ${MATH_LIBRARY})
endif()

install(TARGETS sparse RUNTIME DESTINATION bin)
