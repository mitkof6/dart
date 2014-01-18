###############################################
# apps/performance
file(GLOB performance_srcs "*.cpp")
file(GLOB performance_hdrs "*.h")
add_executable(performance ${performance_srcs} ${performance_hdrs})
target_link_libraries(performance dart)
set_target_properties(performance PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
