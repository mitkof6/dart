###############################################
# apps/urdfLoader
file(GLOB urdfLoader_srcs "*.cpp")
file(GLOB urdfLoader_hdrs "*.h")
add_executable(urdfLoader ${urdfLoader_srcs} ${urdfLoader_hdrs})
target_link_libraries(urdfLoader dart)
set_target_properties(urdfLoader PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
