###############################################
# apps/hoop
file(GLOB hoop_srcs "*.cpp")
file(GLOB hoop_hdrs "*.h")
add_executable(hoop ${hoop_srcs} ${hoop_hdrs})
target_link_libraries(hoop dart)
set_target_properties(hoop PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
