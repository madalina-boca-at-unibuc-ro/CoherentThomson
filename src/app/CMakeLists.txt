add_executable(coherent_thomson_solver main.cpp)

target_link_libraries(coherent_thomson_solver PRIVATE coherent_thomson_core)

#Set the runtime output directory to the root 'bin' folder
set_target_properties(coherent_thomson_solver PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")