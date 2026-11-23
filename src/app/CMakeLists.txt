add_executable(coherent_thomson_solver main.cpp)

target_link_libraries(coherent_thomson_solver PRIVATE coherent_thomson_core)

# Release goes to the root 'bin' folder; any other build type to bin/<type>/, so a Debug build never overwrites the
# Release binary that the scripts and configs run.
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(COHERENT_THOMSON_BIN_DIR "${CMAKE_SOURCE_DIR}/bin")
else()
    set(COHERENT_THOMSON_BIN_DIR "${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}")
endif()
set_target_properties(coherent_thomson_solver PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${COHERENT_THOMSON_BIN_DIR}")
