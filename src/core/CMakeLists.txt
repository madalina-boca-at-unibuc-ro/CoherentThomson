#Build provenance (git commit, dirty state, patch, compiler, flags) for run_log.txt: regenerated on every build by
#cmake/build_info.cmake, and once now so the source exists before the first build.
string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_UPPER)
set(BUILD_INFO_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/build_info_data.cpp")
set(BUILD_INFO_ARGS
    "-DSOURCE_DIR=${CMAKE_SOURCE_DIR}"
    "-DOUTPUT=${BUILD_INFO_SOURCE}"
    "-DCOMPILER=${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}"
    "-DCXX_FLAGS=${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BUILD_TYPE_UPPER}}"
    "-DBUILD_TYPE=${CMAKE_BUILD_TYPE}"
)
execute_process(COMMAND ${CMAKE_COMMAND} ${BUILD_INFO_ARGS} -P "${CMAKE_SOURCE_DIR}/cmake/build_info.cmake")
add_custom_target(coherent_thomson_build_info
    COMMAND ${CMAKE_COMMAND} ${BUILD_INFO_ARGS} -P "${CMAKE_SOURCE_DIR}/cmake/build_info.cmake"
    BYPRODUCTS "${BUILD_INFO_SOURCE}"
    COMMENT "Recording git state for build provenance"
    VERBATIM
)

add_library(coherent_thomson_core STATIC
    detector/detector.cpp
    detector/detector_factory.cpp
    detector/detector_plotter.cpp
    laser/laser_field.cpp
    laser/laser_factory.cpp
    laser/laser_plotter.cpp
    particle/electron.cpp
    particle/electron_factory.cpp
    particle/electron_plotter.cpp
    simulation/simulation.cpp
    simulation/form_factor.cpp
    radiation/radiation.cpp
    radiation/radiation_plotter.cpp
    logging/run_log.cpp
    build_info/build_info.cpp
    ${BUILD_INFO_SOURCE}
)

add_dependencies(coherent_thomson_core coherent_thomson_build_info)

target_include_directories(coherent_thomson_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

find_package(Threads REQUIRED)
target_link_libraries(coherent_thomson_core PUBLIC Threads::Threads)

