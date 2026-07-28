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
    radiation/radiation.cpp
    radiation/radiation_plotter.cpp
)

target_include_directories(coherent_thomson_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

find_package(Threads REQUIRED)
target_link_libraries(coherent_thomson_core PUBLIC Threads::Threads)

