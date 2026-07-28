cmake_minimum_required(VERSION 3.22)
project(
    CoherentThomson
    LANGUAGES CXX
    VERSION 0.1.0
    DESCRIPTION "Coherent Thomson Scattering"
#    HOMEPAGE "https://github.com/VindictaLab/CoherentThomson"
)


#Enforce out of source builds
if(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    message(FATAL_ERROR "Build directory must be outside the source directory.")
endif()

#Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

#Add optimization flags
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wpedantic    -march=native -mtune=native")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "-g -Wall -Wextra -Wpedantic")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_CXX_FLAGS "-O2 -g")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
    set(CMAKE_CXX_FLAGS "-Os")
endif()


#Add subdirectories
add_subdirectory(src/core)
add_subdirectory(src/app)

