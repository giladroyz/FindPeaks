# FindPeaks CMake Configuration
#
# This file configures the build for the find_peaks project, which provides peak finding
# implementations in both C and C++ with Python bindings for both versions.
#
# Project structure:
# - Core libraries in C and C++
# - Example executables for both implementations
# - Python bindings using pybind11 for both implementations

# Basic Project Configuration
# Specifies the minimum CMake version required to build the project
cmake_minimum_required(VERSION 3.14)
# Sets the project name to "find_peaks"
project(find_peaks)

# Sets C++11 as the standard for all C++ targets in the project
set(CMAKE_CXX_STANDARD 11)

# Config a variable with the source directory
set(FIND_PEAKS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

### Options ###
option(BUILD_PYTHON_BINDINGS "build python bindings (help with comparison to scipy)" ON)
option(BUILD_TESTS "build tests" ON)

# Core C++ Library
# Creates a static library named find_peaks_cpp from the source file
add_library(find_peaks_cpp STATIC ${FIND_PEAKS_SRC_DIR}/cpp/find_peaks.cpp)
# Adds the src directory to the include path for this library and any target that links to it
target_include_directories(find_peaks_cpp PUBLIC ${FIND_PEAKS_SRC_DIR})

# Core C Library
# Creates a static library named find_peaks_c from the source file
add_library(find_peaks_c STATIC ${FIND_PEAKS_SRC_DIR}/c/find_peaks.c)
# Adds the src directory to the include path for this library and any target that links to it
target_include_directories(find_peaks_c PUBLIC ${FIND_PEAKS_SRC_DIR})

# C++ Example Executable
# Creates an executable named example_cpp from the source file
add_executable(example_cpp ${CMAKE_CURRENT_SOURCE_DIR}/examples/usage_cpp.cpp)
# Adds the src directory to the include path for this executable
target_include_directories(example_cpp PUBLIC ${FIND_PEAKS_SRC_DIR}/cpp)
# Links the executable with the find_peaks_cpp library
target_link_libraries(example_cpp PUBLIC find_peaks_cpp)

# C Example Executable
# Creates an executable named example_c from the source file
add_executable(example_c ${CMAKE_CURRENT_SOURCE_DIR}/examples/usage_c.c)
# Adds the src directory to the include path for this executable
target_include_directories(example_c PUBLIC ${FIND_PEAKS_SRC_DIR}/c)
# Links the executable with the find_peaks_c library
target_link_libraries(example_c PUBLIC find_peaks_c)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(bindings)
endif()
