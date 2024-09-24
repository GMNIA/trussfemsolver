cmake_minimum_required(VERSION 3.10)

# Project name
project(SimpleTestProject)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 14)

# Add GoogleTest library directory
add_subdirectory(googletest)

# Add the test executable
add_executable(simulation_test test/test_simulations.cpp src/MeshReader.cpp src/Simulation.cpp)

# Link GoogleTest libraries
target_link_libraries(simulation_test PRIVATE gtest gtest_main)

# Include GoogleTest headers
target_include_directories(simulation_test PRIVATE googletest/googletest/include)

# Include directories for your project
target_include_directories(simulation_test PRIVATE src)

