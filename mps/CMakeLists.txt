cmake_minimum_required(VERSION 3.10)

# set the project name
project(particle-method CXX)

set(CMAKE_CXX_FLAGS "-std=c++20 -Wall")
add_subdirectory(src)
add_subdirectory(benchmarks)