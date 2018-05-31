cmake_minimum_required(VERSION 3.5)

project(OpenMP3)

set(CMAKE_CXX_STANDARD 11)

set(OPENMP3_INCLUDE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include")
add_subdirectory(src)
