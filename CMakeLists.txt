cmake_minimum_required(VERSION 3.29)
project(AIS4104-assignments)

set(CMAKE_CXX_STANDARD 20)
add_subdirectory(math)
add_subdirectory(assignment-1)
add_subdirectory(Assignment-2)

add_dependencies(Assignment-2 math)