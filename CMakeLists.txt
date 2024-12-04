cmake_minimum_required (VERSION 3.5)
project(svd3x3)

add_definitions(-std=c++14)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(svd_demo svd_demo.cpp)
