cmake_minimum_required (VERSION 3.5)
project(svd3x3)

# set standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-O3 -march=native -Wall")

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(svd_demo svd_demo.cpp)
add_executable(svd_bench svd_bench.cpp)
