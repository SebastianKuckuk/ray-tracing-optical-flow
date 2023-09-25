###!!! THIS CMAKE LISTS FILE EXISTS MAINLY FOR IDE INTEGRATION
###!!!   PLEASE COMPILE USING THE PROVIDED MAKEFILES

cmake_minimum_required(VERSION 3.20)

project(rt-of LANGUAGES CXX)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -std=c++17")

add_executable(ray-tracing-base
        ${CMAKE_SOURCE_DIR}/src/ray-tracing-base.cpp ${CMAKE_SOURCE_DIR}/src/ray-tracing-util.h
        ${CMAKE_SOURCE_DIR}/src/vec3.h ${CMAKE_SOURCE_DIR}/src/intersect.h ${CMAKE_SOURCE_DIR}/src/conjugate-gradient-base.h
        ${CMAKE_SOURCE_DIR}/src/multigrid-base.h)
