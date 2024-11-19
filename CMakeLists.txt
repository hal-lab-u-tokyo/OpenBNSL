cmake_minimum_required(VERSION 3.22)
message(STATUS "CMAKE_VERSION:${CMAKE_VERSION}" )
project(openbn LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(OpenMP REQUIRED)  

add_subdirectory(external/pybind11)
include_directories(${CMAKE_SOURCE_DIR}/include)

file(GLOB_RECURSE SOURCES src/*.cpp)
add_library(openbn_lib SHARED ${SOURCES})
target_link_libraries(openbn_lib PUBLIC OpenMP::OpenMP_CXX)  # OpenMPをリンク

add_library(openbn MODULE bindings/pybind_module.cpp)
target_link_libraries(openbn PRIVATE pybind11::module openbn_lib)
set_target_properties(openbn PROPERTIES PREFIX "")