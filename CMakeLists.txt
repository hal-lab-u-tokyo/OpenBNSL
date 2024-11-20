cmake_minimum_required(VERSION 3.22)
message(STATUS "CMake version: ${CMAKE_VERSION}")

project(openbn_project LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# find required packages: Python3, OpenMP and pybind11
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(OpenMP REQUIRED)
# find_package(Threads REQUIRED)
find_package(GTest REQUIRED)
add_subdirectory(external/pybind11)
message(STATUS "pybind11 version: ${pybind11_VERSION}")

# compile openbn_lib as shared object, which depends on OpenMP, pybind11, and Python3
file(GLOB_RECURSE SOURCES src/*.cpp)
message(STATUS "SOURCES: ${SOURCES}")
add_library(openbn_lib STATIC ${SOURCES})
set_target_properties(openbn_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(openbn_lib PUBLIC include)
target_link_libraries(openbn_lib 
    PUBLIC 
        OpenMP::OpenMP_CXX 
        pybind11::headers 
        Python3::Python
)

# compile openbn python bindings, which depends on openbn_lib, OpenMP and pybind11
pybind11_add_module(openbn bindings/pybind_module.cpp)
target_include_directories(openbn PUBLIC include)
target_link_libraries(openbn 
    PUBLIC 
        openbn_lib 
        OpenMP::OpenMP_CXX 
        pybind11::headers 
        Python3::Python
)
set_target_properties(openbn PROPERTIES PREFIX "")

# compile test executable, which depends on openbn_lib, OpenMP, GTest and GTest::Main
file(GLOB TEST_SOURCES tests/*.cpp)
message(STATUS "TEST_SOURCES: ${TEST_SOURCES}")
add_executable(openbn_test ${TEST_SOURCES})
target_link_libraries(openbn_test 
    PUBLIC 
        openbn_lib 
        OpenMP::OpenMP_CXX 
        pybind11::embed 
        Python3::Python 
        GTest::GTest 
        GTest::Main
) 
