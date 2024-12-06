cmake_minimum_required(VERSION 3.22)
message(STATUS "CMake version: ${CMAKE_VERSION}")

project(openbnsl_project LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# https://pybind11.readthedocs.io/en/stable/faq.html#someclass-declared-with-greater-visibility-than-the-type-of-its-field-someclass-member-wattributes
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# find required packages: Python3, OpenMP, GMP and pybind11
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(OpenMP REQUIRED)
find_library(GMP_LIB gmp REQUIRED)
# find_package(Threads REQUIRED)
find_package(GTest REQUIRED)
add_subdirectory(external/pybind11)
message(STATUS "pybind11 version: ${pybind11_VERSION}")

# compile openbnsl_lib as shared object, which depends on OpenMP, pybind11, and Python3
file(GLOB_RECURSE SOURCES src/*.cpp)
message(STATUS "SOURCES: ${SOURCES}")
add_library(openbnsl_lib STATIC ${SOURCES})
set_target_properties(openbnsl_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(openbnsl_lib PUBLIC include)
target_link_libraries(openbnsl_lib 
    PUBLIC 
        OpenMP::OpenMP_CXX 
        ${GMP_LIB} 
        pybind11::headers 
        Python3::Python
)

# compile openbnsl python bindings, which depends on openbnsl_lib, OpenMP and pybind11
pybind11_add_module(openbnsl bindings/pybind_module.cpp)
target_include_directories(openbnsl PUBLIC include)
target_link_libraries(openbnsl 
    PUBLIC 
        openbnsl_lib 
        OpenMP::OpenMP_CXX 
        ${GMP_LIB} 
        pybind11::headers 
)
set_target_properties(openbnsl PROPERTIES PREFIX "")

# compile test executable, which depends on openbnsl_lib, OpenMP, GMP, GTest
file(GLOB TEST_SOURCES gtest/*.cpp)
message(STATUS "TEST_SOURCES: ${TEST_SOURCES}")
add_executable(openbnsl_test ${TEST_SOURCES})
target_link_libraries(openbnsl_test 
    PUBLIC 
        openbnsl_lib 
        OpenMP::OpenMP_CXX 
        ${GMP_LIB} 
        GTest::GTest 
        GTest::Main
) 
