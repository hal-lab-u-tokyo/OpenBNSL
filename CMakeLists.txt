cmake_minimum_required(VERSION 3.22)
message(STATUS "CMake version: ${CMAKE_VERSION}")

project(openbnsl_project LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# https://pybind11.readthedocs.io/en/stable/faq.html#someclass-declared-with-greater-visibility-than-the-type-of-its-field-someclass-member-wattributes
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden") 

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# ccache integration
find_program(CCACHE_PROGRAM ccache)
if (CCACHE_PROGRAM)
    message(STATUS "✅ ccache found: ${CCACHE_PROGRAM}")
    set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
else()
    message(WARNING "⚠️ ccache not found: continuing without compiler cache.")
endif()

# find required packages: Python3, OpenMP, GMP and pybind11
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)
# find_package(Threads REQUIRED)
find_package(GTest REQUIRED)
add_subdirectory(external/pybind11)
message(STATUS "pybind11 version: ${pybind11_VERSION}")

# CUDA Toolkit (optional)
find_package(CUDAToolkit)
if (CUDAToolkit_FOUND)
    message(STATUS "✅ CUDA Toolkit found: enabling CUDA support.")
    enable_language(CUDA)
    set(USE_CUDA TRUE)
    add_definitions(-DUSE_CUDA)
else()
    message(WARNING "⚠️ CUDA Toolkit not found: CUDA features will be disabled.")
    set(USE_CUDA FALSE)
endif()

# compile openbnsl_lib as shared object, which depends on OpenMP, pybind11, and Python3
file(GLOB_RECURSE CPP_SOURCES src/*.cpp)
if (USE_CUDA)
    file(GLOB_RECURSE CUDA_SOURCES src/*.cu)
    set(SOURCES ${CPP_SOURCES} ${CUDA_SOURCES})
else()
    set(SOURCES ${CPP_SOURCES})
endif()
message(STATUS "SOURCES: ${SOURCES}")
add_library(openbnsl_lib STATIC ${SOURCES})
set_target_properties(openbnsl_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(openbnsl_lib PUBLIC include)
target_link_libraries(
    openbnsl_lib 
    PUBLIC 
        OpenMP::OpenMP_CXX 
        pybind11::headers 
        Python3::Python 
)
if (USE_CUDA)
    target_link_libraries(openbnsl_lib PUBLIC CUDA::cudart)
endif()

# compile openbnsl python bindings, which depends on openbnsl_lib, OpenMP and pybind11
pybind11_add_module(openbnsl bindings/pybind_module.cpp)
target_include_directories(openbnsl PUBLIC include)
target_link_libraries(
    openbnsl 
    PUBLIC 
        openbnsl_lib 
        OpenMP::OpenMP_CXX 
        pybind11::headers 
)
set_target_properties(openbnsl PROPERTIES PREFIX "" SUFFIX ".so" LINK_FLAGS "-shared")

# compile test executable, which depends on openbnsl_lib, OpenMP, GMP, GTest
file(GLOB_RECURSE TEST_SOURCES gtest/*.cpp)
message(STATUS "TEST_SOURCES: ${TEST_SOURCES}")
add_executable(openbnsl_test ${TEST_SOURCES})
target_link_libraries(
    openbnsl_test 
    PUBLIC 
        openbnsl_lib 
        OpenMP::OpenMP_CXX 
        GTest::GTest 
        GTest::Main
) 
