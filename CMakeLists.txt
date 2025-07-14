cmake_minimum_required(VERSION 3.22)
message(STATUS "CMake version: ${CMAKE_VERSION}")

project(openbnsl_lib LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# https://pybind11.readthedocs.io/en/stable/faq.html#someclass-declared-with-greater-visibility-than-the-type-of-its-field-someclass-member-wattributes
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden") 

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

find_program(CCACHE_PROGRAM ccache)
if (CCACHE_PROGRAM)
    set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
endif()

###########
# Options #
###########
option(BUILD_TESTS "Build test executables" OFF)
option(VERBOSE_OUTPUT "Enable verbose output" OFF)

##########################
# Find required packages #
##########################
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)
if (BUILD_TESTS)
    find_package(GTest REQUIRED)
endif()
add_subdirectory(external/pybind11)
message(STATUS "pybind11 version: ${pybind11_VERSION}")

# CUDA Toolkit (optional)
find_package(CUDAToolkit)
if (CUDAToolkit_FOUND)
    message(STATUS "✅ CUDA Toolkit found: enabling CUDA support.")
    enable_language(CUDA)
    set(USE_CUDA TRUE)
    add_definitions(-DUSE_CUDA)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wall")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3")
else()
    message(WARNING "⚠️ CUDA Toolkit not found: CUDA features will be disabled.")
    set(USE_CUDA FALSE)
endif()

##########################
# Compile core libraries #
##########################
file(GLOB_RECURSE CPP_SOURCES src/*.cpp)
list(APPEND SOURCES ${CPP_SOURCES})
if (USE_CUDA)
    file(GLOB_RECURSE CUDA_SOURCES src/*.cu)
    list(APPEND SOURCES ${CUDA_SOURCES})
endif()
if (VERBOSE_OUTPUT)
    message(STATUS "Source files included in the library:")
    foreach(SRC_FILE IN LISTS SOURCES)
        message(STATUS "  - ${SRC_FILE}")
    endforeach()
endif()
add_library(openbnsl_core STATIC ${SOURCES})
set_target_properties(openbnsl_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(openbnsl_core PUBLIC include)
target_link_libraries(
    openbnsl_core 
    PUBLIC 
        OpenMP::OpenMP_CXX 
        pybind11::headers 
        Python3::Python 
)
if (USE_CUDA)
    target_link_libraries(openbnsl_core PUBLIC CUDA::cudart)
endif()
message(STATUS "✅ Backend library configuration completed")

###########################
# Compile python bindings #
###########################
set(PYTHON_PACKAGE_NAME "openbnsllib")
pybind11_add_module(openbnsl_lib bindings/pybind_module.cpp)
set_target_properties(openbnsl_lib PROPERTIES OUTPUT_NAME ${PYTHON_PACKAGE_NAME})
target_include_directories(openbnsl_lib PUBLIC include)
target_link_libraries(
    openbnsl_lib 
    PUBLIC 
        openbnsl_core 
        OpenMP::OpenMP_CXX 
        pybind11::headers 
)
message(STATUS "✅ Python bindings configuration completed")

###########################
# Compile test executable #
###########################
if (BUILD_TESTS)
    file(GLOB_RECURSE TEST_SOURCES gtest/*.cpp)
    if (VERBOSE_OUTPUT)
        message(STATUS "Test files included in the library:")
        foreach(SRC_FILE IN LISTS TEST_SOURCES)
            message(STATUS "  - ${SRC_FILE}")
        endforeach()
    endif()
    add_executable(openbnsl_test ${TEST_SOURCES})
    target_link_libraries(
        openbnsl_test 
        PUBLIC 
            openbnsl_core 
            OpenMP::OpenMP_CXX 
            GTest::GTest 
            GTest::Main
    )
    message(STATUS "✅ Test executable configuration completed")
endif()