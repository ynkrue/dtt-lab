cmake_minimum_required(VERSION 3.27)

project(dtt-lab VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

option(DTT_ENABLE_CLANG_TIDY "Run clang-tidy if available" OFF)
option(DTT_BUILD_CUDA "Enable CUDA experiments (future work)" OFF)
option(DTT_BUILD_BENCHMARKS "Build benchmarks (Google Benchmark)" OFF)
option(DTT_BUILD_EXAMPLES "Build example executables" OFF)
option(DTT_BUILD_TESTS "Build unit tests" ON)

if(DTT_BUILD_CUDA)
  enable_language(CUDA)
  add_compile_definitions(-DDTT_ENABLE_CUDA)
endif()

add_library(dtt STATIC
  src/core/particles.cpp
  src/tree/bounding-box.cpp
  src/tree/build.cpp
  src/force/force-engine.cpp
  src/force/cpu-tree.cpp
  src/force/mpi-tree.cpp
#  src/app/simulation-runner.cpp
#  src/integrate/euler.cpp
#  src/io/vtk-writer.cpp
#  src/comm/comm-context.cpp
#  src/comm/domain-decomposition.cpp
#  src/comm/halo-exchange.cpp
)

if(DTT_ENABLE_CUDA)
  enable_language(CUDA)
  target_sources(dtt PRIVATE
    cuda/cuda-build.cu
    cuda/cuda-traversal.cu
  )
  target_include_directories(dtt PRIVATE ${PROJECT_SOURCE_DIR}/include) # already set
  # Add any CUDA-specific compile options if needed
endif()

target_include_directories(dtt
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(MSVC)
  target_compile_options(dtt PRIVATE /W4 /permissive- /w44242 /w44061 /w44062 /wd4464)
else()
  target_compile_options(dtt PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()

if(DTT_ENABLE_CLANG_TIDY)
  find_program(CLANG_TIDY_EXE NAMES clang-tidy)
  if(CLANG_TIDY_EXE)
    set_target_properties(dtt PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
  endif()
endif()

find_program(CLANG_FORMAT_EXE NAMES clang-format)
if(CLANG_FORMAT_EXE)
  file(GLOB_RECURSE FORMAT_FILES
    ${PROJECT_SOURCE_DIR}/include/*.h
    ${PROJECT_SOURCE_DIR}/include/*.hpp
    ${PROJECT_SOURCE_DIR}/src/*.cc
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/tests/*.cc
    ${PROJECT_SOURCE_DIR}/tests/*.cpp
    ${PROJECT_SOURCE_DIR}/cuda/*.cu
  )
  add_custom_target(format
    COMMAND ${CLANG_FORMAT_EXE} -i ${FORMAT_FILES}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Running clang-format on sources"
  )
endif()

include(CTest)
if(DTT_BUILD_TESTS)
  enable_testing()

  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
  )
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  add_executable(dtt-tests
    tests/build_quadtree_test.cpp
  )
  target_link_libraries(dtt-tests PRIVATE dtt GTest::gtest_main)

  if(MSVC)
    target_compile_options(dtt-tests PRIVATE /W4 /permissive-)
  else()
    target_compile_options(dtt-tests PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  include(GoogleTest)
  gtest_discover_tests(dtt-tests)
endif()

if(DTT_BUILD_BENCHMARKS)
#  include(FetchContent)
#  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
#  FetchContent_Declare(
#    googlebenchmark
#    URL https://github.com/google/benchmark/archive/refs/tags/v1.8.3.zip
#  )
#  FetchContent_MakeAvailable(googlebenchmark)
#
#  add_executable(dtt-bench
#  )
#  target_link_libraries(dtt-bench PRIVATE dtt benchmark::benchmark benchmark::benchmark_main)
#  if(MSVC)
#    target_compile_options(dtt-bench PRIVATE /W4 /permissive-)
#  else()
#    target_compile_options(dtt-bench PRIVATE -Wall -Wextra -Wpedantic -Werror)
#  endif()
endif()

if(DTT_BUILD_EXAMPLES)
#  add_executable(dtt-sim examples/sim_main.cpp)
endif()
