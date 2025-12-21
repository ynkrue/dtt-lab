cmake_minimum_required(VERSION 3.27)

project(dtt_lab VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

option(DTT_BUILD_TESTS "Build unit tests" ON)
option(DTT_ENABLE_CUDA "Enable CUDA experiments (future work)" OFF)
option(DTT_ENABLE_CLANG_TIDY "Run clang-tidy if available" OFF)
option(DTT_BUILD_BENCHMARKS "Build benchmarks (Google Benchmark)" ON)
option(DTT_BUILD_EXAMPLES "Build example executables" ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(DTT_ENABLE_CUDA)
  enable_language(CUDA)
endif()

find_package(BLAS REQUIRED)

add_library(dtt STATIC
  src/dtt.cpp
  src/sim/particles.cpp
  src/sim/forces.cpp
  src/io/vtk.cpp
)

target_include_directories(dtt
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(TARGET BLAS::BLAS)
  target_link_libraries(dtt PUBLIC BLAS::BLAS)
else()
  target_link_libraries(dtt PUBLIC ${BLAS_LIBRARIES})
endif()

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

  add_executable(dtt_tests
    tests/dummy_test.cpp
    tests/particles_test.cpp
    tests/forces_test.cpp
  )
  target_link_libraries(dtt_tests PRIVATE dtt GTest::gtest_main)

  if(MSVC)
    target_compile_options(dtt_tests PRIVATE /W4 /permissive-)
  else()
    target_compile_options(dtt_tests PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  include(GoogleTest)
  gtest_discover_tests(dtt_tests)
endif()

if(DTT_BUILD_BENCHMARKS)
  include(FetchContent)
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(
    googlebenchmark
    URL https://github.com/google/benchmark/archive/refs/tags/v1.8.3.zip
  )
  FetchContent_MakeAvailable(googlebenchmark)

  add_executable(dtt_bench
    benchmarks/naive_bench.cpp
  )
  target_link_libraries(dtt_bench PRIVATE dtt benchmark::benchmark benchmark::benchmark_main)
  if(MSVC)
    target_compile_options(dtt_bench PRIVATE /W4 /permissive-)
  else()
    target_compile_options(dtt_bench PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()
endif()

if(DTT_BUILD_EXAMPLES)
  add_executable(dtt_sim examples/sim_main.cpp)
  target_link_libraries(dtt_sim PRIVATE dtt)
  if(MSVC)
    target_compile_options(dtt_sim PRIVATE /W4 /permissive-)
  else()
    target_compile_options(dtt_sim PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()
endif()
