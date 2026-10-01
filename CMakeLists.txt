# Standalone build of the lattice Boltzmann river solver: the library, in double
# precision for validation and single precision to match the website build, and the
# validation tests. The WebAssembly build is a separate project in wasm/.
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build -j4
#   ctest --test-dir build -L lbm -E cylinder -j4
cmake_minimum_required(VERSION 3.14)
project(lbm_river LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
find_package(Threads REQUIRED)

include(sources.cmake)
add_library(lbm STATIC ${LBM_CORE_SOURCES})
target_include_directories(lbm PUBLIC include)
target_link_libraries(lbm PUBLIC Threads::Threads)
add_library(lbm_f32 STATIC ${LBM_CORE_SOURCES})
target_include_directories(lbm_f32 PUBLIC include)
target_compile_definitions(lbm_f32 PUBLIC LBM_SINGLE_PRECISION)
target_link_libraries(lbm_f32 PUBLIC Threads::Threads)

option(LBM_BUILD_TESTS "Build the validation tests" ON)
if(LBM_BUILD_TESTS)
    enable_testing()
    foreach(test_name IN ITEMS
        kernel taylor_green poiseuille cylinder cavity smagorinsky stochastic river threads
    )
        add_executable(test_lbm_${test_name} tests/test_lbm_${test_name}.cpp)
        target_include_directories(test_lbm_${test_name} PRIVATE tools)
        target_link_libraries(test_lbm_${test_name} PRIVATE lbm)
        add_test(NAME lbm_${test_name} COMMAND test_lbm_${test_name})
        set_tests_properties(lbm_${test_name} PROPERTIES TIMEOUT 900 LABELS "lbm")
    endforeach()
    # The river again in single precision (the website's storage), over a long run.
    add_executable(test_lbm_river_f32 tests/test_lbm_river.cpp)
    target_include_directories(test_lbm_river_f32 PRIVATE tools)
    target_link_libraries(test_lbm_river_f32 PRIVATE lbm_f32)
    add_test(NAME lbm_river_f32 COMMAND test_lbm_river_f32)
    set_tests_properties(lbm_river_f32 PROPERTIES TIMEOUT 600 LABELS "lbm")
    # The two-bank lake river (HomeRiver.hpp): conservation, dots and determinism in
    # double precision; ten simulated minutes in single precision with stirring.
    add_test(NAME lbm_lake COMMAND test_lbm_river lake)
    set_tests_properties(lbm_lake PROPERTIES TIMEOUT 900 LABELS "lbm")
    add_test(NAME lbm_lake_f32 COMMAND test_lbm_river_f32 lake 600)
    set_tests_properties(lbm_lake_f32 PROPERTIES TIMEOUT 1800 LABELS "lbm")
    add_executable(test_lbm_threads_f32 tests/test_lbm_threads.cpp)
    target_include_directories(test_lbm_threads_f32 PRIVATE tools)
    target_link_libraries(test_lbm_threads_f32 PRIVATE lbm_f32)
    add_test(NAME lbm_threads_f32 COMMAND test_lbm_threads_f32)
    set_tests_properties(lbm_threads_f32 PROPERTIES TIMEOUT 120 LABELS "lbm")
endif()
