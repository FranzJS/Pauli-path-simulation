cmake_minimum_required(VERSION 3.20)
project(pauli_path_benchmark LANGUAGES CXX)

option(PAULI_NATIVE_OPT "Enable -march=native for local benchmark builds" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(pauli_bench
    src/types.cpp
    src/kernel.cpp
    src/circuits.cpp
    src/bfs.cpp
    src/truncation.cpp
    src/hybrid.cpp
    src/reference.cpp
    src/reference_registry.cpp)

target_include_directories(pauli_bench PUBLIC include)
target_compile_options(pauli_bench PRIVATE -Wall -Wextra -Wpedantic)
target_compile_definitions(pauli_bench PRIVATE
    PAULI_REFERENCE_REGISTRY_DEFAULT_PATH="${CMAKE_CURRENT_SOURCE_DIR}/references/reference_registry.json")

if(PAULI_NATIVE_OPT)
    target_compile_options(pauli_bench PRIVATE -march=native)
endif()

add_executable(pauli_benchmark benchmarks/main.cpp)
target_link_libraries(pauli_benchmark PRIVATE pauli_bench)

add_executable(pauli_magnitude_pps_benchmark benchmarks/magnitude_pps_main.cpp)
target_link_libraries(pauli_magnitude_pps_benchmark PRIVATE pauli_bench)

add_executable(pauli_frontier_spectrum benchmarks/frontier_spectrum_main.cpp)
target_link_libraries(pauli_frontier_spectrum PRIVATE pauli_bench)

add_executable(pauli_identity_reference benchmarks/identity_reference_main.cpp)
target_link_libraries(pauli_identity_reference PRIVATE pauli_bench)

include(CTest)
if(BUILD_TESTING)
    add_executable(pauli_kernel_tests tests/test_kernel.cpp)
    target_link_libraries(pauli_kernel_tests PRIVATE pauli_bench)
    add_test(NAME pauli_kernel_tests COMMAND pauli_kernel_tests)

    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(
            NAME pauli_reference_registry_validation
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_reference_registry.py
                    ${CMAKE_CURRENT_SOURCE_DIR}/references/reference_registry.json)
    endif()
endif()
