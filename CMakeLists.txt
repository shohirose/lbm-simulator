cmake_minimum_required(VERSION 3.14)
project(lbm-simulator)

option(LBM_DOWNLOAD_EIGEN "Automatically download and build Eigen library." ON)
option(LBM_DOWNLOAD_FMT "Automatically download and build fmt library." ON)
option(LBM_DOWNLOAD_JSON "Automatically download and build nlohmann::json library." ON)
option(LBM_DOWNLOAD_CLI11 "Automatically download and build CLI11 library." ON)
option(LBM_DOWNLOAD_BENCHMARK "Automatically download and build Google benchmark library." ON)
option(LBM_BUILD_BENCHMARKS "Build benchmark tests." OFF)

include(FetchContent)

if (LBM_DOWNLOAD_EIGEN)
    FetchContent_Declare(
        eigen
        URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
        DOWNLOAD_EXTRACT_TIMESTAMP OFF
    )
    list(APPEND download_libs eigen)
else()
    find_package(Eigen3 REQUIRED)
endif()

if (LBM_DOWNLOAD_FMT)
    FetchContent_Declare(
        fmt
        URL https://github.com/fmtlib/fmt/releases/download/9.1.0/fmt-9.1.0.zip
        DOWNLOAD_EXTRACT_TIMESTAMP OFF
    )
    list(APPEND download_libs fmt)
else()
    find_package(fmt REQUIRED)
endif()

if (LBM_DOWNLOAD_JSON)
    FetchContent_Declare(
        nlohmann_json
        URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
        DOWNLOAD_EXTRACT_TIMESTAMP OFF
    )
    list(APPEND download_libs nlohmann_json)
else()
    find_package(nlohmann_json REQUIRED)
endif()

if (LBM_DOWNLOAD_CLI11)
    FetchContent_Declare(
        cli11
        URL https://github.com/CLIUtils/CLI11/archive/refs/tags/v2.3.2.zip
        DOWNLOAD_EXTRACT_TIMESTAMP OFF
    )
    list(APPEND download_libs cli11)
else()
    find_package(CLI11 REQUIRED)
endif()

if (LBM_BUILD_BENCHMARKS)
    if (LBM_DOWNLOAD_BENCHMARK)
        FetchContent_Declare(
            benchmark
            URL https://github.com/google/benchmark/archive/refs/tags/v1.7.1.zip
            DOWNLOAD_EXTRACT_TIMESTAMP OFF
        )
        list(APPEND download_libs benchmark)
    else()
        find_package(benchmark REQUIRED)
    endif()
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
endif()

FetchContent_MakeAvailable(${download_libs})

add_library(lbm INTERFACE)
target_compile_features(lbm
INTERFACE
    cxx_std_17
)
target_compile_options(lbm
INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:/W3 /utf-8>
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wpedantic -Wextra>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wpedantic -Wextra>
)
target_compile_definitions(lbm
INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:NOMINMAX _USE_MATH_DEFINES>
    EIGEN_NO_AUTOMATIC_RESIZING
)
target_link_libraries(lbm
INTERFACE
    Eigen3::Eigen
    fmt::fmt
)
target_include_directories(lbm
INTERFACE
    include
)

# ------------ Simulator ---------------
add_executable(lbm-simulator
    src/lbm_simulator.cpp
)
target_link_libraries(lbm-simulator
PRIVATE
    lbm
    CLI11::CLI11
    nlohmann_json::nlohmann_json
)

# ------------ Benchmark Tests -------------
if (LBM_BUILD_BENCHMARKS)
    function(add_benchmark name)
        add_executable(${name} ${ARGN})
        target_link_libraries(${name}
        PRIVATE
            Eigen3::Eigen
            benchmark::benchmark
            benchmark::benchmark_main
        )
        target_compile_options(${name}
        PRIVATE
            $<$<CXX_COMPILER_ID:MSVC>:/W3 /utf-8>
            $<$<CXX_COMPILER_ID:GNU>:-Wall -Wpedantic -Wextra>
            $<$<CXX_COMPILER_ID:Clang>:-Wall -Wpedantic -Wextra>
        )
    endfunction()
    
    add_benchmark(run-feq-benchmarks
        benchmark/run_feq_benchmarks.cpp
    )

    add_benchmark(run-propagation-benchmarks
        benchmark/run_propagation_benchmarks.cpp
    )
    target_include_directories(run-propagation-benchmarks
    PRIVATE
        include
    )

    add_benchmark(run-collision-benchmarks
        benchmark/run_collision_benchmarks.cpp
    )
endif()