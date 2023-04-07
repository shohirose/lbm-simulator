cmake_minimum_required(VERSION 3.14)
project(lbm-simulator)

# ---------------------------------- Options -----------------------------------
option(LBM_DOWNLOAD_Eigen3
    "Automatically download and build Eigen3 library."
    ON
)
option(LBM_DOWNLOAD_fmt
    "Automatically download and build fmt library."
    ON
)
option(LBM_DOWNLOAD_nlohmann_json
    "Automatically download and build nlohmann_json library."
    ON
)
option(LBM_DOWNLOAD_CLI11
    "Automatically download and build CLI11 library."
    ON
)
option(LBM_DOWNLOAD_benchmark
    "Automatically download and build Google benchmark library."
    ON
)
option(LBM_BUILD_benchmark
    "Build benchmark tests."
    OFF
)

# -------------------------------- Dependencies --------------------------------
include(FetchContent)

function(download_or_find_package package_name download_url)
    if (LBM_DOWNLOAD_${package_name})
        FetchContent_Declare(
            ${package_name}
            URL ${download_url}
            DOWNLOAD_EXTRACT_TIMESTAMP OFF
        )
        set(download_libs ${download_libs} ${package_name} PARENT_SCOPE)
    else()
        find_package(${package_name} REQUIRED)
    endif()
endfunction()

download_or_find_package(Eigen3
    https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
)

download_or_find_package(fmt
    https://github.com/fmtlib/fmt/releases/download/9.1.0/fmt-9.1.0.zip
)

download_or_find_package(nlohmann_json
    https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
)

download_or_find_package(CLI11
    https://github.com/CLIUtils/CLI11/archive/refs/tags/v2.3.2.zip
)

if (LBM_BUILD_benchmark)
    download_or_find_package(benchmark
        https://github.com/google/benchmark/archive/refs/tags/v1.7.1.zip
    )    
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
endif()

FetchContent_MakeAvailable(${download_libs})

# ------------------------------- lbm library ----------------------------------
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

# ---------------------------------- Simulator ---------------------------------
add_executable(lbm-simulator
    src/lbm_simulator.cpp
)
target_link_libraries(lbm-simulator
PRIVATE
    lbm
    CLI11::CLI11
    nlohmann_json::nlohmann_json
)

# ------------------------------- Benchmark Tests ------------------------------
if (LBM_BUILD_benchmark)
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
        target_compile_features(${name}
        PRIVATE
            cxx_std_17
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