cmake_minimum_required(VERSION 3.14)
project(lbm-simulator)

include(FetchContent)
FetchContent_Declare(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP OFF
)
FetchContent_Declare(
    fmt
    URL https://github.com/fmtlib/fmt/releases/download/9.1.0/fmt-9.1.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP OFF
)
FetchContent_Declare(
    nlohmann_json
    URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
    DOWNLOAD_EXTRACT_TIMESTAMP OFF
)
FetchContent_Declare(
    cli11
    URL https://github.com/CLIUtils/CLI11/archive/refs/tags/v2.3.2.zip
    DOWNLOAD_EXTRACT_TIMESTAMP OFF
)
FetchContent_MakeAvailable(eigen fmt nlohmann_json cli11)

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
    eigen
    fmt
)
target_include_directories(lbm
INTERFACE
    include
)

add_executable(poiseuille-flow-2d
    src/poiseuille-flow-2d.cpp
)
target_link_libraries(poiseuille-flow-2d
PRIVATE
    lbm
    CLI11
    nlohmann_json::nlohmann_json
)

add_executable(cavity-flow-2d
    src/cavity-flow-2d.cpp
)
target_link_libraries(cavity-flow-2d
PRIVATE
    lbm
    CLI11
    nlohmann_json::nlohmann_json
)