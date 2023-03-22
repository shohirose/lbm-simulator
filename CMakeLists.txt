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
FetchContent_MakeAvailable(eigen fmt)

add_library(lbm INTERFACE)
target_compile_features(lbm
INTERFACE
    cxx_std_17
)
target_compile_options(lbm
INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:/W3 /utf-8>
)
target_compile_definitions(lbm
INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:NOMINMAX _USE_MATH_DEFINES>
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
PUBLIC
    lbm
)