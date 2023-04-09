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
        run_feq_benchmarks.cpp
    )

    add_benchmark(run-propagation-benchmarks
        run_propagation_benchmarks.cpp
    )
    target_link_libraries(run-propagation-benchmarks
    PRIVATE
        lbm
    )

    add_benchmark(run-collision-benchmarks
        run_collision_benchmarks.cpp
    )
endif()