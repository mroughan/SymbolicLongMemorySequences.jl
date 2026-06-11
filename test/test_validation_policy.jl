@testset "Validation policy and benchmark infrastructure" begin
    root = normpath(joinpath(@__DIR__, ".."))

    policy = read(joinpath(root, "VALIDATION_POLICY.md"), String)
    @test contains(policy, "Fast Path")
    @test contains(policy, "Manual Validation")
    @test contains(policy, "S5_VALIDATION_LARGE")
    @test contains(policy, "S5_BENCHMARK_LARGE")
    @test contains(policy, "Future Trigram Validation")
    @test contains(policy, "centered one-hot")

    benchmark_project = read(joinpath(root, "benchmark", "Project.toml"), String)
    @test contains(benchmark_project, "BenchmarkTools")
    @test contains(benchmark_project, "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf")
    @test contains(benchmark_project, "StableRNGs")

    benchmark_readme = read(joinpath(root, "benchmark", "README.md"), String)
    @test contains(benchmark_readme, "S5_BENCHMARK_LARGE")
    @test contains(benchmark_readme, "machine-specific evidence")

    benchmark_script = read(joinpath(root, "benchmark", "benchmarks.jl"), String)
    @test contains(benchmark_script, "BenchmarkGroup")
    @test contains(benchmark_script, "S5_BENCHMARK_LARGE")
    @test contains(benchmark_script, "PB1_SpectralFGN_fft=n")
    @test contains(benchmark_script, "MB1_LAMP_d=")
    @test contains(benchmark_script, "MB3_FSS_streams=")

    validation_project = read(joinpath(root, "validation", "Project.toml"), String)
    @test contains(validation_project, "LongMemory")

    comparison_script = read(joinpath(root, "validation", "longmemory_comparison.jl"), String)
    @test contains(comparison_script, "LongMemory.autocovariance")
    @test contains(comparison_script, "LongMemory.autocorrelation")
    @test contains(comparison_script, "LongMemory.periodogram")
end
