include(joinpath(@__DIR__, "..", "validation", "lrd_symbol_diagnostics.jl"))

@testset "LRD symbol diagnostics" begin
    seq = [:a, :b, :a, :c, :a, :b]
    x = indicator_series(seq, :a; center = false)
    @test x == [1.0, 0.0, 1.0, 0.0, 1.0, 0.0]

    xc = indicator_series(seq, :a; center = true)
    @test sum(xc) ≈ 0.0
    @test mean(abs2, xc) > 0

    series = centered_indicator_series(seq, [:a, :b, :c])
    @test length(series) == 3
    @test all(s -> isapprox(sum(s), 0.0; atol = 1e-12), series)

    y = collect(1.0:6.0)
    acv = longmemory_compatible_autocovariance(y, 3)
    μ = mean(y)
    @test acv[1] ≈ sum((y .- μ) .* (y .- μ)) / length(y)
    @test acv[2] ≈ sum((y[1:5] .- μ) .* (y[2:6] .- μ)) / length(y)

    acf = longmemory_compatible_autocorrelation(y, 3)
    @test acf[1] ≈ 1.0
    @test acf ≈ acv ./ acv[1]

    f1, p1 = fft_periodogram_cycles(y .- mean(y))
    f2, p2 = longmemory_compatible_periodogram(y .- mean(y))
    @test f1 == f2
    @test p1 ≈ p2

    acf_visual, freqs, power = indicator_diagnostics(seq, [:a, :b, :c]; maxlag = 2)
    @test length(acf_visual) == 2
    @test length(freqs) == length(power) == length(seq) ÷ 2

    acf_lm, lm_freqs, lm_power = longmemory_indicator_diagnostics(seq, [:a, :b, :c];
                                                                  maxlag = 2)
    @test length(acf_lm) == 2
    @test length(lm_freqs) == length(lm_power) == length(seq) ÷ 2
end
