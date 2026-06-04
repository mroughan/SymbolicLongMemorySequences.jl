@testset "Utils" begin

    rng = MersenneTwister(42)   # local to this testset

    @testset "quantize_to_symbols — type and length" begin
        x    = randn(rng, 1_000)
        syms = S5.quantize_to_symbols(x, [:a, :b, :c], [1/3, 1/3, 1/3])
        @test length(syms) == 1_000
        @test eltype(syms) == Symbol
        @test all(s ∈ (:a, :b, :c) for s in syms)
    end

    @testset "quantize_to_symbols — uniform marginal" begin
        x = randn(rng, 10_000)
        for s in [:a, :b, :c]
            freq = count(==(s), S5.quantize_to_symbols(x, [:a,:b,:c], [1/3,1/3,1/3])) / 10_000
            @test isapprox(freq, 1/3; atol = 0.02)
        end
    end

    @testset "quantize_to_symbols — non-uniform marginal" begin
        x   = randn(rng, 20_000)
        mar = [0.1, 0.4, 0.5]
        s   = S5.quantize_to_symbols(x, [1, 2, 3], mar)
        for (sym, p) in zip([1,2,3], mar)
            @test isapprox(count(==(sym), s) / 20_000, p; atol = 0.02)
        end
    end

    @testset "quantize_to_symbols — single symbol" begin
        x = randn(rng, 500)
        s = S5.quantize_to_symbols(x, [:z], [1.0])
        @test all(==(:z), s)
    end

    @testset "quantize_to_symbols — argument errors" begin
        x = randn(rng, 100)
        @test_throws ArgumentError S5.quantize_to_symbols(x, [:a,:b], [0.4, 0.4])
        @test_throws ArgumentError S5.quantize_to_symbols(x, [:a,:b,:c], [0.5, 0.5])
    end

    @testset "weighted_sample — distribution" begin
        weights = [0.1, 0.6, 0.3]
        counts  = zeros(Int, 3)
        N       = 60_000
        for _ in 1:N
            counts[S5.weighted_sample(rng, weights)] += 1
        end
        @test isapprox(counts[1] / N, 0.1; atol = 0.01)
        @test isapprox(counts[2] / N, 0.6; atol = 0.01)
        @test isapprox(counts[3] / N, 0.3; atol = 0.01)
    end

    @testset "weighted_sample — near-degenerate" begin
        cnt = sum(S5.weighted_sample(rng, [1e-15, 1.0, 1e-15]) for _ in 1:1_000)
        @test isapprox(cnt / 1_000, 2.0; atol = 0.01)
    end

end
