struct ToyAlphabetSymbol
    id::Int
end

@testset "Utils" begin

    rng = StableRNG(42)   # local to this testset

    @testset "bin_counts" begin
        @test bin_counts([0.2, 0.3, 0.5], 10) == [2, 3, 5]
        @test bin_counts([0.2, 0.3, 0.5], 11) == [2, 3, 6]
        @test sum(bin_counts([1/3, 1/3, 1/3], 1_000)) == 1_000
        @test_throws ArgumentError bin_counts([0.5, -0.1, 0.6], 100)
        @test_throws ArgumentError bin_counts([0.5, Inf], 100)
        @test_throws ArgumentError bin_counts([0.0, 0.0], 100)
    end

    @testset "quantize_to_symbols — type and length" begin
        x    = randn(rng, 1_000)
        syms = S5.quantize_to_symbols(x, [:a, :b, :c], [1/3, 1/3, 1/3])
        @test length(syms) == 1_000
        @test eltype(syms) == Symbol
        @test all(s ∈ (:a, :b, :c) for s in syms)
    end

    @testset "quantize_to_symbols — custom alphabets" begin
        x = randn(rng, 100)

        strings = S5.quantize_to_symbols(x, ["low", "high"], [0.4, 0.6])
        @test eltype(strings) == String
        @test all(s ∈ ("low", "high") for s in strings)

        custom_alphabet = [ToyAlphabetSymbol(1), ToyAlphabetSymbol(2)]
        custom = S5.quantize_to_symbols(x, custom_alphabet, [0.5, 0.5])
        @test eltype(custom) == ToyAlphabetSymbol
        @test all(s ∈ custom_alphabet for s in custom)
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
        @test [count(==(sym), s) for sym in [1, 2, 3]] == bin_counts(mar, length(x))
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
        @test_throws ArgumentError S5.quantize_to_symbols(x, [:a,:a], [0.5, 0.5])
        @test_throws ArgumentError S5.quantize_to_symbols(x, [:a,:b], [0.5, NaN])
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

    @testset "empirical marginal, bigram, and trigram" begin
        seq = [:a, :b, :a, :a, :b]
        @test empirical_marginal(seq, [:a, :b]) == [0.6, 0.4]

        bigram = empirical_bigram(seq, [:a, :b])
        @test bigram[1, :] ≈ [1/3, 2/3]
        @test bigram[2, :] ≈ [1.0, 0.0]

        trigram = empirical_trigram(seq, [:a, :b])
        @test trigram[1, 1, :] ≈ [0.0, 1.0]
        @test trigram[1, 2, :] ≈ [1.0, 0.0]
        @test trigram[2, 1, :] ≈ [1.0, 0.0]

        @test total_variation([0.6, 0.4], [0.5, 0.5]) ≈ 0.1
        @test rowwise_total_variation(bigram, [0.5 0.5; 1.0 0.0]) ≈ [1/6, 0.0]
        @test_throws ArgumentError empirical_marginal(seq, [:a, :a])
        @test_throws ArgumentError empirical_marginal([:a, :c], [:a, :b])
        @test_throws ArgumentError rowwise_total_variation(zeros(2, 2), zeros(2, 3))
    end

end
