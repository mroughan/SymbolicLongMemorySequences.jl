@testset "WaveletMarkov (PB3)" begin

    P1 = [0.90 0.10;
          0.25 0.75]
    P2 = [0.30 0.70;
          0.60 0.40]

    @testset "Constructor — valid" begin
        g = WaveletMarkov(0.8, [:a, :b], [P1, P2]; regime_weights = [0.4, 0.6])
        @test g.H == 0.8
        @test length(g.transition_matrices) == 2
        @test g.regime_weights == [0.4, 0.6]
        @test g.driver == :spectral

        g2 = WaveletMarkov(0.7, ["x", "y"], [P1]; cascade_depth = 4,
                           driver = :haar)
        @test g2.cascade_depth == 4
        @test g2.driver == :haar

        specs = [MarkovSpec([:a, :b], P1), MarkovSpec([:a, :b], P2)]
        g3 = WaveletMarkov(0.8, specs; regime_weights = [0.4, 0.6],
                           driver = :haar)
        @test g3.alphabet == [:a, :b]
        @test g3.transition_matrices == [P1, P2]
        @test g3.driver == :haar
    end

    @testset "Constructor — argument errors" begin
        @test_throws ArgumentError WaveletMarkov(0.5, [:a, :b], [P1, P2])
        @test_throws ArgumentError WaveletMarkov(1.0, [:a, :b], [P1, P2])
        @test_throws ArgumentError WaveletMarkov(0.8, [:a, :a], [P1, P2])
        @test_throws ArgumentError WaveletMarkov(0.8, [:a, :b], Matrix{Float64}[])
        @test_throws ArgumentError WaveletMarkov(0.8, [:a, :b], [P1, [1.0 0.2; 0.0 0.8]])
        @test_throws ArgumentError WaveletMarkov(0.8, [:a, :b], [P1, P2]; regime_weights = [1.0])
        @test_throws ArgumentError WaveletMarkov(0.8, [:a, :b], [P1, P2]; regime_weights = [0.4, 0.4])
        @test_throws ArgumentError WaveletMarkov(0.8, [:a, :b], [P1, P2]; cascade_depth = -1)
        @test_throws ArgumentError WaveletMarkov(0.8, [:a, :b], [P1, P2]; driver = :unknown)
        @test_throws ArgumentError WaveletMarkov(0.8, MarkovSpec[])
        @test_throws ArgumentError WaveletMarkov(
            0.8, [MarkovSpec([:a, :b], P1), MarkovSpec(["a", "b"], P2)])
    end

    @testset "generate — output type and length" begin
        g = WaveletMarkov(0.8, ['a', 'b'], [P1, P2])
        seq = generate(g, 2_000; rng = StableRNG(320))
        @test length(seq) == 2_000
        @test eltype(seq) == Char
        @test all(c ∈ ('a', 'b') for c in seq)
    end

    @testset "generate — iid regime matrices preserve target marginal" begin
        alphabet = [:a, :b, :c]
        marginal = [0.2, 0.3, 0.5]
        P = repeat(reshape(marginal, 1, 3), 3, 1)
        g = WaveletMarkov(0.8, alphabet, [P, P]; regime_weights = [0.25, 0.75])
        seq = generate(g, 20_000; rng = StableRNG(321))
        @test total_variation(empirical_marginal(seq, alphabet), marginal) < 0.03
        @test target_marginal(g) ≈ marginal
    end

    @testset "generate — rejects n < 2" begin
        g = WaveletMarkov(0.8, [:a, :b], [P1, P2])
        @test_throws ArgumentError generate(g, 1)
    end

    @testset "generate — both latent drivers are reproducible" begin
        g_spectral = WaveletMarkov(0.8, [:a, :b], [P1, P2]; driver = :spectral)
        g_haar = WaveletMarkov(0.8, [:a, :b], [P1, P2]; driver = :haar)

        @test generate(g_spectral, 256; rng = StableRNG(322)) ==
              generate(g_spectral, 256; rng = StableRNG(322))
        @test generate(g_haar, 256; rng = StableRNG(323)) ==
              generate(g_haar, 256; rng = StableRNG(323))
        @test generate(g_spectral, 2; rng = StableRNG(324)) isa Vector{Symbol}
    end

end
