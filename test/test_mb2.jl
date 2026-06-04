@testset "OnOffMarkov (MB2)" begin

    P1 = [0.90 0.10;
          0.20 0.80]
    P2 = [0.30 0.70;
          0.60 0.40]
    Q = [0.10 0.90;
         0.80 0.20]

    @testset "Constructor — valid" begin
        g = OnOffMarkov(1.5, [:a, :b], [P1, P2], Q; L_min = 2.0)
        @test g.alpha == 1.5
        @test g.L_min == 2.0
        @test length(g.transition_matrices) == 2
        @test g.switching_matrix == Q
    end

    @testset "Constructor — argument errors" begin
        @test_throws ArgumentError OnOffMarkov(1.0, [:a, :b], [P1, P2], Q)
        @test_throws ArgumentError OnOffMarkov(2.0, [:a, :b], [P1, P2], Q)
        @test_throws ArgumentError OnOffMarkov(1.5, [:a, :a], [P1, P2], Q)
        @test_throws ArgumentError OnOffMarkov(1.5, [:a, :b], Matrix{Float64}[], Q)
        @test_throws ArgumentError OnOffMarkov(1.5, [:a, :b], [P1], Q)
        @test_throws ArgumentError OnOffMarkov(1.5, [:a, :b], [P1, [1.0 0.2; 0.0 0.8]], Q)
        @test_throws ArgumentError OnOffMarkov(1.5, [:a, :b], [P1, P2], [1.0 0.0])
        @test_throws ArgumentError OnOffMarkov(1.5, [:a, :b], [P1, P2], Q; L_min = 0.0)
    end

    @testset "generate — output type and length" begin
        g = OnOffMarkov(1.5, ['a', 'b'], [P1, P2], Q)
        seq = generate(g, 3_000; rng = StableRNG(220))
        @test length(seq) == 3_000
        @test eltype(seq) == Char
        @test all(c ∈ ('a', 'b') for c in seq)
    end

    @testset "target_marginal and empirical marginal" begin
        g = OnOffMarkov(1.7, [:a, :b], [P1, P2], Q; L_min = 1.0)
        p = target_marginal(g)
        @test length(p) == 2
        @test isapprox(sum(p), 1.0)
        @test all(≥(0), p)

        seq = generate(g, 40_000; rng = StableRNG(221))
        @test total_variation(empirical_marginal(seq, [:a, :b]), p) < 0.15
    end

    @testset "local Markov helpers" begin
        @test isapprox(stationary_distribution([0.5 0.5; 0.25 0.75]),
                       [1/3, 2/3]; atol = 1e-8)
        @test_throws ArgumentError validate_transition_matrix([0.5 0.6; 0.2 0.8])
    end

end
