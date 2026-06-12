@testset "HawkesSymbol (MB4)" begin

    @testset "Constructor -- valid" begin
        g = HawkesSymbol(0.6, [:a, :b])
        @test g.beta == 0.6
        @test g.d == 1000
        @test g.c == 1.0
        @test g.baseline == [1.0, 1.0]
        @test g.excitation == [1.0 0.0; 0.0 1.0]
        @test isapprox(sum(g.weights), 1.0)

        E = [2.0 0.5; 0.25 3.0]
        g2 = HawkesSymbol(0.4, [:x, :y]; baseline = [1.0, 3.0],
                           excitation = E, d = 50, c = 2.0)
        @test g2.baseline == [1.0, 3.0]
        @test g2.excitation == E
        @test g2.d == 50
        @test g2.c == 2.0
    end

    @testset "Constructor -- argument errors" begin
        @test_throws ArgumentError HawkesSymbol(0.0, [:a, :b])
        @test_throws ArgumentError HawkesSymbol(1.0, [:a, :b])
        @test_throws ArgumentError HawkesSymbol(0.5, [:a, :a])
        @test_throws ArgumentError HawkesSymbol(0.5, [:a, :b]; baseline = [1.0])
        @test_throws ArgumentError HawkesSymbol(0.5, [:a, :b]; baseline = [0.0, 1.0])
        @test_throws ArgumentError HawkesSymbol(0.5, [:a, :b]; baseline = [1.0, Inf])
        @test_throws ArgumentError HawkesSymbol(0.5, [:a, :b]; excitation = [1.0 0.0])
        @test_throws ArgumentError HawkesSymbol(0.5, [:a, :b];
                                                excitation = [1.0 -0.1; 0.0 1.0])
        @test_throws ArgumentError HawkesSymbol(0.5, [:a, :b];
                                                excitation = [1.0 Inf; 0.0 1.0])
        @test_throws ArgumentError HawkesSymbol(0.5, [:a, :b]; d = 0)
        @test_throws ArgumentError HawkesSymbol(0.5, [:a, :b]; c = 0.0)
    end

    @testset "generate -- output type and length" begin
        g = HawkesSymbol(0.6, ['x', 'y', 'z']; d = 20)
        seq = generate(g, 2_000; rng = StableRNG(41))
        @test length(seq) == 2_000
        @test eltype(seq) == Char
        @test all(c in ('x', 'y', 'z') for c in seq)
    end

    @testset "generate -- rejects n < 1" begin
        g = HawkesSymbol(0.6, [:a])
        @test_throws ArgumentError generate(g, 0)
    end

    @testset "self excitation creates conservative burstiness signal" begin
        iid = HawkesSymbol(0.6, [:a, :b]; excitation = zeros(2, 2), d = 50)
        bursty = HawkesSymbol(0.6, [:a, :b];
                               excitation = [8.0 0.0; 0.0 8.0], d = 50)
        seq_iid = generate(iid, 5_000; rng = StableRNG(42))
        seq_bursty = generate(bursty, 5_000; rng = StableRNG(42))
        repeat_rate(seq) = count(seq[i] == seq[i - 1] for i in 2:length(seq)) /
                           (length(seq) - 1)
        @test repeat_rate(seq_bursty) > repeat_rate(seq_iid) + 0.05
    end

    @testset "target_marginal and capabilities" begin
        g = HawkesSymbol(0.6, [:a, :b]; baseline = [1.0, 3.0], d = 20)
        @test target_marginal(g) == [0.25, 0.75]
        caps = control_capabilities(g)
        @test caps.alphabet == :exact
        @test caps.marginal == :implied
        @test caps.lrd == :finite_history
    end

end
