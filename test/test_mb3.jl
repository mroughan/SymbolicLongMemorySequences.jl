@testset "FSS (MB3)" begin

    @testset "Constructor — valid" begin
        g = FSS(1.5, [:a, :b, :c])
        @test g.alpha == 1.5
        @test isapprox((3 - g.alpha) / 2, 0.75)
        @test g.x_min == 1.0
        @test length(g.rates) == 3

        g2 = FSS(1.2, [:x, :y]; rates = [1.0, 3.0], x_min = 0.5)
        @test g2.rates == [1.0, 3.0]
        @test g2.x_min == 0.5
    end

    @testset "Constructor — argument errors" begin
        @test_throws ArgumentError FSS(1.0, [:a, :b])
        @test_throws ArgumentError FSS(2.0, [:a, :b])
        @test_throws ArgumentError FSS(1.5, [:a, :b]; rates = [-1.0, 1.0])
        @test_throws ArgumentError FSS(1.5, [:a, :b, :c]; rates = [1.0, 1.0])
        @test_throws ArgumentError FSS(1.5, [:a]; x_min = 0.0)
    end

    @testset "generate — output type and length" begin
        g   = FSS(1.5, ['x', 'y', 'z'])
        seq = generate(g, 2_000; rng = MersenneTwister(20))
        @test length(seq) == 2_000
        @test eltype(seq) == Char
        @test all(c ∈ ('x', 'y', 'z') for c in seq)
    end

    @testset "generate — uniform marginal (statistical)" begin
        g   = FSS(1.5, [:a, :b, :c])
        seq = generate(g, 15_000; rng = MersenneTwister(21))
        for s in (:a, :b, :c)
            @test isapprox(count(==(s), seq) / 15_000, 1/3; atol = 0.04)
        end
    end

    @testset "generate — rate-controlled marginal (statistical)" begin
        g   = FSS(1.5, [:a, :b]; rates = [1.0, 3.0])
        seq = generate(g, 20_000; rng = MersenneTwister(22))
        @test isapprox(count(==(:a), seq) / 20_000, 0.25; atol = 0.03)
        @test isapprox(count(==(:b), seq) / 20_000, 0.75; atol = 0.03)
    end

    @testset "generate — rejects n < 1" begin
        g = FSS(1.5, [:a])
        @test_throws ArgumentError generate(g, 0)
    end

    @testset "Pareto sampler — tail behaviour (statistical)" begin
        alpha = 1.5
        x_min = 1.0
        N     = 50_000
        r     = MersenneTwister(300)
        draws = [S5._pareto_sample(r, alpha, x_min) for _ in 1:N]
        @test all(>(x_min), draws)
        p90 = x_min / 0.1^(1 / alpha)         # theoretical 90th percentile
        @test isapprox(count(>(p90), draws) / N, 0.1; atol = 0.02)
    end

end
