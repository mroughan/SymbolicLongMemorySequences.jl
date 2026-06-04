function _acf(x::AbstractVector{<:Real}, k::Int)
    n = length(x)
    μ = sum(x) / n
    v = sum((xi - μ)^2 for xi in x) / n
    v == 0.0 && return 0.0
    return sum((x[i] - μ) * (x[i + k] - μ) for i in 1:(n - k)) / (n * v)
end

@testset "LAMP (MB1)" begin

    @testset "Constructor — valid" begin
        g = LAMP(0.5, [:a, :b, :c])
        @test g.beta == 0.5
        @test g.d    == 1000
        @test isapprox(sum(g.marginal), 1.0)
        @test isapprox(sum(g.weights), 1.0)
        @test g.weights[1] > g.weights[2] > g.weights[end]

        g2 = LAMP(0.3, [:x, :y]; d = 200)
        @test g2.d == 200
    end

    @testset "Constructor — argument errors" begin
        @test_throws ArgumentError LAMP(0.0,  [:a, :b])
        @test_throws ArgumentError LAMP(1.0,  [:a, :b])
        @test_throws ArgumentError LAMP(0.5,  [:a, :b], [0.3, 0.5])
        @test_throws ArgumentError LAMP(0.5,  [:a, :b, :c], [0.5, 0.5])
        @test_throws ArgumentError LAMP(0.5,  [:a, :b]; d = 0)
    end

    @testset "generate — output type and length" begin
        g   = LAMP(0.5, ['a', 'b', 'c']; d = 100)
        seq = generate(g, 2_000; rng = MersenneTwister(10))
        @test length(seq) == 2_000
        @test eltype(seq) == Char
        @test all(c ∈ ('a', 'b', 'c') for c in seq)
    end

    @testset "generate — all symbols reachable" begin
        # LAMP with large d exhibits very long mixing times due to LRD —
        # the sample fraction for individual symbols can deviate wildly from the
        # target marginal in finite sequences.  We test only that all symbols are
        # reachable, i.e., at least one occurrence of each symbol appears across
        # multiple independent short runs.
        g = LAMP(0.7, [:a, :b, :c]; d = 50)
        observed = Set{Symbol}()
        for seed in 1:20
            seq = generate(g, 500; rng = MersenneTwister(seed))
            union!(observed, seq)
        end
        @test :a ∈ observed
        @test :b ∈ observed
        @test :c ∈ observed
    end

    @testset "generate — rejects n < 1" begin
        g = LAMP(0.5, [:a])
        @test_throws ArgumentError generate(g, 0)
    end

    @testset "ACF is positive and decreasing (statistical)" begin
        g   = LAMP(0.3, [:a, :b]; d = 500)
        seq = generate(g, 20_000; rng = MersenneTwister(200))
        x   = Float64.(seq .== :a)

        acf1  = _acf(x, 1)
        acf10 = _acf(x, 10)
        acf50 = _acf(x, 50)

        @test acf1  > 0
        @test acf10 > 0
        @test acf1  > acf10 > acf50
    end

end
