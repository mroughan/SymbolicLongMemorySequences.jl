include(joinpath(@__DIR__, "..", "validation", "marginal_control.jl"))

@testset "Marginal control validation" begin

    @testset "validation script returns aggregate rows" begin
        rows = run_marginal_control(; ns = (500,), ks = (2,),
                                    replicates = 3, seed = 101)
        @test length(rows) == 24
        @test all(row.tv_mean ≥ 0 for row in rows)
        @test all(row.maxabs_mean ≥ 0 for row in rows)
    end

    @testset "SpectralFGN rank binning gives exact integer counts" begin
        alphabet = [:a, :b, :c]
        marginal = [0.2, 0.3, 0.5]
        g = SpectralFGN(0.8, alphabet, marginal)
        seq = generate(g, 1_003; rng = StableRNG(202))
        @test [count(==(s), seq) for s in alphabet] == bin_counts(marginal, 1_003)
        @test total_variation(empirical_marginal(seq, alphabet), target_marginal(g)) <
              1 / length(seq)
    end

    @testset "LAMP innovation improves marginal control in iid limit" begin
        alphabet = [:a, :b]
        marginal = [0.2, 0.8]
        g = LAMP(0.5, alphabet, marginal; d = 20, epsilon = 1.0)
        seq = generate(g, 20_000; rng = StableRNG(303))
        @test total_variation(empirical_marginal(seq, alphabet), marginal) < 0.02
    end

    @testset "FSS rates define target marginal" begin
        alphabet = [:a, :b]
        g = FSS(1.5, alphabet; rates = [1.0, 3.0])
        seq = generate(g, 40_000; rng = StableRNG(404))
        @test total_variation(empirical_marginal(seq, alphabet), target_marginal(g)) < 0.08
    end

end
