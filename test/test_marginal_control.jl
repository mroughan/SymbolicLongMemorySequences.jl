include(joinpath(@__DIR__, "..", "validation", "marginal_control.jl"))

@testset "Marginal control validation" begin

    @testset "validation script returns aggregate rows" begin
        rows = run_marginal_control(; ns = (500,), ks = (2,),
                                    replicates = 3, seed = 101)
        @test length(rows) == 33
        @test all(row.tv_mean ≥ 0 for row in rows)
        @test all(row.maxabs_mean ≥ 0 for row in rows)
    end

    @testset "uniform marginal chi-squared helpers" begin
        @test trimmed_window(100; trim_fraction = 0.1) == 11:90
        @test_throws ArgumentError trimmed_window(10; trim_fraction = 0.5)

        counts = marginal_counts([:a, :b, :a, :c, :b, :a], [:a, :b, :c])
        @test counts == [3, 2, 1]

        exact = chisq_uniform_test([10, 10, 10, 10])
        @test exact.statistic == 0.0
        @test exact.df == 3
        @test exact.pvalue == 1.0

        clustered = vcat(fill(:a, 50), fill(:b, 50), fill(:a, 50), fill(:b, 50))
        ess = categorical_effective_sample_size(clustered, [:a, :b]; maxlag = 20)
        @test 0 < ess.effective_n < length(clustered)
        @test ess.integrated_autocorrelation_time > 1

        raw = chisq_uniform_test([60, 40])
        corrected = ess_corrected_chisq_uniform_test([60, 40], 50)
        @test corrected.statistic ≈ raw.statistic / 2
        @test corrected.pvalue > raw.pvalue

        rows, histogram_rows = run_uniform_marginal_validation(;
            n = 200,
            k = 2,
            replicates = 2,
            seed = 909,
            trim_fraction = 0.1,
        )
        @test length(rows) == 11
        @test length(histogram_rows) == 22
        @test all(row.trimmed_n == 160 for row in rows)
        @test all(0 ≤ row.pvalue_min ≤ 1 for row in rows)
        @test all(0 < row.effective_n_min ≤ row.trimmed_n for row in rows)
        @test all(0 ≤ row.pvalue_ess_median ≤ 1 for row in rows)
        @test all(row.mean_full_total_variation ≥ 0 for row in rows)
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

    @testset "uniform validation separates full and trimmed rank-bin control" begin
        alphabet = Symbol.("s" .* string.(1:8))
        result = uniform_marginal_replicates(
            (p, alphabet) -> SpectralFGN(0.8, alphabet, p),
            1_000,
            alphabet;
            replicates = 2,
            seed = 505,
            trim_fraction = 0.1,
        )
        target = fill(1 / 8, 8)
        @test all(total_variation(result.full_freqs[r, :], target) == 0 for r in 1:2)
        @test any(total_variation(result.freqs[r, :], target) > 0 for r in 1:2)
    end

    @testset "LAMP innovation improves marginal control in iid limit" begin
        alphabet = [:a, :b]
        marginal = [0.2, 0.8]
        g = LAMP(0.5, alphabet, marginal; d = 20, epsilon = 1.0)
        seq = generate(g, 20_000; rng = StableRNG(303))
        @test total_variation(empirical_marginal(seq, alphabet), marginal) < 0.02
    end

    @testset "LAMP marginal validation uses marginal-stationary repeat kernel" begin
        alphabet = Symbol.("s" .* string.(1:8))
        target = fill(1 / 8, 8)
        P = _marginal_lamp_transition(target)
        @test stationary_distribution(P) ≈ target
        @test P[1, 1] > target[1]

        exact = uniform_marginal_replicates(
            (p, alphabet) -> LAMP(0.5, alphabet, p; d = 50, epsilon = 0.05,
                                  transition_matrix = _marginal_lamp_transition(p)),
            10_000,
            alphabet;
            replicates = 3,
            seed = 808,
            trim_fraction = 0.1,
        )
        dyadic = uniform_marginal_replicates(
            (p, alphabet) -> DyadicLAMP(0.5, alphabet, p; d = 10_000,
                                        epsilon = 0.05,
                                        transition_matrix =
                                            _marginal_lamp_transition(p)),
            10_000,
            alphabet;
            replicates = 3,
            seed = 808,
            trim_fraction = 0.1,
        )
        exact_tv = sum(total_variation(exact.freqs[r, :], target) for r in 1:3) / 3
        dyadic_tv = sum(total_variation(dyadic.freqs[r, :], target) for r in 1:3) / 3
        @test exact_tv < 0.03
        @test dyadic_tv < 0.03
    end

    @testset "FSS rates define target marginal" begin
        alphabet = [:a, :b]
        g = FSS(1.5, alphabet; rates = [1.0, 3.0])
        seq = generate(g, 40_000; rng = StableRNG(404))
        @test total_variation(empirical_marginal(seq, alphabet), target_marginal(g)) < 0.08
    end

end
