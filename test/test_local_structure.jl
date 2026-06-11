include(joinpath(@__DIR__, "..", "validation", "local_structure.jl"))

@testset "Local structure validation" begin
    @testset "scenarios are valid transition matrices" begin
        for scenario in (iid_scenario(4), persistent_scenario(4), cyclic_scenario(4))
            @test validate_transition_matrix(scenario.transition_matrix) ==
                  scenario.transition_matrix
        end
        @test_throws ArgumentError iid_scenario(0)
        @test_throws ArgumentError persistent_scenario(0)
        @test_throws ArgumentError cyclic_scenario(1)
    end

    @testset "local structure specification pathway" begin
        spec = MarkovSpec([:a, :b], [0.9 0.1; 0.2 0.8])
        @test spec isa LocalStructureSpec
        @test local_structure_order(spec) == 1
    end

    @testset "validation script returns aggregate rows" begin
        rows = run_local_structure_control(; ns = (1_000,), ks = (2,),
                                           replicates = 3, seed = 505)
        @test length(rows) == 6
        @test all(row.tv_mean ≥ 0 for row in rows)
        @test all(row.row_tv_max ≥ row.row_tv_mean for row in rows)
    end

    @testset "identical regime specifications preserve aggregate bigrams" begin
        alphabet = [:a, :b]
        scenario = persistent_scenario(2; persistence = 0.8)
        for (seed, factory) in ((606, _identical_wavelet_markov),
                                (707, _identical_onoff_markov))
            g = factory(scenario.transition_matrix, alphabet)
            observed = empirical_bigram(generate(g, 30_000; rng = StableRNG(seed)),
                                        alphabet)
            @test maximum(rowwise_total_variation(
                observed, scenario.transition_matrix)) < 0.03
        end
    end
end
