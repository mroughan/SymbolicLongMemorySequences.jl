@testset "LGCM (PB2)" begin

    @testset "Constructor — valid" begin
        g = LGCM(0.75, [:a, :b, :c], [0.2, 0.3, 0.5])
        @test g.H == 0.75
        @test g.marginal == [0.2, 0.3, 0.5]
        @test g.calibration_iters == 25

        g2 = LGCM(0.8, ["x", "y"]; calibration_iters = 5, calibration_rate = 0.4)
        @test g2.calibration_iters == 5
        @test g2.calibration_rate == 0.4
    end

    @testset "Constructor — argument errors" begin
        @test_throws ArgumentError LGCM(0.5, [:a, :b])
        @test_throws ArgumentError LGCM(1.0, [:a, :b])
        @test_throws ArgumentError LGCM(0.8, [:a, :a])
        @test_throws ArgumentError LGCM(0.8, [:a, :b], [0.4, 0.4])
        @test_throws ArgumentError LGCM(0.8, [:a, :b], [-0.1, 1.1])
        @test_throws ArgumentError LGCM(0.8, [:a, :b]; calibration_iters = -1)
        @test_throws ArgumentError LGCM(0.8, [:a, :b]; calibration_rate = 0.0)
    end

    @testset "generate — output type and length" begin
        g = LGCM(0.75, ['a', 'b', 'c'])
        seq = generate(g, 2_000; rng = StableRNG(120))
        @test length(seq) == 2_000
        @test eltype(seq) == Char
        @test all(c ∈ ('a', 'b', 'c') for c in seq)
    end

    @testset "generate — calibrated marginal" begin
        alphabet = [:a, :b, :c]
        marginal = [0.2, 0.3, 0.5]
        g = LGCM(0.8, alphabet, marginal; calibration_iters = 35)
        seq = generate(g, 10_000; rng = StableRNG(121))
        @test total_variation(empirical_marginal(seq, alphabet), marginal) < 0.04
    end

    @testset "generate — rejects n < 4" begin
        g = LGCM(0.8, [:a, :b])
        @test_throws ArgumentError generate(g, 3)
    end

    @testset "target_marginal" begin
        g = LGCM(0.8, [:a, :b], [0.1, 0.9])
        @test target_marginal(g) == [0.1, 0.9]
    end

end
