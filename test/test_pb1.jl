using FFTW: fft
using Statistics: mean, std

@testset "SpectralFGN (PB1)" begin

    @testset "Constructor — valid" begin
        g = SpectralFGN(0.7, [:a, :b])
        @test g.H == 0.7
        @test length(g.alphabet) == 2
        @test isapprox(sum(g.marginal), 1.0)

        g2 = SpectralFGN(0.8, [1, 2, 3], [0.2, 0.3, 0.5])
        @test g2.marginal ≈ [0.2, 0.3, 0.5]
    end

    @testset "Constructor — argument errors" begin
        @test_throws ArgumentError SpectralFGN(0.5,  [:a, :b])              # H = 0.5
        @test_throws ArgumentError SpectralFGN(1.0,  [:a, :b])              # H = 1.0
        @test_throws ArgumentError SpectralFGN(0.8,  [:a, :b], [0.3, 0.5]) # sum ≠ 1
        @test_throws ArgumentError SpectralFGN(0.8,  [:a, :b, :c], [0.5, 0.5]) # length mismatch
    end

    @testset "generate — output type and length" begin
        g   = SpectralFGN(0.75, [:a, :b, :c])
        seq = generate(g, 1_000; rng = MersenneTwister(1))
        @test length(seq) == 1_000
        @test eltype(seq) == Symbol
        @test all(s ∈ (:a, :b, :c) for s in seq)
    end

    @testset "generate — uniform marginal" begin
        g   = SpectralFGN(0.8, [:a, :b, :c, :d])
        seq = generate(g, 8_000; rng = MersenneTwister(2))
        for s in (:a, :b, :c, :d)
            @test isapprox(count(==(s), seq) / 8_000, 0.25; atol = 0.04)
        end
    end

    @testset "generate — non-uniform marginal" begin
        mar = [0.1, 0.4, 0.5]
        g   = SpectralFGN(0.8, [1, 2, 3], mar)
        seq = generate(g, 10_000; rng = MersenneTwister(3))
        for (s, p) in zip([1,2,3], mar)
            @test isapprox(count(==(s), seq) / 10_000, p; atol = 0.03)
        end
    end

    @testset "generate — rejects n < 4" begin
        g = SpectralFGN(0.8, [:a])
        @test_throws ArgumentError generate(g, 3)
        @test_throws ArgumentError generate(g, 0)
    end

    @testset "_fgn_spectral — length and normalisation" begin
        r = MersenneTwister(99)
        for n in [4, 5, 100, 999, 1_000, 4_096, 4_097]
            x = S5._fgn_spectral(n, 0.8, r)
            @test length(x) == n
            @test isapprox(mean(x), 0.0; atol = 1e-10)
            @test isapprox(std(x),  1.0; atol = 1e-10)
        end
    end

    @testset "_fgn_spectral — spectral slope (statistical)" begin
        # Estimate H from log-log periodogram regression; seeded for reproducibility.
        H_target = 0.75
        n        = 16_384
        x        = S5._fgn_spectral(n, H_target, MersenneTwister(101))

        P     = abs2.(fft(x)) ./ n
        k_max = n ÷ 20
        ks    = 2:k_max
        lf    = log.(ks ./ n)
        lP    = log.(real.(P[ks .+ 1]))

        m_lf  = sum(lf) / length(lf)
        m_lP  = sum(lP) / length(lP)
        slope = sum((lf .- m_lf) .* (lP .- m_lP)) / sum((lf .- m_lf).^2)
        H_est = (1 - slope) / 2

        @test isapprox(H_est, H_target; atol = 0.12)
    end

end
