import IncCSV

@testset "INC output (save_sequence)" begin

    @testset "SpectralFGN round-trip" begin
        mktempdir() do dir
            g    = SpectralFGN(0.8, [:a, :b, :c])
            seq  = generate(g, 200; rng = MersenneTwister(1))
            path = joinpath(dir, "pb1.inc")

            ret = save_sequence(path, seq, g)
            @test ret == path
            @test isfile(path)

            inc  = IncCSV.readinc(path)
            meta = IncCSV.metadata(inc)

            @test meta["generator"] == "SpectralFGN"
            @test meta["method"]    == "PB1"
            @test meta["package"]   == "S5"
            @test meta["n"]         == 200

            gp = meta["generator_params"]
            @test gp["H"]             == "0.8"
            @test gp["alphabet_size"] == 3

            rows = collect(IncCSV.table(inc))
            @test length(rows) == 200
            @test rows[1].index == 1
            @test rows[end].index == 200
        end
    end

    @testset "LAMP round-trip" begin
        mktempdir() do dir
            g    = LAMP(0.5, [:x, :y]; d = 50)
            seq  = generate(g, 80; rng = MersenneTwister(2))
            path = joinpath(dir, "mb1.inc")
            save_sequence(path, seq, g)

            inc  = IncCSV.readinc(path)
            meta = IncCSV.metadata(inc)
            @test meta["generator"] == "LAMP"
            @test meta["method"]    == "MB1"
            @test meta["n"]         == 80

            gp = meta["generator_params"]
            @test gp["beta"] == "0.5"
            @test gp["d"]    == 50
        end
    end

    @testset "LGCM round-trip" begin
        mktempdir() do dir
            g    = LGCM(0.8, [:a, :b], [0.25, 0.75]; calibration_iters = 4)
            seq  = generate(g, 80; rng = StableRNG(23))
            path = joinpath(dir, "pb2.inc")
            save_sequence(path, seq, g)

            inc  = IncCSV.readinc(path)
            meta = IncCSV.metadata(inc)
            @test meta["generator"] == "LGCM"
            @test meta["method"]    == "PB2"
            @test meta["n"]         == 80

            gp = meta["generator_params"]
            @test gp["H"] == "0.8"
            @test gp["calibration_iters"] == 4
        end
    end

    @testset "OnOffMarkov round-trip" begin
        mktempdir() do dir
            P1 = [0.9 0.1; 0.2 0.8]
            P2 = [0.3 0.7; 0.6 0.4]
            Q = [0.2 0.8; 0.8 0.2]
            g    = OnOffMarkov(1.5, [:a, :b], [P1, P2], Q)
            seq  = generate(g, 80; rng = StableRNG(24))
            path = joinpath(dir, "mb2.inc")
            save_sequence(path, seq, g)

            inc  = IncCSV.readinc(path)
            meta = IncCSV.metadata(inc)
            @test meta["generator"] == "OnOffMarkov"
            @test meta["method"]    == "MB2"
            @test meta["n"]         == 80

            gp = meta["generator_params"]
            @test gp["alpha"] == "1.5"
            @test gp["n_regimes"] == 2
        end
    end

    @testset "FSS round-trip" begin
        mktempdir() do dir
            g    = FSS(1.5, [:p, :q, :r])
            seq  = generate(g, 60; rng = MersenneTwister(3))
            path = joinpath(dir, "mb3.inc")
            save_sequence(path, seq, g)

            inc  = IncCSV.readinc(path)
            meta = IncCSV.metadata(inc)
            @test meta["generator"] == "FSS"
            @test meta["method"]    == "MB3"
            @test meta["n"]         == 60

            gp = meta["generator_params"]
            @test gp["alpha"]     == "1.5"
            @test gp["H_nominal"] == "0.75"
        end
    end

    @testset "sequence content preserved" begin
        mktempdir() do dir
            g    = FSS(1.4, [:a, :b])
            seq  = generate(g, 50; rng = MersenneTwister(400))
            path = joinpath(dir, "content.inc")
            save_sequence(path, seq, g)

            rows = collect(IncCSV.table(IncCSV.readinc(path)))
            syms = [r.symbol for r in rows]
            @test all(syms[i] == string(seq[i]) for i in 1:50)
        end
    end

end
