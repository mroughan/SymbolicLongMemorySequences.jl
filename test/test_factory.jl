@testset "Uniform method factory" begin

    @testset "method metadata" begin
        ids = method_ids()
        @test ids == (:PB1, :PB2, :PB3, :PB4,
                      :MB1a, :MB1b, :MB1c, :MB2, :MB3, :MB4, :MB5)
        @test method_ids(; family = :property_based) == (:PB1, :PB2, :PB3, :PB4)
        @test method_ids(; family = :model_based) == (:MB1a, :MB1b, :MB1c,
                                                       :MB2, :MB3, :MB4, :MB5)
        @test_throws ArgumentError method_ids(; family = :other)

        info = method_info(:MB5)
        @test info.id == :MB5
        @test info.type_name == :DuplicationMutation
        @test info.family == :model_based
        @test :alpha in keys(info.defaults)
        @test info.parameters === method_parameters(:MB5)
        @test method_parameters(:DuplicationMutation) === info.parameters
        @test any(p -> p.name === :mutation_probability &&
                       occursin("probability", lowercase(p.description)),
                  method_parameters(:MB5))
        @test first(method_parameters(:PB1)).name == :H
        @test first(method_parameters(:PB1)).domain == "0.5 < H < 1"
        for id in ids
            info = method_info(id)
            parameter_names = Set(p.name for p in info.parameters)
            @test issubset(Set(keys(info.defaults)), parameter_names)
            @test all(p.kind === :keyword for p in info.parameters)
            @test all(!isempty(p.domain) for p in info.parameters)
            @test all(!isempty(p.description) for p in info.parameters)
        end
        @test method_info("PB1").type_name == :SpectralFGN
        @test method_info(:SpectralFGN).id == :PB1
        @test length(method_info()) == length(ids)
        @test_throws ArgumentError method_info(:NoSuchMethod)
        @test_throws ArgumentError method_parameters(:NoSuchMethod)
    end

    @testset "construct standard cases" begin
        alphabet = [:a, :b, :c]
        marginal = [0.2, 0.3, 0.5]
        for id in method_ids()
            g = make_generator(id, alphabet; marginal)
            seq = generate(g, 128; rng = StableRNG(100 + findfirst(==(id), method_ids())))
            @test length(seq) == 128
            @test eltype(seq) == Symbol
            @test all(s in alphabet for s in seq)
        end
    end

    @testset "overrides and aliases" begin
        g1 = make_generator(:PB1, [:a, :b]; H = 0.75)
        @test g1 isa SpectralFGN
        @test g1.H == 0.75

        g2 = make_generator("MB1c", [:a, :b]; beta = 0.5, d = 20)
        @test g2 isa CalibratedAdditiveMarkov
        @test g2.beta == 0.5
        @test g2.d == 20

        g3 = make_generator(:LAMP, [:a, :b]; case = :repeat,
                            repeat_probability = 0.8)
        @test g3 isa LAMP
        @test g3.transition_matrix[1, 1] > g3.transition_matrix[1, 2]
    end

    @testset "helpful errors" begin
        @test_throws ArgumentError make_generator(:PB1, [:a, :b]; HH = 0.8)
        @test_throws ArgumentError make_generator(:PB1, [:a, :b]; case = :persistent)
        @test_throws ArgumentError make_generator(:NoSuchMethod, [:a, :b])
    end

end
