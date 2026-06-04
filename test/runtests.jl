using Test
using Random
using S5

@testset "S5.jl" begin
    include("test_utils.jl")
    include("test_pb1.jl")
    include("test_mb1.jl")
    include("test_mb3.jl")
    include("test_io.jl")
end
