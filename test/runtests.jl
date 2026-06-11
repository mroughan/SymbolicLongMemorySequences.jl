using Test
using Random
using StableRNGs
using S5

@testset "S5.jl" begin
    include("test_ci_configuration.jl")
    include("test_governance.jl")
    include("test_validation_policy.jl")
    include("test_utils.jl")
    include("test_lrd_symbol_diagnostics.jl")
    include("test_pb1.jl")
    include("test_pb2.jl")
    include("test_pb3.jl")
    include("test_mb1.jl")
    include("test_mb2.jl")
    include("test_mb3.jl")
    include("test_marginal_control.jl")
    include("test_local_structure.jl")
    include("test_io.jl")
end
