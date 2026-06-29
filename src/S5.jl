"""
    S5

**Self-Similar Symbols Sequence Synthesis.**

S5.jl generates Long-Range Dependent (LRD) sequences of categorical (non-numerical)
symbols for use as ground-truth test data in LRD estimation studies.

# Generators

| ID  | Type          | Mechanism                              |
|-----|---------------|----------------------------------------|
| PB1 | `SpectralFGN` | Spectral fGn synthesis + quantization  |
| PB2 | `LGCM`        | Latent Gaussian categorical model      |
| PB3 | `WaveletMarkov` | Spectral/Haar driver + Markov regimes |
| PB4 | `IntermittentMapSymbols` | Intermittent-map driver + quantization |
| MB1a | `LAMP`       | Linear-Additive Markov Process         |
| MB1b | `DyadicLAMP` | Dyadic-bucket LAMP approximation       |
| MB1c | `CalibratedAdditiveMarkov` | Centered additive memory function |
| MB2 | `OnOffMarkov` | Heavy-tailed regime-switching Markov   |
| MB3 | `FSS`         | Fractal Symbol Sequence via FRP/FSNP   |
| MB4 | `HawkesSymbol` | Power-law self-exciting symbol process |
| MB5 | `DuplicationMutation` | Copy-mutate symbolic growth       |

# Common interface

    generate(g, n; rng = Random.default_rng()) -> Vector
    generate_with_latent(g, n; rng)             -> (Vector, Matrix)
    make_generator(id, alphabet; kwargs...)     -> LRDGenerator
    method_ids()                                -> Tuple
    method_info(id)                             -> MethodInfo
    method_parameters(id)                       -> Tuple
    save_sequence(filepath, seq, g; created)   -> filepath

# Examples
```julia
julia> g = make_generator(:PB1, [:a, :b]; H = 0.75)
SpectralFGN{Vector{Symbol}, Vector{Float64}}(H=0.75, k=2)

julia> length(generate(g, 16; rng = MersenneTwister(1)))
16
```

"""
module S5

using FFTW: ifft
using Dates: today
using Distributions: Pareto
using LinearAlgebra: I, mul!, norm
using Random
using Statistics: mean, std
import IncCSV

export LRDGenerator, generate, generate_with_latent, save_sequence
export MethodInfo, ParameterInfo, method_ids, method_info, method_parameters
export make_generator
export LocalStructureSpec, MarkovSpec, local_structure_order
export ControlCapabilities, control_capabilities
export LatentSource, Symbolizer, PropertyBasedGenerator
export SpectralFGNSource, HaarLRDSource, IntermittentMapSource
export QuantileSymbolizer, ArgmaxSymbolizer, MarkovRegimeSymbolizer
export latent_width, generate_latent, symbolize
export SpectralFGN, LGCM, WaveletMarkov, IntermittentMapSymbols
export LAMP, DyadicLAMP, CalibratedAdditiveMarkov, OnOffMarkov, FSS
export HawkesSymbol, DuplicationMutation
export lamp_repeat_transition
export target_marginal, empirical_marginal, empirical_bigram, empirical_trigram
export bin_counts, total_variation, rowwise_total_variation
export validate_transition_matrix, stationary_distribution

include("interface.jl")
include("utils.jl")
include("pb1.jl")
include("pb2.jl")
include("pb3.jl")
include("pb4.jl")
include("property_based.jl")
include("mb1.jl")
include("mb2.jl")
include("mb3.jl")
include("mb4.jl")
include("mb5.jl")
include("factory.jl")
include("controls.jl")
include("io.jl")

end # module S5
