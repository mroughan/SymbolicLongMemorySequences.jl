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
| MB1a | `LAMP`       | Linear-Additive Markov Process         |
| MB1b | `DyadicLAMP` | Dyadic-bucket LAMP approximation       |
| MB2 | `OnOffMarkov` | Heavy-tailed regime-switching Markov   |
| MB3 | `FSS`         | Fractal Symbol Sequence via FRP/FSNP   |

# Common interface

    generate(g, n; rng = Random.default_rng()) -> Vector
    save_sequence(filepath, seq, g; created)   -> filepath

"""
module S5

using FFTW: ifft
using Dates: today
using Distributions: Pareto
using LinearAlgebra: I, mul!, norm
using Random
using Statistics: mean, std
import IncCSV

export LRDGenerator, generate, save_sequence
export LocalStructureSpec, MarkovSpec, local_structure_order
export ControlCapabilities, control_capabilities
export SpectralFGN, LGCM, WaveletMarkov, LAMP, DyadicLAMP, OnOffMarkov, FSS
export lamp_repeat_transition
export target_marginal, empirical_marginal, empirical_bigram, empirical_trigram
export bin_counts, total_variation, rowwise_total_variation
export validate_transition_matrix, stationary_distribution

include("interface.jl")
include("utils.jl")
include("pb1.jl")
include("pb2.jl")
include("pb3.jl")
include("mb1.jl")
include("mb2.jl")
include("mb3.jl")
include("controls.jl")
include("io.jl")

end # module S5
