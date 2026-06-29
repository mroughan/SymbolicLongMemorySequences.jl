"""
    IntermittentMapSymbols(z, alphabet [, marginal]; burnin = 1000)

Property-based symbolic generator (PB4) using a latent intermittent map.

The latent driver follows the Pomeau-Manneville-style update

    x[t+1] = (x[t] + x[t]^z) mod 1

from a random initial state. Intermittency near zero can create long laminar
episodes and broad finite-scale dependence. The generated real-valued driver is
rank-binned into `alphabet`, so finite-sample symbol counts are as close as
possible to `marginal`.

This is a latent-dynamics generator, not an exact symbolic LRD construction. The
parameter `z` controls the strength of intermittency, but S5.jl does not claim a
closed-form finite-sample Hurst parameter for this model.

# Arguments
- `z::Real`: intermittency exponent, `z > 1`.
- `alphabet`: ordered collection of symbols.
- `marginal::AbstractVector{<:Real}`: target marginal probabilities (default:
  uniform).

# Keyword Arguments
- `burnin::Int = 1000`: number of latent-map iterations discarded before
  collecting the driver.

# Complexity
O(n log n) time from rank binning, O(n) memory.

# References
Provata, A., & Beck, C. (2012). Coupled intermittent maps modelling the
statistics of genomic sequences: a network approach. arXiv:1205.2249.

# Examples
```julia
julia> g = IntermittentMapSymbols(1.6, [:A, :B], [0.4, 0.6])
IntermittentMapSymbols{Vector{Symbol}, Vector{Float64}}(z=1.6, k=2, burnin=1000)

julia> seq = generate(g, 1024; rng = MersenneTwister(1))
julia> length(seq) == 1024 && eltype(seq) == Symbol
true
```
"""
struct IntermittentMapSymbols{A, M <: AbstractVector{<:Real}} <: LRDGenerator
    z        :: Float64
    alphabet :: A
    marginal :: M
    burnin   :: Int

    function IntermittentMapSymbols{A, M}(z::Float64, alphabet::A,
                                          marginal::M,
                                          burnin::Int) where {A, M <: AbstractVector{<:Real}}
        z > 1.0 || throw(ArgumentError("z must be > 1, got $z"))
        validate_alphabet(alphabet)
        k = length(alphabet)
        length(marginal) == k ||
            throw(ArgumentError(
                "marginal length $(length(marginal)) ≠ alphabet length $k"))
        burnin ≥ 0 || throw(ArgumentError("burnin must be ≥ 0, got $burnin"))
        new{A, M}(z, alphabet, marginal, burnin)
    end
end

function IntermittentMapSymbols(z::Real, alphabet,
                                marginal::AbstractVector{<:Real} =
                                    fill(1.0 / length(alphabet), length(alphabet));
                                burnin::Int = 1000)
    m = validate_probability_vector(marginal, "marginal")
    IntermittentMapSymbols{typeof(alphabet), typeof(m)}(Float64(z), alphabet, m,
                                                        burnin)
end

function Base.show(io::IO, g::IntermittentMapSymbols)
    print(io, "IntermittentMapSymbols{$(typeof(g.alphabet)), $(typeof(g.marginal))}",
          "(z=$(g.z), k=$(length(g.alphabet)), burnin=$(g.burnin))")
end

"""
    generate(g::IntermittentMapSymbols, n; rng) -> Vector

Generate `n` symbols from an [`IntermittentMapSymbols`](@ref) generator.
"""
function generate(g::IntermittentMapSymbols, n::Int;
                  rng::AbstractRNG = Random.default_rng())
    n ≥ 4 || throw(ArgumentError(
        "IntermittentMapSymbols requires n ≥ 4 for rank binning, got $n"))
    x = rand(rng)
    for _ in 1:g.burnin
        x = mod(x + x^g.z, 1.0)
        x == 0.0 && (x = eps(Float64))
    end

    driver = Vector{Float64}(undef, n)
    @inbounds for t in 1:n
        x = mod(x + x^g.z, 1.0)
        x == 0.0 && (x = eps(Float64))
        driver[t] = x
    end
    return quantize_to_symbols(driver, g.alphabet, g.marginal)
end

"""
    generate_with_latent(g::IntermittentMapSymbols, n; rng) -> sequence, latent

Generate `n` symbols and return the one-row intermittent-map driver matrix used
before quantization.

# Examples
```julia
julia> g = IntermittentMapSymbols(1.6, [:a, :b]; burnin = 10);
julia> seq, latent = generate_with_latent(g, 16; rng = MersenneTwister(1));
julia> length(seq), size(latent)
(16, (1, 16))
```
"""
function generate_with_latent(g::IntermittentMapSymbols, n::Int;
                              rng::AbstractRNG = Random.default_rng())
    n ≥ 4 || throw(ArgumentError(
        "IntermittentMapSymbols requires n ≥ 4 for rank binning, got $n"))
    latent = generate_latent(IntermittentMapSource(g.z; burnin = g.burnin), n, 1; rng)
    return quantize_to_symbols(vec(@view latent[1, :]), g.alphabet, g.marginal), latent
end
