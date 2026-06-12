"""
    HawkesSymbol(beta, alphabet; baseline = ones(k), excitation = I,
                 d = 1000, c = 1.0)

Model-based symbolic generator (MB4): discrete-time Hawkes-style symbols.

At each step the generator forms a non-negative intensity for every symbol,

```math
lambda_j(t) = baseline_j +
    \\sum_{l=1}^{\\min(d,t-1)} w_l excitation_{X_{t-l},j},
```

where `w_l` is proportional to `(l + c)^(-beta)`. The next symbol is sampled
with probability proportional to these intensities.

This is a finite-history, discrete-time symbolic analogue of the Hawkes-process
word-occurrence model of Ogura, Hanada, Amano, and Kondo (2022). It is useful for
creating bursty word-like symbolic sequences: with an identity excitation matrix,
recent appearances of a symbol raise its chance of appearing again.

# Arguments
- `beta::Real`: power-law memory-kernel exponent, `beta in (0, 1)`.
- `alphabet`: ordered collection of unique symbols.

# Keyword Arguments
- `baseline::AbstractVector{<:Real}`: positive baseline symbol intensities.
- `excitation::AbstractMatrix{<:Real}`: non-negative `k x k` matrix where row
  `i` contributes excitation after observing `alphabet[i]`.
- `d::Integer = 1000`: finite history depth.
- `c::Real = 1.0`: positive kernel offset.

# Complexity
O(n * k * min(d, n)) time, O(n + d) memory.

# Notes
`target_marginal(g)` reports `baseline / sum(baseline)`, but the realized
marginal also depends on `excitation`, `beta`, `d`, and finite-sample effects.
This method should therefore be treated as having an implied marginal, not exact
marginal control.

# References
Ogura, H., Hanada, Y., Amano, H., & Kondo, M. (2022). Modeling long-range dynamic
correlations of words in written texts with Hawkes processes. *Entropy*, 24(7),
858. https://doi.org/10.3390/e24070858

# Examples
```julia
julia> g = HawkesSymbol(0.6, [:a, :b]; baseline = [1.0, 1.0],
...                     excitation = [2.0 0.0; 0.0 2.0], d = 100)
HawkesSymbol{Vector{Symbol}}(beta=0.6, k=2, d=100)

julia> seq = generate(g, 1000; rng = MersenneTwister(7))
julia> length(seq) == 1000 && eltype(seq) == Symbol
true
```
"""
struct HawkesSymbol{A} <: LRDGenerator
    beta       :: Float64
    alphabet   :: A
    baseline   :: Vector{Float64}
    excitation :: Matrix{Float64}
    d          :: Int
    c          :: Float64
    weights    :: Vector{Float64}

    function HawkesSymbol{A}(beta::Float64, alphabet::A,
                             baseline::Vector{Float64},
                             excitation::Matrix{Float64},
                             d::Int, c::Float64) where {A}
        0.0 < beta < 1.0 ||
            throw(ArgumentError("beta must be in (0, 1), got $beta"))
        validate_alphabet(alphabet)
        k = length(alphabet)
        length(baseline) == k ||
            throw(ArgumentError(
                "baseline length $(length(baseline)) != alphabet length $k"))
        size(excitation) == (k, k) ||
            throw(ArgumentError(
                "excitation must have size ($k, $k), got $(size(excitation))"))
        all(isfinite, excitation) ||
            throw(ArgumentError("excitation must contain only finite values"))
        all(>=(0), excitation) ||
            throw(ArgumentError("excitation must be non-negative"))
        d >= 1 || throw(ArgumentError("d must be >= 1, got $d"))
        isfinite(c) && c > 0 ||
            throw(ArgumentError("c must be positive and finite, got $c"))

        weights = [(lag + c)^(-beta) for lag in 1:d]
        weights ./= sum(weights)
        new{A}(beta, alphabet, baseline, excitation, d, c, weights)
    end
end

function HawkesSymbol(beta::Real, alphabet;
                      baseline::AbstractVector{<:Real} =
                          fill(1.0, length(alphabet)),
                      excitation::AbstractMatrix{<:Real} =
                          Matrix{Float64}(I, length(alphabet), length(alphabet)),
                      d::Integer = 1000,
                      c::Real = 1.0)
    b = validate_positive_vector(baseline, "baseline")
    E = Matrix{Float64}(excitation)
    return HawkesSymbol{typeof(alphabet)}(
        Float64(beta), alphabet, b, E, Int(d), Float64(c))
end

function Base.show(io::IO, g::HawkesSymbol)
    print(io, "HawkesSymbol{$(typeof(g.alphabet))}",
          "(beta=$(g.beta), k=$(length(g.alphabet)), d=$(g.d))")
end

"""
    generate(g::HawkesSymbol, n; rng) -> Vector

Generate `n` symbols from a [`HawkesSymbol`](@ref) generator.

The current implementation recomputes finite-history intensities directly. This
keeps the reference implementation simple and auditable; dyadic or recursive
approximations can later reduce the cost for very large `d`.
"""
function generate(g::HawkesSymbol, n::Int; rng::AbstractRNG = Random.default_rng())
    n >= 1 || throw(ArgumentError("n must be >= 1, got $n"))
    k = length(g.alphabet)
    result = Vector{eltype(g.alphabet)}(undef, n)
    indices = Vector{Int}(undef, n)
    intensities = similar(g.baseline)

    @inbounds for t in 1:n
        copyto!(intensities, g.baseline)
        maxlag = min(g.d, t - 1)
        for lag in 1:maxlag
            src = indices[t - lag]
            w = g.weights[lag]
            for j in 1:k
                intensities[j] += w * g.excitation[src, j]
            end
        end
        idx = weighted_sample(rng, intensities)
        indices[t] = idx
        result[t] = g.alphabet[idx]
    end

    return result
end
