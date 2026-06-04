"""
    SpectralFGN(H, alphabet [, marginal])

Property-based LRD symbol-sequence generator (PB1).

Synthesises fractional Gaussian noise (fGn) with Hurst parameter `H` via an
approximate spectral (FFT) method, then maps the real-valued output to symbols
using sample-quantile thresholding targeting a prescribed marginal distribution.

# Arguments
- `H::Real`: Hurst parameter, `H ∈ (0.5, 1.0)`. Higher values give stronger LRD.
- `alphabet`: ordered collection of symbols.
- `marginal::AbstractVector{<:Real}`: target marginal probabilities (default:
  uniform).

# Complexity
O(n log n) time, O(n) memory.

# Notes
Short-range structure (bigrams, etc.) is determined entirely by the quantization
scheme and cannot be prescribed independently. For joint control of LRD and local
structure see [`LAMP`](@ref).

The spectral method is *approximate*: it reproduces the asymptotic spectral slope
correctly but may deviate from the exact fGn autocovariance near lag 0.
Requires `n ≥ 4` (at least one interior frequency pair).

# References
Paxson, V. (1997). Fast, approximate synthesis of fractional Gaussian noise for
generating self-similar network traffic. *Computer Communications Review* 27, 5–18.

Dieker, T. (2004). *Simulation of fractional Brownian motion*. PhD thesis,
University of Twente.

# Examples
```julia
julia> g = SpectralFGN(0.8, [:a, :b, :c])
SpectralFGN{Vector{Symbol}, Vector{Float64}}(H=0.8, k=3)

julia> seq = generate(g, 4096; rng = MersenneTwister(1))
julia> length(seq) == 4096 && eltype(seq) == Symbol
true
```
"""
struct SpectralFGN{A, M <: AbstractVector{<:Real}} <: LRDGenerator
    H        :: Float64
    alphabet :: A
    marginal :: M

    function SpectralFGN{A, M}(H::Float64, alphabet::A,
                                marginal::M) where {A, M <: AbstractVector{<:Real}}
        (0.5 < H < 1.0) ||
            throw(ArgumentError("H must be in (0.5, 1.0), got $H"))
        validate_alphabet(alphabet)
        k = length(alphabet)
        length(marginal) == k ||
            throw(ArgumentError(
                "marginal length $(length(marginal)) ≠ alphabet length $k"))
        new{A, M}(H, alphabet, marginal)
    end
end

function SpectralFGN(H::Real, alphabet,
                     marginal::AbstractVector{<:Real} =
                         fill(1.0 / length(alphabet), length(alphabet)))
    m = validate_probability_vector(marginal, "marginal")
    SpectralFGN{typeof(alphabet), typeof(m)}(Float64(H), alphabet, m)
end

function Base.show(io::IO, g::SpectralFGN)
    print(io, "SpectralFGN{$(typeof(g.alphabet)), $(typeof(g.marginal))}",
          "(H=$(g.H), k=$(length(g.alphabet)))")
end

"""
    generate(g::SpectralFGN, n; rng) -> Vector

Generate `n` symbols from a [`SpectralFGN`](@ref) generator.
"""
function generate(g::SpectralFGN, n::Int; rng::AbstractRNG = Random.default_rng())
    n ≥ 4 || throw(ArgumentError(
        "SpectralFGN requires n ≥ 4 (at least one interior FFT frequency), got $n"))
    x = _fgn_spectral(n, g.H, rng)
    return quantize_to_symbols(x, g.alphabet, g.marginal)
end

"""
    _fgn_spectral(n, H, rng) -> Vector{Float64}

Generate length-`n` fractional Gaussian noise with Hurst parameter `H` using
Paxson's (1997) approximate spectral method.

Builds the target power spectrum S(f) ∝ |f|^(1−2H) on the DFT grid, fills with
scaled complex Gaussian noise with Hermitian symmetry, then inverse-FFTs.
Output is normalised to zero mean and unit standard deviation.
Requires `n ≥ 4`.
"""
function _fgn_spectral(n::Int, H::Float64, rng::AbstractRNG)
    n ≥ 4 || throw(ArgumentError("n must be ≥ 4, got $n"))

    Xhat  = zeros(ComplexF64, n)
    nhalf = n ÷ 2
    exp_  = 1.0 - 2H          # spectral exponent (negative for H > 0.5)

    for k in 1:(nhalf - 1)
        f = k / n
        σ = sqrt(f^exp_ / 2)
        c = complex(σ * randn(rng), σ * randn(rng))
        Xhat[k + 1]     = c
        Xhat[n - k + 1] = conj(c)
    end

    # Nyquist component (real-valued; only present for even n)
    if iseven(n)
        f = 0.5
        Xhat[nhalf + 1] = sqrt(f^exp_) * randn(rng) + 0im
    end

    x  = real(ifft(Xhat))
    σx = std(x)
    σx > 0 && (x ./= σx)
    return x .- mean(x)
end
