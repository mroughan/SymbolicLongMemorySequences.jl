"""
    WaveletMarkov(H, alphabet, transition_matrices;
                  regime_weights = uniform, cascade_depth = auto,
                  driver = :spectral)

Property-based LRD symbol-sequence generator (PB3): multiscale cascade driving
a Markov state machine.

`WaveletMarkov` generates a latent long-memory driver, rank-bins the driver into
regimes, and lets each regime select one Markov transition matrix over
`alphabet`. The default `driver = :spectral` uses the same approximate spectral
fGn synthesis as [`SpectralFGN`](@ref) before rank-binning. The legacy
`driver = :haar` path keeps the original simple Haar-style Gaussian cascade for
comparison and validation studies.

# Arguments
- `H::Real`: Hurst parameter for the latent multiscale driver, `H ∈ (0.5, 1.0)`.
- `alphabet`: ordered collection of unique symbols.
- `transition_matrices`: vector of row-stochastic `k × k` matrices, one per regime.

# Keyword Arguments
- `regime_weights`: target fraction of time spent in each regime. Defaults to
  uniform over regimes.
- `cascade_depth::Int = 0`: number of dyadic cascade levels. `0` means choose
  `floor(log2(n))` at generation time. Used only with `driver = :haar`.
- `driver::Symbol = :spectral`: latent regime driver, either `:spectral` or
  `:haar`.

# Complexity
O(n log n + n k) time with O(n + R k²) memory.

# Notes
This is a pragmatic PB3 implementation: the spectral driver gives the current
default latent LRD pathway, while the Haar-like cascade is retained as a
comparison path rather than a calibrated wavelet synthesis package. The important
interface property is present: local bigram structure is controlled by explicit
Markov matrices while a latent long-memory process controls regime persistence.

Symbol-level ACF and spectrum diagnostics only see this regime persistence when
the regimes have different observable stationary symbol distributions. If every
regime has the same stationary marginal, the latent multiscale process may be
mostly hidden from one-hot symbol diagnostics.

# Examples
```julia
julia> P1 = [0.9 0.1; 0.2 0.8];
julia> P2 = [0.3 0.7; 0.6 0.4];
julia> g = WaveletMarkov(0.8, [:a, :b], [P1, P2])
WaveletMarkov{Vector{Symbol}, Vector{Float64}}(H=0.8, k=2, R=2, driver=spectral)

julia> length(generate(g, 64; rng = MersenneTwister(1)))
64
```
"""
struct WaveletMarkov{A, W <: AbstractVector{<:Real}} <: LRDGenerator
    H                   :: Float64
    alphabet            :: A
    transition_matrices :: Vector{Matrix{Float64}}
    regime_weights      :: W
    cascade_depth       :: Int
    driver              :: Symbol

    function WaveletMarkov{A, W}(H::Float64, alphabet::A,
                                 transition_matrices::Vector{Matrix{Float64}},
                                 regime_weights::W,
                                 cascade_depth::Int,
                                 driver::Symbol) where {A, W <: AbstractVector{<:Real}}
        (0.5 < H < 1.0) ||
            throw(ArgumentError("H must be in (0.5, 1.0), got $H"))
        validate_alphabet(alphabet)
        k = length(alphabet)
        isempty(transition_matrices) &&
            throw(ArgumentError("transition_matrices must be non-empty"))
        all(size(P) == (k, k) for P in transition_matrices) ||
            throw(ArgumentError("each transition matrix must have size ($k, $k)"))
        length(regime_weights) == length(transition_matrices) ||
            throw(ArgumentError("regime_weights length must match number of regimes"))
        cascade_depth ≥ 0 ||
            throw(ArgumentError("cascade_depth must be non-negative"))
        driver ∈ (:spectral, :haar) ||
            throw(ArgumentError("driver must be :spectral or :haar, got $driver"))
        new{A, W}(H, alphabet, transition_matrices, regime_weights, cascade_depth,
                  driver)
    end
end

function WaveletMarkov(H::Real, alphabet,
                       transition_matrices::AbstractVector{<:AbstractMatrix{<:Real}};
                       regime_weights::AbstractVector{<:Real} =
                           fill(1.0 / length(transition_matrices),
                                length(transition_matrices)),
                       cascade_depth::Int = 0,
                       driver::Symbol = :spectral)
    Ps = [validate_transition_matrix(P, "transition_matrices[$i]")
          for (i, P) in enumerate(transition_matrices)]
    w = validate_probability_vector(regime_weights, "regime_weights")
    WaveletMarkov{typeof(alphabet), typeof(w)}(Float64(H), alphabet, Ps, w,
                                               cascade_depth, driver)
end

"""
    WaveletMarkov(H, specs; regime_weights = uniform, cascade_depth = auto,
                  driver = :spectral)

Construct a [`WaveletMarkov`](@ref) generator from one [`MarkovSpec`](@ref) per
latent regime. All specifications must use the same ordered alphabet.

# Examples
```julia
julia> spec = MarkovSpec([:a, :b], [0.9 0.1; 0.2 0.8]);
julia> WaveletMarkov(0.8, [spec, spec])
WaveletMarkov{Vector{Symbol}, Vector{Float64}}(H=0.8, k=2, R=2, driver=spectral)
```
"""
function WaveletMarkov(H::Real, specs::AbstractVector{<:MarkovSpec};
                       regime_weights::AbstractVector{<:Real} =
                           fill(1.0 / length(specs), length(specs)),
                       cascade_depth::Int = 0,
                       driver::Symbol = :spectral)
    alphabet, Ps = unpack_markov_specs(specs)
    return WaveletMarkov(H, alphabet, Ps; regime_weights, cascade_depth, driver)
end

function Base.show(io::IO, g::WaveletMarkov)
    print(io, "WaveletMarkov{$(typeof(g.alphabet)), $(typeof(g.regime_weights))}",
          "(H=$(g.H), k=$(length(g.alphabet)), R=$(length(g.transition_matrices)), ",
          "driver=$(g.driver))")
end

"""
    generate(g::WaveletMarkov, n; rng) -> Vector

Generate `n` symbols from a [`WaveletMarkov`](@ref) generator.
"""
function generate(g::WaveletMarkov, n::Int; rng::AbstractRNG = Random.default_rng())
    n ≥ 2 || throw(ArgumentError("n must be ≥ 2, got $n"))

    driver_values = _wavelet_markov_driver(n, g.H, g.cascade_depth, g.driver, rng)
    regimes = quantize_to_symbols(driver_values, collect(1:length(g.transition_matrices)),
                                  g.regime_weights)
    symbol = weighted_sample(rng, target_marginal(g))
    result = Vector{eltype(g.alphabet)}(undef, n)

    @inbounds for t in 1:n
        P = g.transition_matrices[regimes[t]]
        symbol = weighted_sample(rng, @view P[symbol, :])
        result[t] = g.alphabet[symbol]
    end

    return result
end

"""
    generate_with_latent(g::WaveletMarkov, n; rng) -> sequence, latent

Generate `n` symbols and return the one-row latent regime-driver matrix used
before rank-binning into Markov regimes.

# Examples
```julia
julia> P = [0.9 0.1; 0.2 0.8];
julia> g = WaveletMarkov(0.75, [:a, :b], [P, P]);
julia> seq, latent = generate_with_latent(g, 16; rng = MersenneTwister(1));
julia> length(seq), size(latent)
(16, (1, 16))
```
"""
function generate_with_latent(g::WaveletMarkov, n::Int;
                              rng::AbstractRNG = Random.default_rng())
    n ≥ 2 || throw(ArgumentError("n must be ≥ 2, got $n"))

    driver_values = _wavelet_markov_driver(n, g.H, g.cascade_depth, g.driver, rng)
    regimes = quantize_to_symbols(driver_values, collect(1:length(g.transition_matrices)),
                                  g.regime_weights)
    symbol = weighted_sample(rng, target_marginal(g))
    result = Vector{eltype(g.alphabet)}(undef, n)

    @inbounds for t in 1:n
        P = g.transition_matrices[regimes[t]]
        symbol = weighted_sample(rng, @view P[symbol, :])
        result[t] = g.alphabet[symbol]
    end

    return result, reshape(driver_values, 1, :)
end

function _wavelet_markov_driver(n::Int, H::Float64, cascade_depth::Int,
                                driver::Symbol, rng::AbstractRNG)
    if driver === :spectral
        # The spectral method needs an interior FFT frequency pair. Tiny PB3
        # samples are contract tests rather than meaningful LRD simulations, so
        # retain the Haar path there to preserve generate(g, n >= 2).
        return n ≥ 4 ? _fgn_spectral(n, H, rng) :
               _haar_lrd_driver(n, H, cascade_depth, rng)
    elseif driver === :haar
        return _haar_lrd_driver(n, H, cascade_depth, rng)
    else
        throw(ArgumentError("driver must be :spectral or :haar, got $driver"))
    end
end

function _haar_lrd_driver(n::Int, H::Float64, cascade_depth::Int,
                          rng::AbstractRNG)
    depth = cascade_depth == 0 ? max(1, floor(Int, log2(n))) : cascade_depth
    x = zeros(Float64, n)

    for level in 0:depth
        block = 2^level
        σ = block^(H - 0.5)
        pos = 1
        while pos ≤ n
            val = σ * randn(rng)
            stop = min(n, pos + block - 1)
            @inbounds for i in pos:stop
                x[i] += val
            end
            pos += block
        end
    end

    x .-= mean(x)
    sx = std(x)
    sx > 0 && (x ./= sx)
    return x
end
