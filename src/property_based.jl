"""
    LatentSource

Abstract supertype for numerical latent sources used by
[`PropertyBasedGenerator`](@ref).

A latent source produces one or more real-valued series that carry the
large-scale dependence structure. A [`Symbolizer`](@ref) then maps those
numerical series to symbols.

# Examples
```julia
julia> SpectralFGNSource(0.8) isa LatentSource
true
```
"""
abstract type LatentSource end

"""
    Symbolizer

Abstract supertype for transforms that map numerical latent series to symbols.

Symbolizers declare the number of latent input series they need through
[`latent_width`](@ref). For example, [`QuantileSymbolizer`](@ref) needs one
latent series, while [`ArgmaxSymbolizer`](@ref) needs one latent series per
alphabet symbol.

# Examples
```julia
julia> QuantileSymbolizer([:a, :b]) isa Symbolizer
true
```
"""
abstract type Symbolizer end

"""
    SpectralFGNSource(H)

Numerical latent source using the approximate spectral fGn construction also
used by [`SpectralFGN`](@ref).

The source can generate any positive number of independent latent streams. Each
stream uses Hurst parameter `H`. This makes it compatible with one-stream
symbolizers such as [`QuantileSymbolizer`](@ref) and multi-stream symbolizers
such as [`ArgmaxSymbolizer`](@ref).

# Examples
```julia
julia> src = SpectralFGNSource(0.75)
SpectralFGNSource(H=0.75)

julia> size(generate_latent(src, 8, 2; rng = MersenneTwister(1)))
(2, 8)
```
"""
struct SpectralFGNSource <: LatentSource
    H :: Float64
    function SpectralFGNSource(H::Float64)
        (0.5 < H < 1.0) ||
            throw(ArgumentError("H must be in (0.5, 1.0), got $H"))
        new(H)
    end
end
SpectralFGNSource(H::Real) = SpectralFGNSource(Float64(H))
Base.show(io::IO, s::SpectralFGNSource) = print(io, "SpectralFGNSource(H=$(s.H))")

"""
    HaarLRDSource(H; cascade_depth = 0)

Numerical latent source using the simple Haar-like cascade retained for PB3
comparison studies.

`cascade_depth = 0` means choose `floor(log2(n))` at generation time. This is a
pragmatic finite-sample latent driver, not a calibrated wavelet synthesis
package.

# Examples
```julia
julia> src = HaarLRDSource(0.8; cascade_depth = 3)
HaarLRDSource(H=0.8, cascade_depth=3)

julia> size(generate_latent(src, 8, 1; rng = MersenneTwister(1)))
(1, 8)
```
"""
struct HaarLRDSource <: LatentSource
    H             :: Float64
    cascade_depth :: Int
    function HaarLRDSource(H::Float64, cascade_depth::Int)
        (0.5 < H < 1.0) ||
            throw(ArgumentError("H must be in (0.5, 1.0), got $H"))
        cascade_depth ≥ 0 ||
            throw(ArgumentError("cascade_depth must be non-negative"))
        new(H, cascade_depth)
    end
end
HaarLRDSource(H::Real; cascade_depth::Int = 0) =
    HaarLRDSource(Float64(H), cascade_depth)
Base.show(io::IO, s::HaarLRDSource) =
    print(io, "HaarLRDSource(H=$(s.H), cascade_depth=$(s.cascade_depth))")

"""
    IntermittentMapSource(z; burnin = 1000)

One-stream numerical latent source using the Pomeau-Manneville-style
intermittent map also used by [`IntermittentMapSymbols`](@ref).

Because this source represents one deterministic latent map trajectory, it is
compatible with one-stream symbolizers such as [`QuantileSymbolizer`](@ref) and
[`MarkovRegimeSymbolizer`](@ref), but not with [`ArgmaxSymbolizer`](@ref).

# Examples
```julia
julia> src = IntermittentMapSource(1.6; burnin = 10)
IntermittentMapSource(z=1.6, burnin=10)

julia> size(generate_latent(src, 8, 1; rng = MersenneTwister(1)))
(1, 8)
```
"""
struct IntermittentMapSource <: LatentSource
    z      :: Float64
    burnin :: Int
    function IntermittentMapSource(z::Float64, burnin::Int)
        z > 1.0 || throw(ArgumentError("z must be > 1, got $z"))
        burnin ≥ 0 || throw(ArgumentError("burnin must be ≥ 0, got $burnin"))
        new(z, burnin)
    end
end
IntermittentMapSource(z::Real; burnin::Int = 1000) =
    IntermittentMapSource(Float64(z), burnin)
Base.show(io::IO, s::IntermittentMapSource) =
    print(io, "IntermittentMapSource(z=$(s.z), burnin=$(s.burnin))")

"""
    QuantileSymbolizer(alphabet [, marginal])

One-stream symbolizer using rank/quantile binning.

The sorted latent values are assigned to `alphabet` with finite-sample counts as
close as possible to `marginal`.

# Examples
```julia
julia> sym = QuantileSymbolizer([:a, :b], [0.25, 0.75])
QuantileSymbolizer{Vector{Symbol}, Vector{Float64}}(k=2)

julia> latent_width(sym)
1
```
"""
struct QuantileSymbolizer{A, M <: AbstractVector{<:Real}} <: Symbolizer
    alphabet :: A
    marginal :: M
    function QuantileSymbolizer{A, M}(alphabet::A,
                                      marginal::M) where {A, M <: AbstractVector{<:Real}}
        validate_alphabet(alphabet)
        length(marginal) == length(alphabet) ||
            throw(ArgumentError("marginal length must match alphabet length"))
        new{A, M}(alphabet, marginal)
    end
end

function QuantileSymbolizer(alphabet,
                            marginal::AbstractVector{<:Real} =
                                fill(1.0 / length(alphabet), length(alphabet)))
    m = validate_probability_vector(marginal, "marginal")
    return QuantileSymbolizer{typeof(alphabet), typeof(m)}(alphabet, m)
end

Base.show(io::IO, s::QuantileSymbolizer) =
    print(io, "QuantileSymbolizer{$(typeof(s.alphabet)), $(typeof(s.marginal))}",
          "(k=$(length(s.alphabet)))")

"""
    ArgmaxSymbolizer(alphabet [, marginal]; calibration_iters = 25,
                     calibration_rate = 0.7)

Multi-stream symbolizer using calibrated argmax over one latent series per
alphabet symbol.

This is the symbolization transform used by [`LGCM`](@ref). It requires
`latent_width(symbolizer) == length(alphabet)`.

# Examples
```julia
julia> sym = ArgmaxSymbolizer([:a, :b, :c])
ArgmaxSymbolizer{Vector{Symbol}, Vector{Float64}}(k=3)

julia> latent_width(sym)
3
```
"""
struct ArgmaxSymbolizer{A, M <: AbstractVector{<:Real}} <: Symbolizer
    alphabet          :: A
    marginal          :: M
    calibration_iters :: Int
    calibration_rate  :: Float64
    function ArgmaxSymbolizer{A, M}(alphabet::A, marginal::M,
                                    calibration_iters::Int,
                                    calibration_rate::Float64) where {A, M <: AbstractVector{<:Real}}
        validate_alphabet(alphabet)
        length(marginal) == length(alphabet) ||
            throw(ArgumentError("marginal length must match alphabet length"))
        calibration_iters ≥ 0 ||
            throw(ArgumentError("calibration_iters must be non-negative"))
        calibration_rate > 0 && isfinite(calibration_rate) ||
            throw(ArgumentError("calibration_rate must be positive and finite"))
        new{A, M}(alphabet, marginal, calibration_iters, calibration_rate)
    end
end

function ArgmaxSymbolizer(alphabet,
                          marginal::AbstractVector{<:Real} =
                              fill(1.0 / length(alphabet), length(alphabet));
                          calibration_iters::Int = 25,
                          calibration_rate::Real = 0.7)
    m = validate_probability_vector(marginal, "marginal")
    return ArgmaxSymbolizer{typeof(alphabet), typeof(m)}(
        alphabet, m, calibration_iters, Float64(calibration_rate))
end

Base.show(io::IO, s::ArgmaxSymbolizer) =
    print(io, "ArgmaxSymbolizer{$(typeof(s.alphabet)), $(typeof(s.marginal))}",
          "(k=$(length(s.alphabet)))")

"""
    MarkovRegimeSymbolizer(alphabet, transition_matrices;
                           regime_weights = uniform)

One-stream symbolizer that rank-bins a latent driver into regimes, then emits
symbols from regime-specific Markov transition matrices.

This is the symbolization transform used by [`WaveletMarkov`](@ref). The latent
source supplies regime persistence; the transition matrices supply local
bigram structure.

# Examples
```julia
julia> P1 = [0.9 0.1; 0.2 0.8];
julia> P2 = [0.3 0.7; 0.6 0.4];
julia> sym = MarkovRegimeSymbolizer([:a, :b], [P1, P2])
MarkovRegimeSymbolizer{Vector{Symbol}, Vector{Float64}}(k=2, R=2)

julia> latent_width(sym)
1
```
"""
struct MarkovRegimeSymbolizer{A, W <: AbstractVector{<:Real}} <: Symbolizer
    alphabet            :: A
    transition_matrices :: Vector{Matrix{Float64}}
    regime_weights      :: W
    function MarkovRegimeSymbolizer{A, W}(
        alphabet::A, transition_matrices::Vector{Matrix{Float64}},
        regime_weights::W) where {A, W <: AbstractVector{<:Real}}
        validate_alphabet(alphabet)
        k = length(alphabet)
        isempty(transition_matrices) &&
            throw(ArgumentError("transition_matrices must be non-empty"))
        all(size(P) == (k, k) for P in transition_matrices) ||
            throw(ArgumentError("each transition matrix must have size ($k, $k)"))
        length(regime_weights) == length(transition_matrices) ||
            throw(ArgumentError("regime_weights length must match number of regimes"))
        new{A, W}(alphabet, transition_matrices, regime_weights)
    end
end

function MarkovRegimeSymbolizer(
    alphabet, transition_matrices::AbstractVector{<:AbstractMatrix{<:Real}};
    regime_weights::AbstractVector{<:Real} =
        fill(1.0 / length(transition_matrices), length(transition_matrices)))
    Ps = [validate_transition_matrix(P, "transition_matrices[$i]")
          for (i, P) in enumerate(transition_matrices)]
    w = validate_probability_vector(regime_weights, "regime_weights")
    return MarkovRegimeSymbolizer{typeof(alphabet), typeof(w)}(alphabet, Ps, w)
end

Base.show(io::IO, s::MarkovRegimeSymbolizer) =
    print(io, "MarkovRegimeSymbolizer{$(typeof(s.alphabet)), ",
          "$(typeof(s.regime_weights))}(k=$(length(s.alphabet)), ",
          "R=$(length(s.transition_matrices)))")

"""
    PropertyBasedGenerator(source, symbolizer)

Composable property-based generator.

Property-based synthesis has two layers: a numerical [`LatentSource`](@ref)
that carries the large-scale dependence, and a [`Symbolizer`](@ref) that maps
the latent series to a finite alphabet. Not every source can feed every
symbolizer; construction checks the required [`latent_width`](@ref).

Named generators such as [`SpectralFGN`](@ref), [`LGCM`](@ref),
[`WaveletMarkov`](@ref), and [`IntermittentMapSymbols`](@ref) remain the stable
standard cases. `PropertyBasedGenerator` exposes the lower-level composition
path for controlled experiments.

# Examples
```julia
julia> src = SpectralFGNSource(0.8);
julia> sym = QuantileSymbolizer([:a, :b], [0.25, 0.75]);
julia> g = PropertyBasedGenerator(src, sym)
PropertyBasedGenerator(source=SpectralFGNSource, symbolizer=QuantileSymbolizer)

julia> length(generate(g, 16; rng = MersenneTwister(1)))
16
```
"""
struct PropertyBasedGenerator{S <: LatentSource, T <: Symbolizer} <: LRDGenerator
    source     :: S
    symbolizer :: T
    function PropertyBasedGenerator(source::S,
                                    symbolizer::T) where {S <: LatentSource,
                                                          T <: Symbolizer}
        width = latent_width(symbolizer)
        _supports_latent_width(source, width) ||
            throw(ArgumentError(
                "$(typeof(source)) cannot generate the $width latent stream(s) " *
                "required by $(typeof(symbolizer))"))
        new{S, T}(source, symbolizer)
    end
end

Base.show(io::IO, g::PropertyBasedGenerator) =
    print(io, "PropertyBasedGenerator(source=$(nameof(typeof(g.source))), ",
          "symbolizer=$(nameof(typeof(g.symbolizer))))")

"""
    latent_width(symbolizer) -> Int

Return the number of latent numerical series required by `symbolizer`.

# Examples
```julia
julia> latent_width(QuantileSymbolizer([:a, :b]))
1

julia> latent_width(ArgmaxSymbolizer([:a, :b]))
2
```
"""
latent_width(::QuantileSymbolizer) = 1
latent_width(::MarkovRegimeSymbolizer) = 1
latent_width(s::ArgmaxSymbolizer) = length(s.alphabet)

"""
    generate_latent(source, n, width; rng) -> Matrix{Float64}

Generate a `width × n` matrix of numerical latent series from `source`.

# Examples
```julia
julia> latent = generate_latent(SpectralFGNSource(0.75), 8, 2;
...                            rng = MersenneTwister(1));
julia> size(latent)
(2, 8)
```
"""
function generate_latent(source::SpectralFGNSource, n::Int, width::Int;
                         rng::AbstractRNG = Random.default_rng())
    n ≥ 4 || throw(ArgumentError(
        "SpectralFGNSource requires n ≥ 4 (for latent fGn synthesis), got $n"))
    width ≥ 1 || throw(ArgumentError("width must be positive, got $width"))
    latent = Matrix{Float64}(undef, width, n)
    @inbounds for i in 1:width
        latent[i, :] .= _fgn_spectral(n, source.H, rng)
    end
    return latent
end

function generate_latent(source::HaarLRDSource, n::Int, width::Int;
                         rng::AbstractRNG = Random.default_rng())
    n ≥ 2 || throw(ArgumentError("HaarLRDSource requires n ≥ 2, got $n"))
    width ≥ 1 || throw(ArgumentError("width must be positive, got $width"))
    latent = Matrix{Float64}(undef, width, n)
    @inbounds for i in 1:width
        latent[i, :] .= _haar_lrd_driver(n, source.H, source.cascade_depth, rng)
    end
    return latent
end

function generate_latent(source::IntermittentMapSource, n::Int, width::Int;
                         rng::AbstractRNG = Random.default_rng())
    n ≥ 4 || throw(ArgumentError(
        "IntermittentMapSource requires n ≥ 4 for rank binning, got $n"))
    width == 1 ||
        throw(ArgumentError("IntermittentMapSource supports width 1, got $width"))
    x = rand(rng)
    for _ in 1:source.burnin
        x = mod(x + x^source.z, 1.0)
        x == 0.0 && (x = eps(Float64))
    end
    latent = Matrix{Float64}(undef, 1, n)
    @inbounds for t in 1:n
        x = mod(x + x^source.z, 1.0)
        x == 0.0 && (x = eps(Float64))
        latent[1, t] = x
    end
    return latent
end

"""
    symbolize(symbolizer, latent; rng) -> Vector

Map a `width × n` latent matrix to symbols using `symbolizer`.

# Examples
```julia
julia> latent = reshape([0.1, 0.8, 0.2, 0.9], 1, 4);
julia> symbolize(QuantileSymbolizer([:a, :b]), latent; rng = MersenneTwister(1))
4-element Vector{Symbol}:
 :a
 :b
 :a
 :b
```
"""
function symbolize(symbolizer::QuantileSymbolizer, latent::AbstractMatrix{<:Real};
                   rng::AbstractRNG = Random.default_rng())
    size(latent, 1) == 1 ||
        throw(ArgumentError("QuantileSymbolizer requires one latent stream"))
    return quantize_to_symbols(vec(@view latent[1, :]), symbolizer.alphabet,
                               symbolizer.marginal)
end

function symbolize(symbolizer::ArgmaxSymbolizer, latent::AbstractMatrix{<:Real};
                   rng::AbstractRNG = Random.default_rng())
    k, n = size(latent)
    k == length(symbolizer.alphabet) ||
        throw(ArgumentError("ArgmaxSymbolizer requires one latent stream per symbol"))
    dense_latent = Matrix{Float64}(latent)
    offsets = log.(symbolizer.marginal .+ eps(Float64))
    _calibrate_lgcm_offsets!(offsets, dense_latent, symbolizer.marginal,
                             symbolizer.calibration_iters,
                             symbolizer.calibration_rate)
    result = Vector{eltype(symbolizer.alphabet)}(undef, n)
    @inbounds for t in 1:n
        result[t] = symbolizer.alphabet[_argmax_with_offsets(dense_latent, offsets, t)]
    end
    return result
end

function symbolize(symbolizer::MarkovRegimeSymbolizer,
                   latent::AbstractMatrix{<:Real};
                   rng::AbstractRNG = Random.default_rng())
    size(latent, 1) == 1 ||
        throw(ArgumentError("MarkovRegimeSymbolizer requires one latent stream"))
    n = size(latent, 2)
    n ≥ 2 || throw(ArgumentError("MarkovRegimeSymbolizer requires n ≥ 2, got $n"))
    regimes = quantize_to_symbols(vec(@view latent[1, :]),
                                  collect(1:length(symbolizer.transition_matrices)),
                                  symbolizer.regime_weights)
    symbol = weighted_sample(rng, _target_marginal(symbolizer))
    result = Vector{eltype(symbolizer.alphabet)}(undef, n)
    @inbounds for t in 1:n
        P = symbolizer.transition_matrices[regimes[t]]
        symbol = weighted_sample(rng, @view P[symbol, :])
        result[t] = symbolizer.alphabet[symbol]
    end
    return result
end

function generate(g::PropertyBasedGenerator, n::Int;
                  rng::AbstractRNG = Random.default_rng())
    width = latent_width(g.symbolizer)
    latent = generate_latent(g.source, n, width; rng)
    return symbolize(g.symbolizer, latent; rng)
end

"""
    generate_with_latent(g, n; rng) -> sequence, latent

Generate a property-based symbolic sequence and return the numerical latent
series used before symbolization.

The returned `latent` value is a `width × n` matrix, where `width` is
[`latent_width`](@ref) of the symbolizer. This helper is intended for
validation and research workflows where the numerical long-memory driver should
be diagnosed alongside the final symbolic sequence. It is additive to the common
[`generate`](@ref) contract; ordinary callers can continue to use
`generate(g, n; rng)`.

# Examples
```julia
julia> g = PropertyBasedGenerator(SpectralFGNSource(0.75),
...                               QuantileSymbolizer([:a, :b]));

julia> seq, latent = generate_with_latent(g, 16; rng = MersenneTwister(1));

julia> length(seq), size(latent)
(16, (1, 16))
```
"""
function generate_with_latent(g::PropertyBasedGenerator, n::Int;
                              rng::AbstractRNG = Random.default_rng())
    width = latent_width(g.symbolizer)
    latent = generate_latent(g.source, n, width; rng)
    return symbolize(g.symbolizer, latent; rng), latent
end

_supports_latent_width(::SpectralFGNSource, width::Int) = width ≥ 1
_supports_latent_width(::HaarLRDSource, width::Int) = width ≥ 1
_supports_latent_width(::IntermittentMapSource, width::Int) = width == 1

_target_marginal(s::QuantileSymbolizer) = copy(Float64.(s.marginal))
_target_marginal(s::ArgmaxSymbolizer) = copy(Float64.(s.marginal))
function _target_marginal(s::MarkovRegimeSymbolizer)
    p = zeros(Float64, length(s.alphabet))
    for (r, P) in enumerate(s.transition_matrices)
        p .+= s.regime_weights[r] .* stationary_distribution(P)
    end
    p ./= sum(p)
    return p
end

_symbolizer_marginal_capability(::QuantileSymbolizer) = :finite_sample
_symbolizer_marginal_capability(::ArgmaxSymbolizer) = :empirical
_symbolizer_marginal_capability(::MarkovRegimeSymbolizer) = :implied
_symbolizer_bigram_capability(::MarkovRegimeSymbolizer) = :per_regime
_symbolizer_bigram_capability(::Symbolizer) = :induced
_symbolizer_alphabet(s::Union{QuantileSymbolizer,ArgmaxSymbolizer,
                              MarkovRegimeSymbolizer}) = s.alphabet
