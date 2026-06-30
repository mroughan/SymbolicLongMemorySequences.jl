"""
    FSS(alpha, alphabet; rates = ones(k), x_min = 1.0)

Model-based LRD symbol-sequence generator (MB3): Fractal Symbol Sequence.

Each symbol is governed by an independent Pareto-distributed renewal process:
inter-arrival times `τ ~ Pareto(α, x_min)`. The output merges all symbol
streams in event-time order — the symbol with the earliest pending event is
emitted at each step.

LRD arises through heavy-tailed inter-arrival times. For tail index `α ∈ (1, 2)`,
the return-time variance is infinite, giving nominal Hurst parameter
`H = (3−α)/2`.

# Arguments
- `alpha::Real`: Pareto tail index, `α ∈ (1, 2)`.
- `alphabet`: ordered collection of symbols.

# Keyword Arguments
- `rates::AbstractVector{<:Real}`: per-symbol base arrival rates. Symbol `i`
  appears with long-run frequency proportional to `rates[i]`. Default: uniform.
- `x_min::Real = 1.0`: Pareto scale parameter (minimum inter-arrival time).

# Complexity
O(n·k) time, O(k + n) memory (`k` = alphabet size).

# Notes
Because each symbol stream is independent, joint symbol statistics (bigrams, etc.)
cannot be prescribed independently of the marginal.

*Missing-scales pitfall* (Roughan, Yates & Veitch 1999): if the mean inter-arrival
time `x_min · α/(α−1) / rateᵢ` is large relative to `n`, the observable LRD
scale range is reduced. Keep rates such that each symbol appears O(√n) or more
times.

# References
Lowen, S. B., & Teich, M. C. (1995). Estimation and simulation of fractal
stochastic point processes. *Fractals* 3(1), 183–210.

Roughan, M., Yates, J., & Veitch, D. (1999). The mystery of the missing scales:
pitfalls in the use of fractal renewal processes. *Applications of Heavy Tailed
Distributions in Economics, Engineering and Statistics*.

# Examples
```julia
julia> g = FSS(1.4, [:a, :b, :c])
FSS{Vector{Symbol}, Vector{Float64}}(α=1.4, H≈0.8, k=3)

julia> seq = generate(g, 5000; rng = MersenneTwister(7))
julia> length(seq) == 5000 && eltype(seq) == Symbol
true
```
"""
struct FSS{A, R <: AbstractVector{<:Real}} <: LRDGenerator
    alpha    :: Float64
    alphabet :: A
    rates    :: R
    x_min    :: Float64

    function FSS{A, R}(alpha::Float64, alphabet::A, rates::R,
                        x_min::Float64) where {A, R <: AbstractVector{<:Real}}
        (1.0 < alpha < 2.0) ||
            throw(ArgumentError("alpha must be in (1, 2), got $alpha"))
        validate_alphabet(alphabet)
        k = length(alphabet)
        length(rates) == k ||
            throw(ArgumentError(
                "rates length $(length(rates)) ≠ alphabet length $k"))
        all(isfinite, rates) ||
            throw(ArgumentError("rates must contain only finite values"))
        all(>(0), rates) ||
            throw(ArgumentError("rates must be positive"))
        x_min > 0 ||
            throw(ArgumentError("x_min must be positive, got $x_min"))
        isfinite(x_min) ||
            throw(ArgumentError("x_min must be finite, got $x_min"))
        new{A, R}(alpha, alphabet, rates, x_min)
    end
end

function FSS(alpha::Real, alphabet;
             rates::AbstractVector{<:Real} = fill(1.0, length(alphabet)),
             x_min::Real = 1.0)
    r = validate_positive_vector(rates, "rates")
    FSS{typeof(alphabet), typeof(r)}(Float64(alpha), alphabet, r, Float64(x_min))
end

function Base.show(io::IO, g::FSS)
    H = round((3 - g.alpha) / 2; digits = 4)
    print(io, "FSS{$(typeof(g.alphabet)), $(typeof(g.rates))}",
          "(α=$(g.alpha), H≈$H, k=$(length(g.alphabet)))")
end

"""
    generate(g::FSS, n; rng) -> Vector

Generate `n` symbols from a [`FSS`](@ref) generator via a Pareto renewal merge.

Each symbol maintains an independent clock advanced by a Pareto-distributed
inter-arrival time after each emission. At each step the symbol with the smallest
clock value is emitted.
"""
function generate(g::FSS, n::Int; rng::AbstractRNG = Random.default_rng())
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1, got $n"))
    k = length(g.alphabet)

    dist = Pareto(g.alpha, g.x_min)

    # Initialise per-symbol clocks with one draw each
    times  = [rand(rng, dist) / g.rates[i] for i in 1:k]
    result = Vector{eltype(g.alphabet)}(undef, n)

    @inbounds for t in 1:n
        idx        = argmin(times)
        result[t]  = g.alphabet[idx]
        times[idx] += rand(rng, dist) / g.rates[idx]
    end

    return result
end

"""
    _pareto_sample(rng, alpha, x_min) -> Float64

Draw from a Pareto distribution with shape `alpha` and scale `x_min` using
Distributions.jl.

# Examples
```julia
julia> SymbolicLongMemorySequences._pareto_sample(MersenneTwister(1), 1.4, 1.0) >= 1.0
true
```
"""
@inline function _pareto_sample(rng::AbstractRNG, alpha::Float64, x_min::Float64)
    return rand(rng, Pareto(alpha, x_min))
end
