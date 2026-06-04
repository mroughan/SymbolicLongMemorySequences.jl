"""
    LAMP(beta, alphabet [, marginal]; d = 1000)

Model-based LRD symbol-sequence generator (MB1): Linear-Additive Markov Process.

At each step the probability of the next symbol is a convex combination of
one-hot indicator vectors for the most recent `d` history symbols:

    P(Xₜ = s | Xₜ₋₁, …, Xₜ₋ᵈ) = Σⱼ wⱼ · 𝟏[Xₜ₋ⱼ = s]

with power-law weights `wⱼ ∝ j^{-(1+β)}`, so the autocovariance decays as a
power law with exponent `β`, giving Hurst parameter `H = (2−β)/2`.

# Arguments
- `beta::Real`: ACF decay exponent, `β ∈ (0, 1)`.
- `alphabet`: ordered collection of symbols.
- `marginal::AbstractVector{<:Real}`: stationary marginal (default: uniform).

# Keyword Arguments
- `d::Int = 1000`: history depth. The effective LRD range is bounded by `d`; for
  a sequence of length `n` set `d ≥ n^{1/(1+β)}` to avoid truncation artefacts.

# Complexity
O(n·d) time, O(d + n) memory.

# References
Kumar, R., Raghu, M., Sarlos, T., & Tomkins, A. (2017). Linear additive Markov
processes. *WWW '17*, 411–419.

Singh, M., Greenberg, C., & Klakow, D. (2016). The custom decay language model
for long range dependencies. *TSD*, 343–351.

# Examples
```julia
julia> g = LAMP(0.5, [:a, :b, :c]; d = 500)
LAMP{Vector{Symbol}, Vector{Float64}}(β=0.5, k=3, d=500)

julia> seq = generate(g, 5000; rng = MersenneTwister(42))
julia> length(seq) == 5000 && eltype(seq) == Symbol
true
```
"""
struct LAMP{A, M <: AbstractVector{<:Real}} <: LRDGenerator
    beta     :: Float64
    alphabet :: A
    marginal :: M
    d        :: Int
    weights  :: Vector{Float64}

    function LAMP{A, M}(beta::Float64, alphabet::A, marginal::M,
                         d::Int, weights::Vector{Float64}) where {A, M <: AbstractVector{<:Real}}
        (0.0 < beta < 1.0) ||
            throw(ArgumentError("beta must be in (0, 1), got $beta"))
        k = length(alphabet)
        length(marginal) == k ||
            throw(ArgumentError(
                "marginal length $(length(marginal)) ≠ alphabet length $k"))
        isapprox(sum(marginal), 1.0; atol = 1e-8) ||
            throw(ArgumentError("marginal must sum to 1, got $(sum(marginal))"))
        d ≥ 1 || throw(ArgumentError("d must be ≥ 1, got $d"))
        new{A, M}(beta, alphabet, marginal, d, weights)
    end
end

function LAMP(beta::Real, alphabet,
              marginal::AbstractVector{<:Real} =
                  fill(1.0 / length(alphabet), length(alphabet));
              d::Int = 1000)
    m  = Float64.(marginal)
    w  = [j^(-(1.0 + Float64(beta))) for j in 1:d]
    w ./= sum(w)
    LAMP{typeof(alphabet), typeof(m)}(Float64(beta), alphabet, m, d, w)
end

function Base.show(io::IO, g::LAMP)
    print(io, "LAMP{$(typeof(g.alphabet)), $(typeof(g.marginal))}",
          "(β=$(g.beta), k=$(length(g.alphabet)), d=$(g.d))")
end

"""
    generate(g::LAMP, n; rng) -> Vector

Generate `n` symbols from a [`LAMP`](@ref) generator using a ring-buffer history.

The ring buffer stores the `d` most recent symbol indices. The next symbol is
drawn from a probability vector formed by weighting each history element by its
power-law weight `wⱼ`.
"""
function generate(g::LAMP, n::Int; rng::AbstractRNG = Random.default_rng())
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1, got $n"))
    k = length(g.alphabet)
    d = g.d

    # Ring buffer of symbol indices (1-based).
    # head = next write position. Element j steps back is at mod1(head-j, d).
    buf  = [weighted_sample(rng, g.marginal) for _ in 1:d]
    head = 1

    q      = Vector{Float64}(undef, k)
    result = Vector{eltype(g.alphabet)}(undef, n)

    @inbounds for t in 1:n
        fill!(q, 0.0)
        for j in 1:d
            q[buf[mod1(head - j, d)]] += g.weights[j]
        end

        idx        = weighted_sample(rng, q)
        result[t]  = g.alphabet[idx]
        buf[head]  = idx
        head       = mod1(head + 1, d)
    end

    return result
end
