"""
    OnOffMarkov(alpha, alphabet, transition_matrices, switching_matrix; L_min = 1.0)

Model-based LRD symbol-sequence generator (MB2): Heavy-tailed On/Off
doubly-stochastic Markov chain.

The generator alternates between regimes. Each regime has its own Markov
transition matrix over `alphabet`; regime sojourn lengths are Pareto-distributed
with tail index `alpha`. A row-stochastic `switching_matrix` controls which
regime follows the current one after a sojourn ends.

# Arguments
- `alpha::Real`: Pareto tail index, `alpha ∈ (1, 2)`, with nominal
  `H = (3 - alpha) / 2`.
- `alphabet`: ordered collection of unique symbols.
- `transition_matrices`: vector of row-stochastic `k × k` matrices, one per regime.
- `switching_matrix`: row-stochastic `R × R` regime transition matrix.

# Keyword Arguments
- `L_min::Real = 1.0`: Pareto scale parameter for regime sojourns.

# Complexity
O(n · k) time with the current sequential sampler and O(n + R · k²) memory.

# Notes
This method is the natural first implementation for user-specified bigram
structure: each regime has an explicit Markov transition matrix. Aggregate
marginals and bigrams depend on regime occupancy, switching dynamics, and the
per-regime stationary distributions.
"""
struct OnOffMarkov{A} <: LRDGenerator
    alpha               :: Float64
    alphabet            :: A
    transition_matrices :: Vector{Matrix{Float64}}
    switching_matrix    :: Matrix{Float64}
    L_min               :: Float64

    function OnOffMarkov{A}(alpha::Float64, alphabet::A,
                            transition_matrices::Vector{Matrix{Float64}},
                            switching_matrix::Matrix{Float64},
                            L_min::Float64) where {A}
        (1.0 < alpha < 2.0) ||
            throw(ArgumentError("alpha must be in (1, 2), got $alpha"))
        validate_alphabet(alphabet)
        k = length(alphabet)
        isempty(transition_matrices) &&
            throw(ArgumentError("transition_matrices must be non-empty"))
        all(size(P) == (k, k) for P in transition_matrices) ||
            throw(ArgumentError("each transition matrix must have size ($k, $k)"))
        R = length(transition_matrices)
        size(switching_matrix) == (R, R) ||
            throw(ArgumentError(
                "switching_matrix must have size ($R, $R), got $(size(switching_matrix))"))
        L_min > 0 && isfinite(L_min) ||
            throw(ArgumentError("L_min must be positive and finite, got $L_min"))
        new{A}(alpha, alphabet, transition_matrices, switching_matrix, L_min)
    end
end

function OnOffMarkov(alpha::Real, alphabet,
                     transition_matrices::AbstractVector{<:AbstractMatrix{<:Real}},
                     switching_matrix::AbstractMatrix{<:Real};
                     L_min::Real = 1.0)
    Ps = [validate_transition_matrix(P, "transition_matrices[$i]")
          for (i, P) in enumerate(transition_matrices)]
    Q = validate_transition_matrix(switching_matrix, "switching_matrix")
    OnOffMarkov{typeof(alphabet)}(Float64(alpha), alphabet, Ps, Q, Float64(L_min))
end

function Base.show(io::IO, g::OnOffMarkov)
    H = round((3 - g.alpha) / 2; digits = 4)
    print(io, "OnOffMarkov{$(typeof(g.alphabet))}",
          "(α=$(g.alpha), H≈$H, k=$(length(g.alphabet)), R=$(length(g.transition_matrices)))")
end

"""
    generate(g::OnOffMarkov, n; rng) -> Vector

Generate `n` symbols from an [`OnOffMarkov`](@ref) generator.
"""
function generate(g::OnOffMarkov, n::Int; rng::AbstractRNG = Random.default_rng())
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1, got $n"))

    regime_dist = Pareto(g.alpha, g.L_min)
    regime_stationary = stationary_distribution(g.switching_matrix)
    regime = weighted_sample(rng, regime_stationary)
    symbol = weighted_sample(rng, target_marginal(g))

    result = Vector{eltype(g.alphabet)}(undef, n)
    t = 1
    while t ≤ n
        sojourn = max(1, ceil(Int, rand(rng, regime_dist)))
        P = g.transition_matrices[regime]
        stop = min(n, t + sojourn - 1)
        @inbounds while t ≤ stop
            symbol = weighted_sample(rng, @view P[symbol, :])
            result[t] = g.alphabet[symbol]
            t += 1
        end
        regime = weighted_sample(rng, @view g.switching_matrix[regime, :])
    end

    return result
end
