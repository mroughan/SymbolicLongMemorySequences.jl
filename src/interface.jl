"""
    LRDGenerator

Abstract supertype for all LRD symbol-sequence generators in S5.jl.

Concrete subtypes must implement:

    generate(g::MyGenerator, n::Int; rng::AbstractRNG = Random.default_rng()) -> Vector
"""
abstract type LRDGenerator end

"""
    LocalStructureSpec

Abstract supertype for explicit local-structure specifications.

`MarkovSpec` is the current concrete first-order specification. Future higher-order
specifications, such as sparse trigram controls, should subtype
`LocalStructureSpec` and define [`local_structure_order`](@ref).
"""
abstract type LocalStructureSpec end

"""
    MarkovSpec(alphabet, transition_matrix)

Validated first-order Markov specification over an ordered symbol alphabet.

`transition_matrix[i, j]` is the conditional probability of emitting
`alphabet[j]` after `alphabet[i]`.

# Examples
```julia
julia> spec = MarkovSpec([:a, :b], [0.9 0.1; 0.2 0.8])
MarkovSpec{Vector{Symbol}}(k=2)

julia> spec.transition_matrix[1, :]
2-element Vector{Float64}:
 0.9
 0.1
```
"""
struct MarkovSpec{A} <: LocalStructureSpec
    alphabet          :: A
    transition_matrix :: Matrix{Float64}

    function MarkovSpec{A}(alphabet::A,
                           transition_matrix::Matrix{Float64}) where {A}
        validate_alphabet(alphabet)
        k = length(alphabet)
        size(transition_matrix) == (k, k) ||
            throw(ArgumentError(
                "transition_matrix must have size ($k, $k), got $(size(transition_matrix))"))
        new{A}(alphabet, transition_matrix)
    end
end

function MarkovSpec(alphabet, transition_matrix::AbstractMatrix{<:Real})
    P = validate_transition_matrix(transition_matrix)
    MarkovSpec{typeof(alphabet)}(alphabet, P)
end

"""
    local_structure_order(spec) -> Int

Return the Markov order of a local-structure specification.

`MarkovSpec` has order 1. This function is the extension point for future
higher-order local-structure specifications; S5.jl does not currently expose a
trigram-control specification.

# Examples
```julia
julia> local_structure_order(MarkovSpec([:a, :b], [0.9 0.1; 0.2 0.8]))
1
```
"""
function local_structure_order end

local_structure_order(::MarkovSpec) = 1

function Base.show(io::IO, spec::MarkovSpec)
    print(io, "MarkovSpec{$(typeof(spec.alphabet))}(k=$(length(spec.alphabet)))")
end

"""
    ControlCapabilities

Programmatic description of a generator's user-facing control strengths.

Fields use stable symbolic levels:

- `alphabet`: `:exact`
- `marginal`: `:finite_sample`, `:empirical`, `:implied`,
  `:innovation_target`, or `:asymptotic`
- `bigram`: `:per_regime` or `:induced`
- `trigram`: `:induced`
- `lrd`: `:approximate`, `:latent_approximate`, `:finite_history`, or `:nominal`
"""
struct ControlCapabilities
    alphabet :: Symbol
    marginal :: Symbol
    bigram   :: Symbol
    trigram  :: Symbol
    lrd      :: Symbol
end

"""
    control_capabilities(g) -> ControlCapabilities

Return the declared control strengths of generator `g`.

# Examples
```julia
julia> control_capabilities(SpectralFGN(0.8, [:a, :b])).marginal
:finite_sample

julia> control_capabilities(FSS(1.5, [:a, :b])).bigram
:induced
```
"""
function control_capabilities end

"""
    generate(g, n; rng = Random.default_rng()) -> Vector

Generate a sequence of `n` symbols using LRD generator `g`.

Returns a `Vector` whose element type matches the alphabet of `g`.

# Arguments
- `g::LRDGenerator`: configured generator instance.
- `n::Int`: number of symbols to emit.

# Keyword Arguments
- `rng::AbstractRNG`: random number generator (default: `Random.default_rng()`).

# Examples
```julia
julia> g = SpectralFGN(0.8, ['a', 'b', 'c'])
julia> seq = generate(g, 1000)
julia> length(seq) == 1000
true
```
"""
function generate end
