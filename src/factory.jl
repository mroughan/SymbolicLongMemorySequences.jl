"""
    ParameterInfo

Metadata for one keyword accepted by [`make_generator`](@ref).

Fields:
- `name`: keyword name.
- `kind`: currently `:keyword`; reserved for future input categories.
- `default`: factory default value.
- `domain`: concise domain or accepted-value summary.
- `description`: short human-readable explanation.

# Examples
```julia
julia> p = method_parameters(:PB1)[1]
ParameterInfo(name=H, default=0.8)

julia> p.domain
"0.5 < H < 1"
```
"""
struct ParameterInfo
    name        :: Symbol
    kind        :: Symbol
    default     :: Any
    domain      :: String
    description :: String
end

function Base.show(io::IO, p::ParameterInfo)
    print(io, "ParameterInfo(name=$(p.name), default=$(p.default))")
end

"""
    MethodInfo

Metadata returned by [`method_info`](@ref) for one SymbolicLongMemorySequences synthesis method.

Fields:
- `id`: stable method identifier, such as `:PB1` or `:MB5`.
- `family`: `:property_based` or `:model_based`.
- `type_name`: exported generator type name.
- `defaults`: standard-case keyword defaults used by [`make_generator`](@ref).
- `parameters`: keyword metadata for the factory inputs accepted by this method.
- `standard_cases`: named construction presets accepted by [`make_generator`](@ref).
- `description`: short human-readable summary.

# Examples
```julia
julia> info = method_info(:PB1)
MethodInfo(id=PB1, type=SpectralFGN)

julia> info.defaults.H
0.8

julia> method_parameters(:PB1)[1].name
:H
```
"""
struct MethodInfo
    id             :: Symbol
    family         :: Symbol
    type_name      :: Symbol
    defaults       :: NamedTuple
    parameters     :: Tuple
    standard_cases :: Tuple
    description    :: String
end

function Base.show(io::IO, info::MethodInfo)
    print(io, "MethodInfo(id=$(info.id), type=$(info.type_name))")
end

const _METHOD_IDS = (:PB1, :PB2, :PB3, :PB4, :MB1a, :MB1b, :MB1c,
                     :MB2, :MB3, :MB4, :MB5)

const _METHOD_ALIASES = Dict{Symbol, Symbol}(
    :PB1 => :PB1, :SpectralFGN => :PB1,
    :PB2 => :PB2, :LGCM => :PB2,
    :PB3 => :PB3, :WaveletMarkov => :PB3,
    :PB4 => :PB4, :IntermittentMapSymbols => :PB4,
    :MB1 => :MB1a, :MB1a => :MB1a, :LAMP => :MB1a,
    :MB1b => :MB1b, :DyadicLAMP => :MB1b,
    :MB1c => :MB1c, :CalibratedAdditiveMarkov => :MB1c,
    :MB2 => :MB2, :OnOffMarkov => :MB2,
    :MB3 => :MB3, :FSS => :MB3,
    :MB4 => :MB4, :HawkesSymbol => :MB4,
    :MB5 => :MB5, :DuplicationMutation => :MB5,
)

function _canonical_method_id(id)
    sym = Symbol(id)
    haskey(_METHOD_ALIASES, sym) ||
        throw(ArgumentError("unknown SymbolicLongMemorySequences method id or type name: $id"))
    return _METHOD_ALIASES[sym]
end

"""
    method_ids(; family = :all) -> Tuple{Vararg{Symbol}}

Return the stable method identifiers accepted by [`make_generator`](@ref).

Use `family = :property_based` or `family = :model_based` to filter the list.

# Examples
```julia
julia> method_ids()[1:3]
(:PB1, :PB2, :PB3)

julia> :MB5 in method_ids(family = :model_based)
true
```
"""
function method_ids(; family::Symbol = :all)
    if family === :all
        return _METHOD_IDS
    elseif family === :property_based
        return (:PB1, :PB2, :PB3, :PB4)
    elseif family === :model_based
        return (:MB1a, :MB1b, :MB1c, :MB2, :MB3, :MB4, :MB5)
    else
        throw(ArgumentError("family must be :all, :property_based, or :model_based"))
    end
end

function _param(name::Symbol, default, domain::AbstractString,
                description::AbstractString)
    return ParameterInfo(name, :keyword, default, String(domain), String(description))
end

function _common_marginal(default = :uniform)
    return _param(:marginal, default,
                  "`:uniform` or a probability vector with one entry per symbol",
                  "Target marginal distribution when the method has a marginal-control knob.")
end

function _common_case(default, accepted)
    return _param(:case, default, accepted,
                  "Named standard construction case used when detailed matrices are omitted.")
end

function _method_info_table()
    (
        PB1 = MethodInfo(:PB1, :property_based, :SpectralFGN,
                         (; H = 0.8, marginal = :uniform),
                         (_param(:H, 0.8, "0.5 < H < 1",
                                 "Nominal Hurst parameter for the latent fGn-like process."),
                          _common_marginal()),
                         (:standard,),
                         "Approximate spectral fGn followed by rank quantization."),
        PB2 = MethodInfo(:PB2, :property_based, :LGCM,
                         (; H = 0.8, marginal = :uniform,
                            calibration_iters = 25, calibration_rate = 0.7),
                         (_param(:H, 0.8, "0.5 < H < 1",
                                 "Nominal Hurst parameter for each latent Gaussian channel."),
                          _common_marginal(),
                          _param(:calibration_iters, 25, "integer >= 0",
                                 "Number of empirical marginal-calibration iterations."),
                          _param(:calibration_rate, 0.7, "0 <= calibration_rate <= 1",
                                 "Damping applied during empirical threshold calibration.")),
                         (:standard,),
                         "Latent Gaussian categorical model with marginal calibration."),
        PB3 = MethodInfo(:PB3, :property_based, :WaveletMarkov,
                         (; H = 0.8, marginal = :uniform, transition_matrices = nothing,
                            regime_weights = nothing, cascade_depth = 0,
                            driver = :spectral, case = :persistent_regimes),
                         (_param(:H, 0.8, "0.5 < H < 1",
                                 "Nominal Hurst parameter for the latent regime driver."),
                          _common_marginal(),
                          _param(:transition_matrices, nothing,
                                 "`nothing` or a vector of stochastic matrices",
                                 "Regime-specific Markov transition matrices."),
                          _param(:regime_weights, nothing,
                                 "`nothing` or a probability vector",
                                 "Stationary regime weights used to choose latent regimes."),
                          _param(:cascade_depth, 0, "integer >= 0",
                                 "Optional Haar-cascade refinement depth for the driver."),
                          _param(:driver, :spectral, "`:spectral` or `:haar`",
                                 "Numerical LRD driver used before Markov symbolization."),
                          _common_case(:persistent_regimes,
                                       "`:persistent_regimes` or `:iid_regimes`")),
                         (:persistent_regimes, :iid_regimes),
                         "Latent LRD regime driver selecting Markov transition matrices."),
        PB4 = MethodInfo(:PB4, :property_based, :IntermittentMapSymbols,
                         (; z = 1.6, marginal = :uniform, burnin = 1000),
                         (_param(:z, 1.6, "1 < z < 2",
                                 "Intermittency parameter controlling nominal long memory."),
                          _common_marginal(),
                          _param(:burnin, 1000, "integer >= 0",
                                 "Number of latent map iterations discarded before sampling.")),
                         (:standard,),
                         "Intermittent latent map followed by rank quantization."),
        MB1a = MethodInfo(:MB1a, :model_based, :LAMP,
                          (; beta = 0.5, marginal = :uniform, d = 1000,
                             epsilon = 0.02, transition_matrix = nothing,
                             repeat_probability = 0.9, case = :repeat),
                          (_param(:beta, 0.5, "0 < beta < 1",
                                  "Power-law memory exponent for finite-history weights."),
                           _common_marginal(),
                           _param(:d, 1000, "integer >= 1",
                                  "Explicit finite history cutoff."),
                           _param(:epsilon, 0.02, "0 <= epsilon <= 1",
                                  "Mixture weight for the innovation component."),
                           _param(:transition_matrix, nothing,
                                  "`nothing` or a stochastic matrix",
                                  "Local transition matrix used by the LAMP mechanism."),
                           _param(:repeat_probability, 0.9, "0 <= repeat_probability <= 1",
                                  "Persistence level for the standard repeat case."),
                           _common_case(:repeat, "`:repeat`, `:persistent`, or `:iid`")),
                          (:repeat, :iid),
                          "Exact finite-history LAMP with power-law history weights."),
        MB1b = MethodInfo(:MB1b, :model_based, :DyadicLAMP,
                          (; beta = 0.5, marginal = :uniform, d = 100_000,
                             epsilon = 0.02, transition_matrix = nothing,
                             repeat_probability = 0.9, case = :repeat),
                          (_param(:beta, 0.5, "0 < beta < 1",
                                  "Power-law memory exponent for dyadic history weights."),
                           _common_marginal(),
                           _param(:d, 100_000, "integer >= 1",
                                  "Maximum represented history depth."),
                           _param(:epsilon, 0.02, "0 <= epsilon <= 1",
                                  "Mixture weight for the innovation component."),
                           _param(:transition_matrix, nothing,
                                  "`nothing` or a stochastic matrix",
                                  "Local transition matrix used by the dyadic LAMP mechanism."),
                           _param(:repeat_probability, 0.9, "0 <= repeat_probability <= 1",
                                  "Persistence level for the standard repeat case."),
                           _common_case(:repeat, "`:repeat`, `:persistent`, or `:iid`")),
                          (:repeat, :iid),
                          "Scalable dyadic-bucket approximation to LAMP."),
        MB1c = MethodInfo(:MB1c, :model_based, :CalibratedAdditiveMarkov,
                          (; beta = 0.5, marginal = :uniform, d = 1000,
                             strength = 0.8, case = :standard),
                          (_param(:beta, 0.5, "0 < beta < 1",
                                  "Power-law exponent for the centered additive memory kernel."),
                           _common_marginal(),
                           _param(:d, 1000, "integer >= 1",
                                  "Explicit finite memory cutoff."),
                           _param(:strength, 0.8, "0 <= strength <= 1",
                                  "Dependence strength in the standard case."),
                           _common_case(:standard, "`:standard` or `:iid`")),
                          (:standard, :iid),
                          "Centered additive Markov memory function."),
        MB2 = MethodInfo(:MB2, :model_based, :OnOffMarkov,
                         (; alpha = 1.4, marginal = :uniform,
                            transition_matrices = nothing, switching_matrix = nothing,
                            L_min = 50.0, case = :persistent_regimes),
                         (_param(:alpha, 1.4, "1 < alpha < 2",
                                 "Tail exponent for regime holding times."),
                          _common_marginal(),
                          _param(:transition_matrices, nothing,
                                 "`nothing` or a vector of stochastic matrices",
                                 "Regime-specific Markov transition matrices."),
                          _param(:switching_matrix, nothing,
                                 "`nothing` or a stochastic matrix",
                                 "Transition matrix for regime switches."),
                          _param(:L_min, 50.0, "L_min > 0",
                                 "Minimum scale for heavy-tailed holding times."),
                          _common_case(:persistent_regimes,
                                       "`:persistent_regimes` or `:iid_regimes`")),
                         (:persistent_regimes, :iid_regimes),
                         "Heavy-tailed regime-switching Markov chain."),
        MB3 = MethodInfo(:MB3, :model_based, :FSS,
                         (; alpha = 1.4, marginal = nothing,
                            rates = :uniform, x_min = 1.0),
                         (_param(:alpha, 1.4, "1 < alpha < 2",
                                 "Tail exponent for renewal durations."),
                          _param(:marginal, nothing,
                                 "`nothing` or a positive vector used as renewal rates",
                                 "Compatibility alias for `rates`; not a direct marginal guarantee."),
                          _param(:rates, :uniform,
                                 "`:uniform` or a positive vector with one entry per symbol",
                                 "Relative symbol renewal rates."),
                          _param(:x_min, 1.0, "x_min > 0",
                                 "Minimum Pareto renewal duration.")),
                         (:standard,),
                         "Fractal symbol sequence via independent Pareto renewals."),
        MB4 = MethodInfo(:MB4, :model_based, :HawkesSymbol,
                         (; beta = 0.6, marginal = nothing, baseline = :uniform,
                            excitation = :identity, d = 1000, c = 1.0),
                         (_param(:beta, 0.6, "0 < beta < 1",
                                 "Power-law exponent for finite-history excitation decay."),
                          _param(:marginal, nothing,
                                 "`nothing` or a positive vector used as baseline intensities",
                                 "Compatibility alias for `baseline`; not exact marginal control."),
                          _param(:baseline, :uniform,
                                 "`:uniform` or a positive vector with one entry per symbol",
                                 "Baseline symbol intensities."),
                          _param(:excitation, :identity,
                                 "`:identity` or a nonnegative matrix",
                                 "Symbol-to-symbol self-excitation matrix."),
                          _param(:d, 1000, "integer >= 1",
                                 "Explicit finite history cutoff."),
                          _param(:c, 1.0, "c > 0",
                                 "Scale parameter for the excitation kernel.")),
                         (:identity_excitation,),
                         "Finite-history Hawkes-style symbolic self-excitation."),
        MB5 = MethodInfo(:MB5, :model_based, :DuplicationMutation,
                         (; alpha = 1.4, marginal = :uniform,
                            mutation_probability = 0.02, seed_length = 64,
                            max_block_length = 4096),
                         (_param(:alpha, 1.4, "1 < alpha < 2",
                                 "Tail exponent for power-law copy distances."),
                          _common_marginal(),
                          _param(:mutation_probability, 0.02,
                                 "0 <= mutation_probability <= 1",
                                 "Probability of replacing copied symbols with innovations."),
                          _param(:seed_length, 64, "integer >= 1",
                                 "Initial iid prefix length."),
                          _param(:max_block_length, 4096, "integer >= 1",
                                 "Maximum contiguous copy block length.")),
                         (:standard,),
                         "Power-law lag copy-and-mutate symbolic growth."),
    )
end

"""
    method_info(id) -> MethodInfo
    method_info() -> Tuple{Vararg{MethodInfo}}

Return metadata and standard defaults for one method accepted by
[`make_generator`](@ref). `id` may be a method identifier such as `:PB1` or an
exported type name such as `:SpectralFGN`. With no argument, return metadata for
all methods in [`method_ids`](@ref) order.

# Examples
```julia
julia> method_info(:MB5).type_name
:DuplicationMutation

julia> method_info(:SpectralFGN).id
:PB1

julia> length(method_info()) == length(method_ids())
true
```
"""
function method_info(id)
    canon = _canonical_method_id(id)
    return getproperty(_method_info_table(), canon)
end

function method_info()
    table = _method_info_table()
    return tuple((getproperty(table, id) for id in _METHOD_IDS)...)
end

"""
    method_parameters(id) -> Tuple{Vararg{ParameterInfo}}

Return keyword metadata for the factory inputs accepted by one method.

This is a discovery helper for user interfaces, examples, and benchmark grids.
All methods still require the positional `alphabet` input to
[`make_generator`](@ref), and sequence length `n` is supplied later to
[`generate`](@ref).

# Examples
```julia
julia> first(method_parameters(:PB1)).name
:H

julia> any(p -> p.name === :mutation_probability, method_parameters(:MB5))
true
```
"""
function method_parameters(id)
    return method_info(id).parameters
end

function _uniform_probs(alphabet)
    validate_alphabet(alphabet)
    return fill(1.0 / length(alphabet), length(alphabet))
end

function _resolve_marginal(alphabet, marginal)
    marginal === :uniform && return _uniform_probs(alphabet)
    return validate_probability_vector(marginal, "marginal")
end

function _resolve_rates(alphabet, rates)
    rates === :uniform && return fill(1.0, length(alphabet))
    return validate_positive_vector(rates, "rates")
end

function _persistent_matrix(p::AbstractVector{<:Real}, repeat_probability::Real)
    0.0 <= repeat_probability <= 1.0 ||
        throw(ArgumentError("repeat_probability must be in [0, 1]"))
    k = length(p)
    P = Matrix{Float64}(undef, k, k)
    @inbounds for i in 1:k, j in 1:k
        P[i, j] = (i == j ? repeat_probability : 0.0) +
                  (1 - repeat_probability) * p[j]
    end
    return P
end

function _iid_matrix(p::AbstractVector{<:Real})
    k = length(p)
    P = Matrix{Float64}(undef, k, k)
    @inbounds for i in 1:k, j in 1:k
        P[i, j] = p[j]
    end
    return P
end

function _standard_regime_matrices(alphabet, marginal, case::Symbol)
    p = _resolve_marginal(alphabet, marginal)
    if case === :iid_regimes
        return [_iid_matrix(p), _iid_matrix(p)]
    elseif case === :persistent_regimes
        return [_iid_matrix(p), _persistent_matrix(p, 0.92)]
    else
        throw(ArgumentError("case must be :persistent_regimes or :iid_regimes"))
    end
end

function _standard_switching_matrix(nregimes::Int)
    nregimes >= 2 || throw(ArgumentError("at least two regimes are required"))
    Q = fill(1.0 / (nregimes - 1), nregimes, nregimes)
    @inbounds for i in 1:nregimes
        Q[i, i] = 0.0
    end
    return Q
end

function _standard_excitation(alphabet, excitation)
    k = length(alphabet)
    excitation === :identity && return Matrix{Float64}(I, k, k)
    return Matrix{Float64}(excitation)
end

"""
    make_generator(id, alphabet; kwargs...) -> LRDGenerator

Construct a standard SymbolicLongMemorySequences generator by method identifier.

This is a convenience API for common cases. It does not replace the explicit
scientific constructors: use `method_info(id).defaults` to inspect default
parameters, then pass keyword overrides as needed. `id` may be a method id
(`:PB1`, `:MB1c`) or a type name (`:SpectralFGN`, `:DuplicationMutation`).

# Common Keywords
- `marginal = :uniform`: target marginal where the method has one.
- `case`: standard preset for methods that need local/regime structure.

# Examples
```julia
julia> g = make_generator(:PB1, [:a, :b]; H = 0.75)
SpectralFGN{Vector{Symbol}, Vector{Float64}}(H=0.75, k=2)

julia> generate(g, 4; rng = MersenneTwister(1)) isa Vector{Symbol}
true

julia> make_generator(:MB5, ['A', 'C']; alpha = 1.4, max_block_length = 128)
DuplicationMutation{Vector{Char}, Vector{Float64}}(α=1.4, k=2, μ=0.02, seed=64, max_block=128)
```
"""
function make_generator(id, alphabet; kwargs...)
    canon = _canonical_method_id(id)
    try
        return _make_generator(Val(canon), alphabet; kwargs...)
    catch err
        err isa MethodError ||
            rethrow()
        throw(ArgumentError("invalid keyword argument for method $canon"))
    end
end

function _make_generator(::Val{:PB1}, alphabet; H::Real = 0.8,
                         marginal = :uniform)
    return SpectralFGN(H, alphabet, _resolve_marginal(alphabet, marginal))
end

function _make_generator(::Val{:PB2}, alphabet; H::Real = 0.8,
                         marginal = :uniform,
                         calibration_iters::Int = 25,
                         calibration_rate::Real = 0.7)
    return LGCM(H, alphabet, _resolve_marginal(alphabet, marginal);
                calibration_iters, calibration_rate)
end

function _make_generator(::Val{:PB3}, alphabet; H::Real = 0.8,
                         marginal = :uniform,
                         transition_matrices = nothing,
                         regime_weights = nothing,
                         cascade_depth::Int = 0,
                         driver::Symbol = :spectral,
                         case::Symbol = :persistent_regimes)
    Ps = transition_matrices === nothing ?
         _standard_regime_matrices(alphabet, marginal, case) :
         [validate_transition_matrix(P, "transition_matrices[$i]")
          for (i, P) in enumerate(transition_matrices)]
    weights = regime_weights === nothing ?
              fill(1.0 / length(Ps), length(Ps)) :
              validate_probability_vector(regime_weights, "regime_weights")
    return WaveletMarkov(H, alphabet, Ps; regime_weights = weights,
                         cascade_depth, driver)
end

function _make_generator(::Val{:PB4}, alphabet; z::Real = 1.6,
                         marginal = :uniform,
                         burnin::Int = 1000)
    return IntermittentMapSymbols(z, alphabet, _resolve_marginal(alphabet, marginal);
                                  burnin)
end

function _make_generator(::Val{:MB1a}, alphabet; beta::Real = 0.5,
                         marginal = :uniform, d::Int = 1000,
                         epsilon::Real = 0.02,
                         transition_matrix = nothing,
                         repeat_probability::Real = 0.9,
                         case::Symbol = :repeat)
    p = _resolve_marginal(alphabet, marginal)
    P = transition_matrix === nothing ?
        (case === :iid ? _iid_matrix(p) :
         (case === :repeat || case === :persistent) ?
         lamp_repeat_transition(p; repeat_probability) :
         throw(ArgumentError("case must be :repeat, :persistent, or :iid"))) :
        validate_transition_matrix(transition_matrix)
    return LAMP(beta, alphabet, p; d, epsilon, transition_matrix = P)
end

function _make_generator(::Val{:MB1b}, alphabet; beta::Real = 0.5,
                         marginal = :uniform, d::Int = 100_000,
                         epsilon::Real = 0.02,
                         transition_matrix = nothing,
                         repeat_probability::Real = 0.9,
                         case::Symbol = :repeat)
    p = _resolve_marginal(alphabet, marginal)
    P = transition_matrix === nothing ?
        (case === :iid ? _iid_matrix(p) :
         (case === :repeat || case === :persistent) ?
         lamp_repeat_transition(p; repeat_probability) :
         throw(ArgumentError("case must be :repeat, :persistent, or :iid"))) :
        validate_transition_matrix(transition_matrix)
    return DyadicLAMP(beta, alphabet, p; d, epsilon, transition_matrix = P)
end

function _make_generator(::Val{:MB1c}, alphabet; beta::Real = 0.5,
                         marginal = :uniform, d::Int = 1000,
                         strength::Real = 0.8,
                         case::Symbol = :standard)
    actual_strength = case === :iid ? 0.0 :
                      case === :standard ? strength :
                      throw(ArgumentError("case must be :standard or :iid"))
    return CalibratedAdditiveMarkov(beta, alphabet,
                                    _resolve_marginal(alphabet, marginal);
                                    d, strength = actual_strength)
end

function _make_generator(::Val{:MB2}, alphabet; alpha::Real = 1.4,
                         marginal = :uniform,
                         transition_matrices = nothing,
                         switching_matrix = nothing,
                         L_min::Real = 50.0,
                         case::Symbol = :persistent_regimes)
    Ps = transition_matrices === nothing ?
         _standard_regime_matrices(alphabet, marginal, case) :
         [validate_transition_matrix(P, "transition_matrices[$i]")
          for (i, P) in enumerate(transition_matrices)]
    Q = switching_matrix === nothing ?
        _standard_switching_matrix(length(Ps)) :
        validate_transition_matrix(switching_matrix, "switching_matrix")
    return OnOffMarkov(alpha, alphabet, Ps, Q; L_min)
end

function _make_generator(::Val{:MB3}, alphabet; alpha::Real = 1.4,
                         marginal = nothing, rates = :uniform,
                         x_min::Real = 1.0)
    actual_rates = marginal === nothing ? rates : marginal
    return FSS(alpha, alphabet; rates = _resolve_rates(alphabet, actual_rates),
               x_min)
end

function _make_generator(::Val{:MB4}, alphabet; beta::Real = 0.6,
                         marginal = nothing,
                         baseline = :uniform,
                         excitation = :identity,
                         d::Int = 1000,
                         c::Real = 1.0)
    actual_baseline = marginal === nothing ? baseline : marginal
    b = actual_baseline === :uniform ? fill(1.0, length(alphabet)) :
        validate_positive_vector(actual_baseline, "baseline")
    E = _standard_excitation(alphabet, excitation)
    return HawkesSymbol(beta, alphabet; baseline = b, excitation = E, d, c)
end

function _make_generator(::Val{:MB5}, alphabet; alpha::Real = 1.4,
                         marginal = :uniform,
                         mutation_probability::Real = 0.02,
                         seed_length::Int = 64,
                         max_block_length::Int = 4096)
    return DuplicationMutation(alpha, alphabet, _resolve_marginal(alphabet, marginal);
                               mutation_probability, seed_length, max_block_length)
end
