"""
    LGCM(H, alphabet [, marginal]; calibration_iters = 25, calibration_rate = 0.7)

Property-based LRD symbol-sequence generator (PB2): Latent Gaussian Categorical
Model.

For each symbol, `LGCM` generates a latent fGn stream with Hurst parameter `H`.
At each time step the emitted symbol is the argmax of the latent streams after
adding per-symbol offsets. The offsets are calibrated on the generated latent
sample so the output marginal is close to the requested `marginal`.

# Arguments
- `H::Real`: Hurst parameter, `H ∈ (0.5, 1.0)`.
- `alphabet`: ordered collection of unique symbols.
- `marginal::AbstractVector{<:Real}`: target marginal probabilities.

# Keyword Arguments
- `calibration_iters::Int = 25`: number of mean-offset calibration passes.
- `calibration_rate::Real = 0.7`: step size for the log-ratio calibration update.

# Complexity
O(calibration_iters · n · k) time and O(n · k) memory.

# Notes
The marginal calibration is empirical: it adjusts offsets for the realised latent
sample rather than solving the exact multivariate Gaussian choice probabilities.
It gives practical marginal control while preserving the latent argmax mechanism.
For exact finite-sample marginal counts, use [`SpectralFGN`](@ref).
"""
struct LGCM{A, M <: AbstractVector{<:Real}} <: LRDGenerator
    H                 :: Float64
    alphabet          :: A
    marginal          :: M
    calibration_iters :: Int
    calibration_rate  :: Float64

    function LGCM{A, M}(H::Float64, alphabet::A, marginal::M,
                        calibration_iters::Int,
                        calibration_rate::Float64) where {A, M <: AbstractVector{<:Real}}
        (0.5 < H < 1.0) ||
            throw(ArgumentError("H must be in (0.5, 1.0), got $H"))
        validate_alphabet(alphabet)
        k = length(alphabet)
        length(marginal) == k ||
            throw(ArgumentError(
                "marginal length $(length(marginal)) ≠ alphabet length $k"))
        calibration_iters ≥ 0 ||
            throw(ArgumentError("calibration_iters must be non-negative"))
        calibration_rate > 0 && isfinite(calibration_rate) ||
            throw(ArgumentError("calibration_rate must be positive and finite"))
        new{A, M}(H, alphabet, marginal, calibration_iters, calibration_rate)
    end
end

function LGCM(H::Real, alphabet,
              marginal::AbstractVector{<:Real} =
                  fill(1.0 / length(alphabet), length(alphabet));
              calibration_iters::Int = 25,
              calibration_rate::Real = 0.7)
    m = validate_probability_vector(marginal, "marginal")
    LGCM{typeof(alphabet), typeof(m)}(Float64(H), alphabet, m,
                                      calibration_iters, Float64(calibration_rate))
end

function Base.show(io::IO, g::LGCM)
    print(io, "LGCM{$(typeof(g.alphabet)), $(typeof(g.marginal))}",
          "(H=$(g.H), k=$(length(g.alphabet)))")
end

"""
    generate(g::LGCM, n; rng) -> Vector

Generate `n` symbols from a [`LGCM`](@ref) generator.
"""
function generate(g::LGCM, n::Int; rng::AbstractRNG = Random.default_rng())
    n ≥ 4 || throw(ArgumentError(
        "LGCM requires n ≥ 4 (for latent fGn synthesis), got $n"))

    k = length(g.alphabet)
    latent = Matrix{Float64}(undef, k, n)
    @inbounds for i in 1:k
        latent[i, :] .= _fgn_spectral(n, g.H, rng)
    end

    offsets = log.(g.marginal .+ eps(Float64))
    _calibrate_lgcm_offsets!(offsets, latent, g.marginal,
                             g.calibration_iters, g.calibration_rate)

    result = Vector{eltype(g.alphabet)}(undef, n)
    @inbounds for t in 1:n
        idx = _argmax_with_offsets(latent, offsets, t)
        result[t] = g.alphabet[idx]
    end
    return result
end

function _calibrate_lgcm_offsets!(offsets::Vector{Float64},
                                  latent::Matrix{Float64},
                                  target::Vector{Float64},
                                  iters::Int,
                                  rate::Float64)
    k, n = size(latent)
    counts = zeros(Int, k)
    for _ in 1:iters
        fill!(counts, 0)
        @inbounds for t in 1:n
            counts[_argmax_with_offsets(latent, offsets, t)] += 1
        end
        observed = (counts .+ 0.5) ./ (n + 0.5 * k)
        offsets .+= rate .* (log.(target .+ eps(Float64)) .- log.(observed))
        offsets .-= mean(offsets)
    end
    return offsets
end

@inline function _argmax_with_offsets(latent::Matrix{Float64},
                                      offsets::Vector{Float64},
                                      t::Int)
    best_i = 1
    best_v = latent[1, t] + offsets[1]
    @inbounds for i in 2:size(latent, 1)
        v = latent[i, t] + offsets[i]
        if v > best_v
            best_i = i
            best_v = v
        end
    end
    return best_i
end
