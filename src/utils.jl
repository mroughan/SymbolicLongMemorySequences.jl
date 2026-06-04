"""
    validate_probability_vector(p, name = "probabilities") -> Vector{Float64}

Convert `p` to `Vector{Float64}` and check that it is finite, non-negative,
non-empty, and sums to one.
"""
function validate_probability_vector(p::AbstractVector{<:Real},
                                     name::AbstractString = "probabilities")
    isempty(p) && throw(ArgumentError("$name must be non-empty"))
    q = Float64.(p)
    all(isfinite, q) || throw(ArgumentError("$name must contain only finite values"))
    all(≥(0), q) || throw(ArgumentError("$name must be non-negative"))
    s = sum(q)
    s > 0 || throw(ArgumentError("$name must have positive total mass"))
    isapprox(s, 1.0; atol = 1e-8) ||
        throw(ArgumentError("$name must sum to 1, got $s"))
    return q
end

"""
    validate_positive_vector(x, name = "values") -> Vector{Float64}

Convert `x` to `Vector{Float64}` and check that it is finite, positive, and
non-empty.
"""
function validate_positive_vector(x::AbstractVector{<:Real},
                                  name::AbstractString = "values")
    isempty(x) && throw(ArgumentError("$name must be non-empty"))
    y = Float64.(x)
    all(isfinite, y) || throw(ArgumentError("$name must contain only finite values"))
    all(>(0), y) || throw(ArgumentError("$name must be positive"))
    return y
end

"""
    validate_alphabet(alphabet) -> alphabet

Check that `alphabet` is non-empty and contains no duplicate entries.
"""
function validate_alphabet(alphabet)
    length(alphabet) ≥ 1 || throw(ArgumentError("alphabet must be non-empty"))
    length(unique(collect(alphabet))) == length(alphabet) ||
        throw(ArgumentError("alphabet entries must be unique"))
    return alphabet
end

"""
    bin_counts(marginal, n) -> Vector{Int}

Return integer bin counts for `n` observations with proportions as close as
possible to `marginal`.

The counts are obtained by flooring `n .* marginal` and distributing the
remaining observations to the largest fractional remainders. Ties are broken by
alphabet order, which makes the result deterministic.
"""
function bin_counts(marginal::AbstractVector{<:Real}, n::Int)
    n ≥ 0 || throw(ArgumentError("n must be non-negative, got $n"))
    p = validate_probability_vector(marginal, "marginal")
    raw = n .* p
    counts = floor.(Int, raw)
    remaining = n - sum(counts)
    if remaining > 0
        order = sortperm(raw .- counts; rev = true, alg = Base.Sort.MergeSort)
        @inbounds for idx in @view(order[1:remaining])
            counts[idx] += 1
        end
    end
    return counts
end

"""
    quantize_to_symbols(x, alphabet, marginal) -> Vector

Map a real-valued sequence `x` to symbols from `alphabet` using rank binning.

The sorted values of `x` are partitioned into integer bin counts from
[`bin_counts`](@ref), then mapped back to the original order. This avoids
threshold edge cases and makes each finite-sample marginal as close as possible
to the requested `marginal`.

# Arguments
- `x::AbstractVector{<:Real}`: real-valued input sequence.
- `alphabet`: ordered collection of `k` unique symbols.
- `marginal::AbstractVector{<:Real}`: target symbol probabilities.

# Examples
```julia
julia> x = randn(1000)
julia> s = quantize_to_symbols(x, [:L, :M, :H], [0.25, 0.5, 0.25])
julia> count(==(:M), s) / 1000 ≈ 0.5
true
```
"""
function quantize_to_symbols(x::AbstractVector{<:Real}, alphabet,
                             marginal::AbstractVector{<:Real})
    k = length(alphabet)
    validate_alphabet(alphabet)
    length(marginal) == k ||
        throw(ArgumentError(
            "marginal length $(length(marginal)) ≠ alphabet length $k"))
    counts = bin_counts(marginal, length(x))
    result = Vector{eltype(alphabet)}(undef, length(x))
    order = sortperm(x, alg = Base.Sort.MergeSort)

    pos = 1
    @inbounds for (idx, c) in enumerate(counts)
        sym = alphabet[idx]
        for j in pos:(pos + c - 1)
            result[order[j]] = sym
        end
        pos += c
    end
    return result
end

"""
    weighted_sample(rng, weights) -> Int

Draw an index from `1:length(weights)` with probability proportional to `weights`.
Uses the sequential inverse-CDF method; O(k) per call.
"""
function weighted_sample(rng::AbstractRNG, weights::AbstractVector{<:Real})
    u    = rand(rng) * sum(weights)
    cumw = 0.0
    for (i, w) in enumerate(weights)
        cumw += w
        u ≤ cumw && return i
    end
    return length(weights)   # numerical safety for floating-point edge cases
end

"""
    target_marginal(g) -> Vector{Float64}

Return the marginal probabilities a generator claims to target.
"""
target_marginal(g::LRDGenerator) =
    throw(MethodError(target_marginal, (g,)))

"""
    empirical_marginal(seq, alphabet) -> Vector{Float64}

Estimate the marginal distribution of `seq` over `alphabet`.
"""
function empirical_marginal(seq::AbstractVector, alphabet)
    validate_alphabet(alphabet)
    counts = zeros(Int, length(alphabet))
    index = Dict(s => i for (i, s) in enumerate(alphabet))
    @inbounds for s in seq
        i = get(index, s, 0)
        i == 0 && throw(ArgumentError("sequence contains symbol outside alphabet: $s"))
        counts[i] += 1
    end
    return counts ./ length(seq)
end

"""
    empirical_bigram(seq, alphabet) -> Matrix{Float64}

Estimate row-normalised bigram transition probabilities over `alphabet`.
Rows with no observations are left as zeros.
"""
function empirical_bigram(seq::AbstractVector, alphabet)
    validate_alphabet(alphabet)
    k = length(alphabet)
    counts = zeros(Int, k, k)
    index = Dict(s => i for (i, s) in enumerate(alphabet))
    for t in 1:(length(seq) - 1)
        i = get(index, seq[t], 0)
        j = get(index, seq[t + 1], 0)
        (i == 0 || j == 0) &&
            throw(ArgumentError("sequence contains symbol outside alphabet"))
        counts[i, j] += 1
    end
    probs = zeros(Float64, k, k)
    for i in 1:k
        rowsum = sum(@view counts[i, :])
        rowsum > 0 && (probs[i, :] .= counts[i, :] ./ rowsum)
    end
    return probs
end

"""
    empirical_trigram(seq, alphabet) -> Array{Float64,3}

Estimate trigram probabilities `P(X[t+2] | X[t], X[t+1])` over `alphabet`.
Slices with no observations are left as zeros.
"""
function empirical_trigram(seq::AbstractVector, alphabet)
    validate_alphabet(alphabet)
    k = length(alphabet)
    counts = zeros(Int, k, k, k)
    index = Dict(s => i for (i, s) in enumerate(alphabet))
    for t in 1:(length(seq) - 2)
        i = get(index, seq[t], 0)
        j = get(index, seq[t + 1], 0)
        l = get(index, seq[t + 2], 0)
        (i == 0 || j == 0 || l == 0) &&
            throw(ArgumentError("sequence contains symbol outside alphabet"))
        counts[i, j, l] += 1
    end
    probs = zeros(Float64, k, k, k)
    for i in 1:k, j in 1:k
        rowsum = sum(@view counts[i, j, :])
        rowsum > 0 && (probs[i, j, :] .= counts[i, j, :] ./ rowsum)
    end
    return probs
end

"""
    total_variation(p, q) -> Float64

Return the total variation distance between two probability arrays.
"""
function total_variation(p::AbstractArray{<:Real}, q::AbstractArray{<:Real})
    size(p) == size(q) ||
        throw(ArgumentError("probability arrays must have matching sizes"))
    return 0.5 * sum(abs.(Float64.(p) .- Float64.(q)))
end

"""
    rowwise_total_variation(observed, target) -> Vector{Float64}

Return total variation distance for each row of two transition matrices.
"""
function rowwise_total_variation(observed::AbstractMatrix{<:Real},
                                 target::AbstractMatrix{<:Real})
    size(observed) == size(target) ||
        throw(ArgumentError("transition matrices must have matching sizes"))
    return [total_variation(@view(observed[i, :]), @view(target[i, :]))
            for i in axes(observed, 1)]
end
