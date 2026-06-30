using FFTW
using Statistics

"""
    indicator_series(seq, symbol; center = true) -> Vector{Float64}

Convert a symbolic sequence into a numeric one-hot indicator series for `symbol`.
Values are `1.0` where `seq[t] == symbol` and `0.0` elsewhere. With
`center = true`, subtract the sample mean so long-memory diagnostics are not
dominated by the marginal symbol frequency.
"""
function indicator_series(seq::AbstractVector, symbol; center::Bool = true)
    x = Float64.(seq .== symbol)
    center && (x .-= mean(x))
    return x
end

"""
    centered_indicator_series(seq, alphabet) -> Vector{Vector{Float64}}

Return one centered one-hot indicator series per symbol in `alphabet`.

Constant series are retained by this function so callers can decide whether to
skip undefined autocorrelation values or report the zero-variance condition.
"""
function centered_indicator_series(seq::AbstractVector, alphabet)
    return [indicator_series(seq, symbol; center = true) for symbol in alphabet]
end

function nextpow2int(n::Int)
    p = 1
    while p < n
        p <<= 1
    end
    return p
end

"""
    fft_unbiased_autocorrelation(x, maxlag) -> Vector{Float64}

Return autocorrelations for lags `1:maxlag` using the convention previously used
by the visual SymbolicLongMemorySequences validation plots: FFT convolution, sample mean removed by the
caller, lag covariance divided by `n - lag`, and normalization by lag-zero
covariance divided by `n`.
"""
function fft_unbiased_autocorrelation(x::AbstractVector{<:Real}, maxlag::Int)
    n = length(x)
    0 ≤ maxlag < n || throw(ArgumentError("maxlag must satisfy 0 ≤ maxlag < length(x)"))
    c0 = mean(abs2, x)
    c0 > 0 || throw(ArgumentError("autocorrelation is undefined for zero variance"))

    pad = nextpow2int(2n)
    padded = zeros(Float64, pad)
    padded[1:n] .= x
    F = fft(padded)
    raw = real(ifft(abs2.(F)))

    acf = Vector{Float64}(undef, maxlag)
    @inbounds for lag in 1:maxlag
        acf[lag] = (raw[lag + 1] / (n - lag)) / c0
    end
    return acf
end

"""
    fft_periodogram_cycles(x) -> frequencies, power

Return the positive-frequency periodogram used by the SymbolicLongMemorySequences validation plots.
Frequencies are cycles per observation and exclude the zero-frequency term.
"""
function fft_periodogram_cycles(x::AbstractVector{<:Real})
    n = length(x)
    G = fft(Float64.(x))
    power = [abs2(G[i + 1]) / n for i in 1:(n ÷ 2)]
    freqs = collect(1:(n ÷ 2)) ./ n
    return freqs, power
end

"""
    longmemory_compatible_autocovariance(x, k) -> Vector{Float64}

Return lag `0:k-1` autocovariances using the convention in LongMemory.jl:
subtract the sample mean and divide every lag sum by `T = length(x)`.
"""
function longmemory_compatible_autocovariance(x::AbstractVector{<:Real}, k::Int)
    T = length(x)
    1 ≤ k < T || throw(ArgumentError("k must satisfy 1 ≤ k < length(x)"))
    y = Float64.(x)
    μ = mean(y)
    acv = Vector{Float64}(undef, k)
    @inbounds for lag in 0:(k - 1)
        total = 0.0
        for t in 1:(T - lag)
            total += (y[t] - μ) * (y[t + lag] - μ)
        end
        acv[lag + 1] = total / T
    end
    return acv
end

"""
    longmemory_compatible_autocorrelation(x, k) -> Vector{Float64}

Return lag `0:k-1` autocorrelations using the LongMemory.jl autocovariance
normalization. This is useful for testing SymbolicLongMemorySequences transformations without making
LongMemory.jl a runtime dependency of SymbolicLongMemorySequences.jl.
"""
function longmemory_compatible_autocorrelation(x::AbstractVector{<:Real}, k::Int)
    acv = longmemory_compatible_autocovariance(x, k)
    acv[1] > 0 || throw(ArgumentError("autocorrelation is undefined for zero variance"))
    return acv ./ acv[1]
end

"""
    longmemory_compatible_periodogram(x) -> frequencies, power

Return the same positive-frequency periodogram scale used by LongMemory.jl,
adapted to SymbolicLongMemorySequences's plotting convention. LongMemory.jl reports angular frequencies
including zero; this helper reports cycles per observation and drops zero.
"""
function longmemory_compatible_periodogram(x::AbstractVector{<:Real})
    n = length(x)
    Iω = abs.(rfft(Float64.(x))) .^ 2 ./ n
    ω = collect(2π .* (0:(n - 1)) ./ n)
    lastidx = iseven(n) ? n ÷ 2 + 1 : cld(n, 2)
    freqs = ω[2:lastidx] ./ (2π)
    power = Iω[2:lastidx]
    return freqs, power
end

function _valid_indicator_series(seq, alphabet)
    series = Vector{Float64}[]
    for x in centered_indicator_series(seq, alphabet)
        mean(abs2, x) > 0 && push!(series, x)
    end
    return series
end

"""
    indicator_diagnostics(seq, alphabet; maxlag = length(seq) ÷ 2)

Compute the visual SymbolicLongMemorySequences diagnostics: signed one-hot autocorrelation averaged across
symbols, and the positive-frequency one-hot periodogram averaged across symbols.

The transformation is explicit: each symbol becomes a centered one-hot series,
zero-variance series are skipped, and diagnostics are averaged across the
remaining symbol indicators.
"""
function indicator_diagnostics(seq, alphabet; maxlag::Int = length(seq) ÷ 2)
    series = _valid_indicator_series(seq, alphabet)
    isempty(series) && throw(ArgumentError("all indicator series have zero variance"))

    acf = zeros(Float64, maxlag)
    freqs = Float64[]
    power = Float64[]

    for x in series
        acf .+= fft_unbiased_autocorrelation(x, maxlag)
        f, pxx = fft_periodogram_cycles(x)
        if isempty(power)
            freqs = f
            power = zeros(Float64, length(pxx))
        end
        power .+= pxx
    end

    acf ./= length(series)
    power ./= length(series)
    return acf, freqs, power
end

"""
    longmemory_indicator_diagnostics(seq, alphabet; maxlag = length(seq) ÷ 2)

Compute one-hot diagnostics using LongMemory.jl-compatible conventions. The
returned autocorrelation covers lags `1:maxlag` after dropping LongMemory's lag
zero value; the periodogram uses cycles per observation after dropping zero
frequency.
"""
function longmemory_indicator_diagnostics(seq, alphabet; maxlag::Int = length(seq) ÷ 2)
    series = _valid_indicator_series(seq, alphabet)
    isempty(series) && throw(ArgumentError("all indicator series have zero variance"))

    acf = zeros(Float64, maxlag)
    freqs = Float64[]
    power = Float64[]

    for x in series
        acf .+= longmemory_compatible_autocorrelation(x, maxlag + 1)[2:end]
        f, pxx = longmemory_compatible_periodogram(x)
        if isempty(power)
            freqs = f
            power = zeros(Float64, length(pxx))
        end
        power .+= pxx
    end

    acf ./= length(series)
    power ./= length(series)
    return acf, freqs, power
end
