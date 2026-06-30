const ROOT = normpath(joinpath(@__DIR__, ".."))
ROOT in LOAD_PATH || pushfirst!(LOAD_PATH, ROOT)

using LongMemory
using SymbolicLongMemorySequences
using StableRNGs
using Statistics

include(joinpath(@__DIR__, "lrd_symbol_diagnostics.jl"))

function longmemory_package_periodogram_cycles(x::AbstractVector{<:Real})
    power, angular_freqs = LongMemory.periodogram(Float64.(x))
    return angular_freqs[2:end] ./ (2π), power[2:end]
end

function compare_longmemory_series(x::AbstractVector{<:Real}; k::Int)
    local_acv = longmemory_compatible_autocovariance(x, k)
    local_acf = longmemory_compatible_autocorrelation(x, k)
    local_freqs, local_power = longmemory_compatible_periodogram(x)

    package_acv = vec(LongMemory.autocovariance(Float64.(x), k))
    package_acf = vec(LongMemory.autocorrelation(Float64.(x), k))
    package_freqs, package_power = longmemory_package_periodogram_cycles(x)

    return (;
        autocovariance_maxdiff = maximum(abs.(local_acv .- package_acv)),
        autocorrelation_maxdiff = maximum(abs.(local_acf .- package_acf)),
        frequency_maxdiff = maximum(abs.(local_freqs .- package_freqs)),
        periodogram_maxdiff = maximum(abs.(local_power .- package_power)),
    )
end

function compare_longmemory_indicators(seq, alphabet; maxlag::Int = min(30, length(seq) ÷ 4))
    rows = NamedTuple[]
    for symbol in alphabet
        x = indicator_series(seq, symbol; center = true)
        mean(abs2, x) == 0 && continue
        diffs = compare_longmemory_series(x; k = maxlag + 1)
        push!(rows, (; symbol = string(symbol), maxlag, diffs...))
    end
    return rows
end

function print_comparison(rows)
    println(join((
        "symbol",
        "maxlag",
        "autocovariance_maxdiff",
        "autocorrelation_maxdiff",
        "frequency_maxdiff",
        "periodogram_maxdiff",
    ), ","))
    for row in rows
        println(join((
            row.symbol,
            row.maxlag,
            row.autocovariance_maxdiff,
            row.autocorrelation_maxdiff,
            row.frequency_maxdiff,
            row.periodogram_maxdiff,
        ), ","))
    end
end

function run_longmemory_comparison(; n::Int = 2_048,
                                    maxlag::Int = 64,
                                    seed::Int = 20260611)
    alphabet = ['A', 'B', 'C', 'D']
    generator = SpectralFGN(0.8, alphabet, fill(1 / length(alphabet), length(alphabet)))
    seq = generate(generator, n; rng = StableRNG(seed))
    rows = compare_longmemory_indicators(seq, alphabet; maxlag)
    print_comparison(rows)
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_longmemory_comparison()
end
