module S5Benchmarks

using BenchmarkTools
using StableRNGs

const ROOT = normpath(joinpath(@__DIR__, ".."))
ROOT in LOAD_PATH || pushfirst!(LOAD_PATH, ROOT)

using S5

const DEFAULT_NS = (10_000, 100_000)
const LARGE_NS = (10_000, 100_000, 1_000_000)
const DEFAULT_KS = (2, 8, 64)

function _alphabet(k::Int)
    k ≥ 1 || throw(ArgumentError("alphabet size must be positive"))
    return Symbol.("s" .* string.(1:k))
end

function _marginal(k::Int)
    p = [i^(-1.1) for i in 1:k]
    p ./= sum(p)
    return p
end

function _iid_matrix(p::AbstractVector{<:Real})
    k = length(p)
    return repeat(reshape(Float64.(p), 1, k), k, 1)
end

function _persistent_matrix(k::Int; persistence::Float64 = 0.75)
    k ≥ 2 || throw(ArgumentError("persistent matrix requires k ≥ 2"))
    offdiag = (1 - persistence) / (k - 1)
    P = fill(offdiag, k, k)
    for i in 1:k
        P[i, i] = persistence
    end
    return P
end

function _cases(k::Int)
    alphabet = _alphabet(k)
    p = _marginal(k)
    P_iid = _iid_matrix(p)
    P_persistent = _persistent_matrix(k)
    switch = [0.2 0.8; 0.8 0.2]
    lamp_d = min(1_000, max(100, 10 * k))
    dyadic_d = 100_000

    return (
        "PB1_SpectralFGN_fft=n" => SpectralFGN(0.8, alphabet, p),
        "PB2_LGCM_iters=8" => LGCM(0.8, alphabet, p; calibration_iters = 8),
        "PB3_WaveletMarkov_spectral_regimes=2" =>
            WaveletMarkov(0.8, alphabet, [P_iid, P_persistent]; driver = :spectral),
        "MB1a_LAMP_d=$lamp_d" => LAMP(0.5, alphabet, p; d = lamp_d, epsilon = 0.02),
        "MB1b_DyadicLAMP_d=$dyadic_d" =>
            DyadicLAMP(0.5, alphabet, p; d = dyadic_d, epsilon = 0.02),
        "MB2_OnOffMarkov_regimes=2_Lmin=10" =>
            OnOffMarkov(1.5, alphabet, [P_iid, P_persistent], switch; L_min = 10.0),
        "MB3_FSS_streams=$k" => FSS(1.5, alphabet; rates = p),
        "MB4_HawkesSymbol_d=$lamp_d" =>
            HawkesSymbol(0.6, alphabet; baseline = p, excitation = P_persistent,
                         d = lamp_d),
    )
end

"""
    make_suite(; ns = DEFAULT_NS, ks = DEFAULT_KS, samples = 5, seconds = 2)

Create a `BenchmarkGroup` for all implemented generators across sequence lengths
and alphabet sizes. Each sample uses a fresh `StableRNG`, so benchmark trials do
not mutate a shared RNG state.
"""
function make_suite(; ns = DEFAULT_NS, ks = DEFAULT_KS,
                      samples::Int = 5, seconds::Real = 2)
    suite = BenchmarkGroup()
    for k in ks
        suite["k=$k"] = BenchmarkGroup()
        for (method, generator) in _cases(k)
            suite["k=$k"][method] = BenchmarkGroup()
            for n in ns
                seed = 20260611 + n + k
                suite["k=$k"][method]["n=$n"] =
                    @benchmarkable generate($generator, $n; rng = rng) setup = begin
                        rng = StableRNG($seed)
                    end samples = samples seconds = seconds evals = 1
            end
        end
    end
    return suite
end

function run_suite(; large::Bool = get(ENV, "S5_BENCHMARK_LARGE", "false") == "true",
                     samples::Int = parse(Int, get(ENV, "S5_BENCHMARK_SAMPLES", "5")),
                     seconds::Real = parse(Float64, get(ENV, "S5_BENCHMARK_SECONDS", "2")))
    ns = large ? LARGE_NS : DEFAULT_NS
    suite = make_suite(; ns, samples, seconds)
    results = run(suite; verbose = true)
    display(results)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_suite()
end

end # module S5Benchmarks
