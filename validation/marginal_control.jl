using S5
using StableRNGs
using Statistics

struct MarginalScenario
    name::String
    marginal::Vector{Float64}
end

function uniform_scenario(k::Int)
    MarginalScenario("uniform", fill(1 / k, k))
end

function moderate_skew_scenario(k::Int)
    p = collect(range(1.0, 2.0; length = k))
    p ./= sum(p)
    MarginalScenario("moderate_skew", p)
end

function zipf_scenario(k::Int; exponent::Float64 = 1.1)
    p = [i^(-exponent) for i in 1:k]
    p ./= sum(p)
    MarginalScenario("zipf_$(exponent)", p)
end

function marginal_errors(generator_factory, n::Int, p::Vector{Float64},
                         alphabet; replicates::Int, seed::Int)
    tv = Vector{Float64}(undef, replicates)
    maxabs = Vector{Float64}(undef, replicates)

    for r in 1:replicates
        rng = StableRNG(seed + r)
        g = generator_factory(p, alphabet)
        seq = generate(g, n; rng)
        observed = empirical_marginal(seq, alphabet)
        tv[r] = total_variation(observed, target_marginal(g))
        maxabs[r] = maximum(abs.(observed .- target_marginal(g)))
    end

    return (;
        tv_mean = mean(tv),
        tv_max = maximum(tv),
        maxabs_mean = mean(maxabs),
        maxabs_max = maximum(maxabs),
    )
end

function run_marginal_control(; ns = (1_000, 5_000),
                                ks = (2, 8, 32),
                                replicates = 20,
                                seed = 20260604)
    rows = NamedTuple[]

    factories = [
        ("SpectralFGN", (p, alphabet) -> SpectralFGN(0.8, alphabet, p)),
        ("LGCM", (p, alphabet) -> LGCM(0.8, alphabet, p; calibration_iters = 8)),
        ("WaveletMarkov", (p, alphabet) -> _iid_wavelet_markov(p, alphabet)),
        ("LAMP", (p, alphabet) -> LAMP(0.5, alphabet, p; d = 200, epsilon = 0.05)),
        ("DyadicLAMP", (p, alphabet) -> DyadicLAMP(0.5, alphabet, p;
                                                   d = 10_000, epsilon = 0.05)),
        ("OnOffMarkov", (p, alphabet) -> _iid_onoff_markov(p, alphabet)),
        ("FSS", (p, alphabet) -> FSS(1.5, alphabet; rates = p)),
        ("HawkesSymbol", (p, alphabet) -> _baseline_hawkes_symbol(p, alphabet)),
    ]

    for k in ks
        alphabet = Symbol.("s" .* string.(1:k))
        scenarios = (uniform_scenario(k), moderate_skew_scenario(k), zipf_scenario(k))

        for scenario in scenarios, n in ns, (method, factory) in factories
            errors = marginal_errors(factory, n, scenario.marginal, alphabet;
                                     replicates, seed)
            push!(rows, (;
                method,
                scenario = scenario.name,
                n,
                k,
                replicates,
                errors...,
            ))
        end
    end

    return rows
end

function _iid_wavelet_markov(p, alphabet)
    k = length(alphabet)
    P = repeat(reshape(p, 1, k), k, 1)
    WaveletMarkov(0.8, alphabet, [P, P]; regime_weights = [0.5, 0.5])
end

function _iid_onoff_markov(p, alphabet)
    k = length(alphabet)
    P = repeat(reshape(p, 1, k), k, 1)
    Q = [0.2 0.8; 0.8 0.2]
    OnOffMarkov(1.5, alphabet, [P, P], Q)
end

function _baseline_hawkes_symbol(p, alphabet)
    k = length(alphabet)
    HawkesSymbol(0.6, alphabet; baseline = p, excitation = zeros(k, k), d = 200)
end

function print_results(rows)
    println("method,scenario,n,k,replicates,tv_mean,tv_max,maxabs_mean,maxabs_max")
    for row in rows
        println(join((
            row.method,
            row.scenario,
            row.n,
            row.k,
            row.replicates,
            row.tv_mean,
            row.tv_max,
            row.maxabs_mean,
            row.maxabs_max,
        ), ","))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    rows = run_marginal_control()
    print_results(rows)
end
