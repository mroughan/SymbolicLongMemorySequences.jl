using SymbolicLongMemorySequences
using StableRNGs
using Statistics

struct LocalStructureScenario
    name::String
    transition_matrix::Matrix{Float64}
end

function iid_scenario(k::Int)
    k ≥ 1 || throw(ArgumentError("iid scenario requires k ≥ 1"))
    return LocalStructureScenario("iid", fill(1.0 / k, k, k))
end

function persistent_scenario(k::Int; persistence::Float64 = 0.8)
    k ≥ 1 || throw(ArgumentError("persistent scenario requires k ≥ 1"))
    0.0 ≤ persistence ≤ 1.0 ||
        throw(ArgumentError("persistence must be in [0, 1]"))
    off_diagonal = k == 1 ? 0.0 : (1 - persistence) / (k - 1)
    P = fill(off_diagonal, k, k)
    for i in 1:k
        P[i, i] = persistence
    end
    return LocalStructureScenario("persistent", P)
end

function cyclic_scenario(k::Int; forward_probability::Float64 = 0.8)
    k ≥ 2 || throw(ArgumentError("cyclic scenario requires k ≥ 2"))
    0.0 ≤ forward_probability ≤ 1.0 ||
        throw(ArgumentError("forward_probability must be in [0, 1]"))
    P = fill((1 - forward_probability) / (k - 1), k, k)
    for i in 1:k
        P[i, mod1(i + 1, k)] = forward_probability
    end
    return LocalStructureScenario("cyclic", P)
end

function local_structure_errors(generator_factory, n::Int,
                                scenario::LocalStructureScenario, alphabet;
                                replicates::Int, seed::Int)
    tv = Vector{Float64}(undef, replicates)
    row_mean = Vector{Float64}(undef, replicates)
    row_max = Vector{Float64}(undef, replicates)

    for r in 1:replicates
        rng = StableRNG(seed + r)
        g = generator_factory(scenario.transition_matrix, alphabet)
        observed = empirical_bigram(generate(g, n; rng), alphabet)
        row_errors = rowwise_total_variation(observed, scenario.transition_matrix)
        tv[r] = total_variation(observed, scenario.transition_matrix) / length(alphabet)
        row_mean[r] = mean(row_errors)
        row_max[r] = maximum(row_errors)
    end

    return (;
        tv_mean = mean(tv),
        tv_max = maximum(tv),
        row_tv_mean = mean(row_mean),
        row_tv_max = maximum(row_max),
    )
end

function run_local_structure_control(; ns = (2_000, 10_000),
                                       ks = (2, 4),
                                       replicates = 20,
                                       seed = 20260611)
    rows = NamedTuple[]
    factories = (
        ("WaveletMarkov", _identical_wavelet_markov),
        ("OnOffMarkov", _identical_onoff_markov),
    )

    for k in ks
        alphabet = Symbol.("s" .* string.(1:k))
        scenarios = (iid_scenario(k), persistent_scenario(k), cyclic_scenario(k))
        for scenario in scenarios, n in ns, (method, factory) in factories
            errors = local_structure_errors(factory, n, scenario, alphabet;
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

function _identical_wavelet_markov(P, alphabet)
    spec = MarkovSpec(alphabet, P)
    return WaveletMarkov(0.8, [spec, spec])
end

function _identical_onoff_markov(P, alphabet)
    spec = MarkovSpec(alphabet, P)
    Q = [0.2 0.8; 0.8 0.2]
    return OnOffMarkov(1.5, [spec, spec], Q)
end

function print_local_structure_results(rows)
    println("method,scenario,n,k,replicates,tv_mean,tv_max,row_tv_mean,row_tv_max")
    for row in rows
        println(join((
            row.method,
            row.scenario,
            row.n,
            row.k,
            row.replicates,
            row.tv_mean,
            row.tv_max,
            row.row_tv_mean,
            row.row_tv_max,
        ), ","))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    print_local_structure_results(run_local_structure_control())
end
