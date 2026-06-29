module S5Benchmarks

using BenchmarkTools
using Dates
using StableRNGs

const ROOT = normpath(joinpath(@__DIR__, ".."))
ROOT in LOAD_PATH || pushfirst!(LOAD_PATH, ROOT)

using S5

const DEFAULT_NS = (10_000, 100_000)
const LARGE_NS = (10_000, 100_000, 1_000_000)
const SCALING_NS = (100, 1_000, 10_000, 100_000, 1_000_000)
const DEFAULT_KS = (2, 8, 64)
const SCALING_KS = (2, 8)
const RESULTS_DIR = joinpath(@__DIR__, "results")

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
        "PB4_IntermittentMapSymbols_z=1.6" =>
            IntermittentMapSymbols(1.6, alphabet, p; burnin = 1000),
        "MB1a_LAMP_d=$lamp_d" => LAMP(0.5, alphabet, p; d = lamp_d, epsilon = 0.02),
        "MB1b_DyadicLAMP_d=$dyadic_d" =>
            DyadicLAMP(0.5, alphabet, p; d = dyadic_d, epsilon = 0.02),
        "MB1c_CalibratedAdditiveMarkov_d=$lamp_d" =>
            CalibratedAdditiveMarkov(0.5, alphabet, p; d = lamp_d, strength = 0.75),
        "MB2_OnOffMarkov_regimes=2_Lmin=10" =>
            OnOffMarkov(1.5, alphabet, [P_iid, P_persistent], switch; L_min = 10.0),
        "MB3_FSS_streams=$k" => FSS(1.5, alphabet; rates = p),
        "MB4_HawkesSymbol_d=$lamp_d" =>
            HawkesSymbol(0.6, alphabet; baseline = p, excitation = P_persistent,
                         d = lamp_d),
        "MB5_DuplicationMutation_alpha=1.5" =>
            DuplicationMutation(1.5, alphabet, p; mutation_probability = 0.02,
                                max_block_length = 10_000),
    )
end

"""
    make_suite(; ns = DEFAULT_NS, ks = DEFAULT_KS, samples = 5, seconds = 2)

Create a `BenchmarkGroup` for all implemented generators across sequence lengths
and alphabet sizes. Each sample uses a fresh `StableRNG`, so benchmark trials do
not mutate a shared RNG state.
"""
function make_suite(; ns = DEFAULT_NS, ks = DEFAULT_KS,
                      samples::Int = 5, seconds::Real = 2,
                      syntheses_per_trial::Int = 1)
    syntheses_per_trial ≥ 1 ||
        throw(ArgumentError("syntheses_per_trial must be positive"))
    suite = BenchmarkGroup()
    for k in ks
        suite["k=$k"] = BenchmarkGroup()
        for (method, generator) in _cases(k)
            suite["k=$k"][method] = BenchmarkGroup()
            for n in ns
                seeds = [20260611 + 1_000_000 * k + n + i
                         for i in 1:syntheses_per_trial]
                suite["k=$k"][method]["n=$n"] =
                    @benchmarkable _generate_repeated($generator, $n, $seeds) samples = samples seconds = seconds evals = 1
            end
        end
    end
    return suite
end

function _generate_repeated(generator, n::Int, seeds::AbstractVector{<:Integer})
    for seed in seeds
        generate(generator, n; rng = StableRNG(seed))
    end
    return nothing
end

function _parse_int_suffix(label::AbstractString, prefix::AbstractString)
    startswith(label, prefix) ||
        throw(ArgumentError("label $label does not start with $prefix"))
    return parse(Int, label[(lastindex(prefix) + 1):end])
end

function _trial_row(trial, method::AbstractString, k::Int, n::Int,
                    syntheses_per_trial::Int)
    estimate = minimum(trial)
    total_time_ns = Float64(estimate.time)
    return (;
        method = String(method),
        k,
        n,
        syntheses_per_trial,
        total_time_ns,
        total_time_ms = total_time_ns / 1.0e6,
        time_ns = total_time_ns / syntheses_per_trial,
        time_ms = total_time_ns / syntheses_per_trial / 1.0e6,
        memory_bytes = Int(estimate.memory),
        allocations = Int(estimate.allocs),
    )
end

function collect_rows(results; syntheses_per_trial::Int = 1)
    rows = NamedTuple[]
    for klabel in sort(collect(keys(results)); by = x -> _parse_int_suffix(x, "k="))
        k = _parse_int_suffix(klabel, "k=")
        group = results[klabel]
        for method in sort(collect(keys(group)))
            method_group = group[method]
            for nlabel in sort(collect(keys(method_group)); by = x -> _parse_int_suffix(x, "n="))
                n = _parse_int_suffix(nlabel, "n=")
                push!(rows, _trial_row(method_group[nlabel], method, k, n,
                                       syntheses_per_trial))
            end
        end
    end
    return rows
end

function _csv_cell(x)
    s = string(x)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

function write_csv(path::AbstractString, rows)
    open(path, "w") do io
        println(io, "method,k,n,syntheses_per_trial,total_time_ns,total_time_ms,time_ns,time_ms,memory_bytes,allocations")
        for row in rows
            println(io, join((
                _csv_cell(row.method),
                row.k,
                row.n,
                row.syntheses_per_trial,
                row.total_time_ns,
                row.total_time_ms,
                row.time_ns,
                row.time_ms,
                row.memory_bytes,
                row.allocations,
            ), ","))
        end
    end
    return path
end

function svg_escape(s)
    replace(replace(replace(string(s), "&" => "&amp;"), "<" => "&lt;"), ">" => "&gt;")
end

function _method_short(label::AbstractString)
    return first(split(label, "_"))
end

function _sanitize(label::AbstractString)
    replace(label, r"[^A-Za-z0-9]+" => "_") |> x -> strip(x, '_')
end

function _bar_color(i::Int)
    palette = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
               "#4c78a8")
    return palette[mod1(i, length(palette))]
end

function write_relative_time_svg(path::AbstractString, rows, k::Int, n::Int)
    selected = [row for row in rows if row.k == k && row.n == n]
    isempty(selected) && throw(ArgumentError("no benchmark rows for k=$k, n=$n"))
    sort!(selected; by = row -> row.time_ms)
    fastest = first(selected).time_ms
    ratios = [row.time_ms / fastest for row in selected]
    maxratio = maximum(ratios)

    width = 1200
    row_h = 34
    top = 70
    left = 300
    right = 130
    height = top + row_h * length(selected) + 70
    plot_w = width - left - right

    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<text x="$(width / 2)" y="30" text-anchor="middle" font-family="sans-serif" font-size="22">Relative generation time, k=$k, n=$n</text>""")
        println(io, """<text x="$(width / 2)" y="52" text-anchor="middle" font-family="sans-serif" font-size="13">Bars are normalized to the fastest method in this benchmark slice</text>""")
        for (i, row) in enumerate(selected)
            y = top + (i - 1) * row_h
            ratio = ratios[i]
            bar_w = plot_w * ratio / maxratio
            println(io, """<text x="$(left - 12)" y="$(y + 21)" text-anchor="end" font-family="sans-serif" font-size="12">$(svg_escape(_method_short(row.method)))</text>""")
            println(io, """<rect x="$left" y="$(y + 6)" width="$(round(bar_w; digits = 2))" height="20" fill="$(_bar_color(i))"/>""")
            println(io, """<text x="$(left + bar_w + 8)" y="$(y + 21)" font-family="sans-serif" font-size="12">$(round(ratio; digits = 2))x ($(round(row.time_ms; digits = 3)) ms)</text>""")
        end
        println(io, """<text x="$(left + plot_w / 2)" y="$(height - 22)" text-anchor="middle" font-family="sans-serif" font-size="14">relative time, fastest = 1x</text>""")
        println(io, """</svg>""")
    end
    return path
end

function write_scaling_svg(path::AbstractString, rows, k::Int)
    selected = [row for row in rows if row.k == k]
    isempty(selected) && throw(ArgumentError("no benchmark rows for k=$k"))
    methods = sort(unique(row.method for row in selected))
    ns = sort(unique(row.n for row in selected))
    xmin, xmax = extrema(Float64.(ns))
    ymin, ymax = extrema(row.time_ms for row in selected)
    lxmin, lxmax = log10(xmin), log10(xmax)
    lymin, lymax = log10(ymin), log10(ymax)
    lymin == lymax && (lymax += 1)

    width, height = 1200, 760
    ml, mr, mt, mb = 90, 260, 70, 85
    plot_left, plot_right = ml, width - mr
    plot_top, plot_bottom = mt, height - mb
    plot_w = plot_right - plot_left
    plot_h = plot_bottom - plot_top
    sx(v) = plot_left + (log10(v) - lxmin) / (lxmax - lxmin) * plot_w
    sy(v) = plot_bottom - (log10(v) - lymin) / (lymax - lymin) * plot_h

    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<text x="$(width / 2)" y="30" text-anchor="middle" font-family="sans-serif" font-size="22">Generation-time scaling, k=$k</text>""")
        println(io, """<text x="$(width / 2)" y="52" text-anchor="middle" font-family="sans-serif" font-size="13">Log-log plot of BenchmarkTools minimum time against sequence length</text>""")
        println(io, """<line x1="$plot_left" y1="$plot_bottom" x2="$plot_right" y2="$plot_bottom" stroke="black"/>""")
        println(io, """<line x1="$plot_left" y1="$plot_top" x2="$plot_left" y2="$plot_bottom" stroke="black"/>""")
        for n in ns
            x = sx(n)
            println(io, """<line x1="$x" y1="$plot_bottom" x2="$x" y2="$(plot_bottom + 5)" stroke="black"/>""")
            println(io, """<text x="$x" y="$(plot_bottom + 24)" text-anchor="middle" font-family="sans-serif" font-size="11">$(n)</text>""")
        end
        for e in floor(Int, lymin):ceil(Int, lymax)
            val = 10.0^e
            ymin <= val <= ymax || continue
            y = sy(val)
            println(io, """<line x1="$(plot_left - 5)" y1="$y" x2="$plot_left" y2="$y" stroke="black"/>""")
            println(io, """<text x="$(plot_left - 10)" y="$(y + 4)" text-anchor="end" font-family="sans-serif" font-size="11">10^$e</text>""")
        end
        for (i, method) in enumerate(methods)
            pts = sort([row for row in selected if row.method == method]; by = row -> row.n)
            poly = join(["$(round(sx(row.n); digits = 2)),$(round(sy(row.time_ms); digits = 2))" for row in pts], " ")
            color = _bar_color(i)
            println(io, """<polyline fill="none" stroke="$color" stroke-width="2" points="$poly"/>""")
            for row in pts
                println(io, """<circle cx="$(round(sx(row.n); digits = 2))" cy="$(round(sy(row.time_ms); digits = 2))" r="3" fill="$color"/>""")
            end
            ly = mt + 20 + (i - 1) * 22
            println(io, """<line x1="$(plot_right + 28)" y1="$ly" x2="$(plot_right + 58)" y2="$ly" stroke="$color" stroke-width="2"/>""")
            println(io, """<text x="$(plot_right + 66)" y="$(ly + 4)" font-family="sans-serif" font-size="12">$(svg_escape(_method_short(method)))</text>""")
        end
        println(io, """<text x="$(plot_left + plot_w / 2)" y="$(height - 28)" text-anchor="middle" font-family="sans-serif" font-size="15">sequence length n</text>""")
        println(io, """<text x="28" y="$(plot_top + plot_h / 2)" transform="rotate(-90 28 $(plot_top + plot_h / 2))" text-anchor="middle" font-family="sans-serif" font-size="15">time (ms)</text>""")
        println(io, """</svg>""")
    end
    return path
end

function write_results_markdown(path::AbstractString, rows; mode::Symbol, samples::Int,
                                seconds::Real, syntheses_per_trial::Int,
                                created::DateTime = now())
    ks = sort(unique(row.k for row in rows))
    ns = sort(unique(row.n for row in rows))
    maxn = maximum(ns)
    maxk = maximum(ks)
    slice = sort([row for row in rows if row.k == maxk && row.n == maxn];
                 by = row -> row.time_ms)
    fastest = first(slice).time_ms

    open(path, "w") do io
        println(io, "# Benchmark Results")
        println(io)
        timestamp = Dates.format(created, dateformat"yyyy-mm-dd HH:MM")
        println(io, "Generated: $timestamp")
        println(io)
        println(io, "These results are machine-specific BenchmarkTools measurements of generator hot paths after construction. Times are the minimum observed generation time for each benchmarkable, so they are useful for relative comparison on this machine rather than as portable guarantees.")
        println(io)
        println(io, "- Suite: `$mode`")
        println(io, "- Sequence lengths: `$(join(ns, "`, `"))`")
        println(io, "- Alphabet sizes: `$(join(ks, "`, `"))`")
        println(io, "- Samples per benchmark: `$samples`")
        println(io, "- Seconds budget per benchmark: `$seconds`")
        println(io, "- Syntheses per BenchmarkTools trial: `$syntheses_per_trial`")
        println(io, "- Reported times: per-synthesis averages from each trial")
        println(io)
        println(io, "## Relative Time")
        println(io)
        println(io, "The relative-time plots are histogram-style horizontal bar charts. They normalize each method to the fastest method for the same `k` and largest `n` in this run.")
        println(io)
        for k in ks
            println(io, "![Relative generation time, k=$k](results/relative_times_k$(k)_n$(maxn).svg)")
            println(io)
        end
        println(io, "At `k = $maxk` and `n = $maxn`:")
        println(io)
        println(io, "| Method | Time (ms) | Relative | Trial allocations | Trial memory (bytes) |")
        println(io, "|---|---:|---:|---:|---:|")
        for row in slice
            println(io, "| `$(row.method)` | $(round(row.time_ms; digits = 3)) | $(round(row.time_ms / fastest; digits = 2))x | $(row.allocations) | $(row.memory_bytes) |")
        end
        println(io)
        println(io, "## Scaling With Sequence Length")
        println(io)
        println(io, "The scaling plots show time versus generated sequence length on log-log axes for each fixed alphabet size.")
        println(io)
        for k in ks
            println(io, "![Generation-time scaling, k=$k](results/scaling_k$(k).svg)")
            println(io)
        end
        println(io, "## Interpretation")
        println(io)
        if mode == :scaling
            println(io, "- This rare scaling run omits `k = 64` to focus on sequence-length behavior for `k = 2` and `k = 8`.")
        end
        println(io, "- `OnOffMarkov`, `FSS`, and `DuplicationMutation` are the fastest cases in this grid because their hot paths avoid scanning long histories.")
        println(io, "- `SpectralFGN`, `IntermittentMapSymbols`, and spectral-driver `WaveletMarkov` scale well with `n`; their visible cost is mostly FFT/rank-binning and sequential emission work.")
        println(io, "- `LGCM` grows strongly with alphabet size because it generates one latent stream per symbol and performs calibration/argmax work.")
        println(io, "- `LAMP`, `CalibratedAdditiveMarkov`, `HawkesSymbol`, and especially `DyadicLAMP` pay for explicit or approximate history handling. Their relative cost rises when configured memory depth or alphabet size increases.")
    end
    return path
end

function write_artifacts(rows; mode::Symbol, samples::Int, seconds::Real,
                         syntheses_per_trial::Int,
                         outdir::AbstractString = RESULTS_DIR)
    mkpath(outdir)
    for path in readdir(outdir; join = true)
        if endswith(path, ".csv") || endswith(path, ".svg")
            rm(path; force = true)
        end
    end
    write_csv(joinpath(outdir, "benchmarks.csv"), rows)
    ks = sort(unique(row.k for row in rows))
    maxn = maximum(row.n for row in rows)
    for k in ks
        write_relative_time_svg(joinpath(outdir, "relative_times_k$(k)_n$(maxn).svg"),
                                rows, k, maxn)
        write_scaling_svg(joinpath(outdir, "scaling_k$(k).svg"), rows, k)
    end
    write_results_markdown(joinpath(@__DIR__, "RESULTS.md"), rows;
                           mode, samples, seconds, syntheses_per_trial)
end

function run_suite(; large::Bool = get(ENV, "S5_BENCHMARK_LARGE", "false") == "true",
                     scaling::Bool = get(ENV, "S5_BENCHMARK_SCALING", "false") == "true",
                     samples::Int = parse(Int, get(ENV, "S5_BENCHMARK_SAMPLES",
                                                   scaling ? "3" : "5")),
                     seconds::Real = parse(Float64, get(ENV, "S5_BENCHMARK_SECONDS",
                                                        scaling ? "0.1" : "2")),
                     syntheses_per_trial::Int =
                         parse(Int, get(ENV, "S5_BENCHMARK_SYNTH_REPEATS",
                                        scaling ? "10" : "1")),
                     write_results::Bool = get(ENV, "S5_BENCHMARK_WRITE_RESULTS", "true") == "true")
    mode = scaling ? :scaling : (large ? :large : :default)
    ns = scaling ? SCALING_NS : (large ? LARGE_NS : DEFAULT_NS)
    ks = scaling ? SCALING_KS : DEFAULT_KS
    suite = make_suite(; ns, ks, samples, seconds, syntheses_per_trial)
    results = run(suite; verbose = true)
    display(results)
    if write_results
        rows = collect_rows(results; syntheses_per_trial)
        write_artifacts(rows; mode, samples, seconds, syntheses_per_trial)
        println("Wrote benchmark results to $(joinpath(@__DIR__, "RESULTS.md"))")
    end
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_suite()
end

end # module S5Benchmarks
