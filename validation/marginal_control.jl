using SymbolicLongMemorySequences
using Distributions
using StableRNGs
using Statistics

const DEFAULT_MARGINAL_OUTDIR = joinpath(@__DIR__, "results", "marginal_control")
const MARGINAL_VALIDATION_LAMP_REPEAT_PROBABILITY = 0.4

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

function marginal_generator_factories()
    return [
        ("PB1_SpectralFGN", (p, alphabet) -> SpectralFGN(0.8, alphabet, p)),
        ("PB2_LGCM", (p, alphabet) -> LGCM(0.8, alphabet, p; calibration_iters = 8)),
        ("PB3_WaveletMarkov", (p, alphabet) -> _iid_wavelet_markov(p, alphabet)),
        ("PB4_IntermittentMapSymbols", (p, alphabet) ->
            IntermittentMapSymbols(1.6, alphabet, p; burnin = 100)),
        ("MB1a_LAMP", (p, alphabet) -> LAMP(0.5, alphabet, p; d = 200,
                                            epsilon = 0.05,
                                            transition_matrix =
                                                _marginal_lamp_transition(p))),
        ("MB1b_DyadicLAMP", (p, alphabet) -> DyadicLAMP(0.5, alphabet, p;
                                                         d = 10_000,
                                                         epsilon = 0.05,
                                                         transition_matrix =
                                                             _marginal_lamp_transition(p))),
        ("MB1c_CalibratedAdditiveMarkov", (p, alphabet) ->
            CalibratedAdditiveMarkov(0.5, alphabet, p; d = 200, strength = 0.0)),
        ("MB2_OnOffMarkov", (p, alphabet) -> _iid_onoff_markov(p, alphabet)),
        ("MB3_FSS", (p, alphabet) -> FSS(1.5, alphabet; rates = p)),
        ("MB4_HawkesSymbol", (p, alphabet) -> _baseline_hawkes_symbol(p, alphabet)),
        ("MB5_DuplicationMutation", (p, alphabet) ->
            DuplicationMutation(1.5, alphabet, p; mutation_probability = 1.0,
                                seed_length = 32, max_block_length = 100)),
    ]
end

function run_marginal_control(; ns = (1_000, 5_000),
                                ks = (2, 8, 32),
                                replicates = 20,
                                seed = 20260604)
    rows = NamedTuple[]

    factories = marginal_generator_factories()

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

function trimmed_window(n::Int; trim_fraction::Real = 0.1)
    n > 0 || throw(ArgumentError("n must be positive"))
    0 ≤ trim_fraction < 0.5 ||
        throw(ArgumentError("trim_fraction must be in [0, 0.5)"))
    trim = floor(Int, n * Float64(trim_fraction))
    lo = trim + 1
    hi = n - trim
    lo ≤ hi || throw(ArgumentError("trim_fraction leaves no observations"))
    return lo:hi
end

function marginal_counts(seq, alphabet)
    index = Dict{eltype(alphabet), Int}()
    for (i, symbol) in pairs(alphabet)
        index[symbol] = i
    end
    counts = zeros(Int, length(alphabet))
    for symbol in seq
        counts[index[symbol]] += 1
    end
    return counts
end

function chisq_uniform_test(counts::AbstractVector{<:Integer})
    k = length(counts)
    k ≥ 2 || throw(ArgumentError("chi-squared test requires at least two bins"))
    m = sum(counts)
    m > 0 || throw(ArgumentError("chi-squared test requires observations"))
    expected = m / k
    statistic = sum((Float64(c) - expected)^2 / expected for c in counts)
    pvalue = ccdf(Chisq(k - 1), statistic)
    return (; statistic, df = k - 1, pvalue)
end

"""
    indicator_effective_sample_size(seq, symbol; maxlag=min(1_000, length(seq) ÷ 4))

Estimate the effective sample size for the centered one-hot indicator of
`symbol` in `seq`. The estimate uses an integrated autocorrelation time with
initial-positive truncation, so only positive autocorrelations before the first
non-positive lag contribute.

Example output:

```julia
julia> indicator_effective_sample_size(vcat(fill(:a, 50), fill(:b, 50)), :a; maxlag=10)
(effective_n = ..., integrated_autocorrelation_time = ..., used_lags = ..., maxlag = 10)
```
"""
function indicator_effective_sample_size(seq, symbol; maxlag::Int = min(1_000, length(seq) ÷ 4))
    n = length(seq)
    n > 1 || throw(ArgumentError("effective sample size requires at least two observations"))
    0 ≤ maxlag < n || throw(ArgumentError("maxlag must be in [0, n - 1]"))

    p = count(==(symbol), seq) / n
    x = Vector{Float64}(undef, n)
    @inbounds for i in eachindex(seq)
        x[i] = (seq[i] == symbol ? 1.0 : 0.0) - p
    end

    variance = sum(abs2, x) / n
    if variance ≤ eps(Float64)
        return (; effective_n = 1.0, integrated_autocorrelation_time = Float64(n),
                used_lags = 0, maxlag)
    end

    positive_sum = 0.0
    used_lags = 0
    @inbounds for lag in 1:maxlag
        covariance = 0.0
        for i in 1:(n - lag)
            covariance += x[i] * x[i + lag]
        end
        rho = covariance / ((n - lag) * variance)
        (!isfinite(rho) || rho ≤ 0) && break
        positive_sum += rho
        used_lags = lag
    end

    iact = max(1.0, 1.0 + 2.0 * positive_sum)
    effective_n = n / iact
    return (; effective_n, integrated_autocorrelation_time = iact, used_lags, maxlag)
end

"""
    categorical_effective_sample_size(seq, alphabet; maxlag=min(1_000, length(seq) ÷ 4))

Estimate a conservative effective sample size for a categorical sequence by
taking the minimum symbol-level ESS across centered one-hot indicators. This is
used to calibrate marginal-frequency diagnostics under dependence.

Example output:

```julia
julia> categorical_effective_sample_size([:a, :a, :b, :b], [:a, :b]; maxlag=1)
(effective_n = ..., integrated_autocorrelation_time = ..., maxlag = 1, ...)
```
"""
function categorical_effective_sample_size(seq, alphabet;
                                           maxlag::Int = min(1_000, length(seq) ÷ 4))
    n = length(seq)
    symbol_results =
        [indicator_effective_sample_size(seq, symbol; maxlag) for symbol in alphabet]
    symbol_effective_ns = [result.effective_n for result in symbol_results]
    valid = filter(x -> isfinite(x) && x > 0, symbol_effective_ns)
    effective_n = isempty(valid) ? 1.0 : min(Float64(n), minimum(valid))
    iact = n / effective_n
    return (;
        effective_n,
        integrated_autocorrelation_time = iact,
        maxlag,
        symbol_effective_ns,
        symbol_integrated_autocorrelation_times =
            [result.integrated_autocorrelation_time for result in symbol_results],
        symbol_used_lags = [result.used_lags for result in symbol_results],
    )
end

"""
    ess_corrected_chisq_uniform_test(counts, effective_n)

Return an approximate ESS-adjusted uniform chi-squared diagnostic. The raw
chi-squared statistic is scaled by `min(1, effective_n / sum(counts))`, and the
scaled statistic is compared with the same `Chisq(k - 1)` reference as the raw
diagnostic.

Example output:

```julia
julia> ess_corrected_chisq_uniform_test([60, 40], 50)
(statistic = 2.0, df = 1, pvalue = ..., scale = 0.5)
```
"""
function ess_corrected_chisq_uniform_test(counts::AbstractVector{<:Integer},
                                          effective_n::Real)
    effective_n > 0 || throw(ArgumentError("effective_n must be positive"))
    raw = chisq_uniform_test(counts)
    m = sum(counts)
    scale = min(1.0, Float64(effective_n) / m)
    statistic = raw.statistic * scale
    pvalue = ccdf(Chisq(raw.df), statistic)
    return (; statistic, df = raw.df, pvalue, scale)
end

function uniform_marginal_replicates(generator_factory, n::Int, alphabet;
                                     replicates::Int, seed::Int,
                                     trim_fraction::Real)
    k = length(alphabet)
    p = fill(1 / k, k)
    full_freqs = Matrix{Float64}(undef, replicates, k)
    freqs = Matrix{Float64}(undef, replicates, k)
    stats = Vector{Float64}(undef, replicates)
    pvalues = Vector{Float64}(undef, replicates)
    effective_ns = Vector{Float64}(undef, replicates)
    iacts = Vector{Float64}(undef, replicates)
    ess_stats = Vector{Float64}(undef, replicates)
    ess_pvalues = Vector{Float64}(undef, replicates)
    trimmed_n = 0

    for r in 1:replicates
        rng = StableRNG(seed + r)
        g = generator_factory(p, alphabet)
        seq = generate(g, n; rng)
        full_counts = marginal_counts(seq, alphabet)
        full_freqs[r, :] .= full_counts ./ length(seq)
        window = trimmed_window(length(seq); trim_fraction)
        trimmed = @view seq[window]
        counts = marginal_counts(trimmed, alphabet)
        trimmed_n = length(trimmed)
        freqs[r, :] .= counts ./ trimmed_n
        test = chisq_uniform_test(counts)
        stats[r] = test.statistic
        pvalues[r] = test.pvalue
        ess = categorical_effective_sample_size(trimmed, alphabet)
        effective_ns[r] = ess.effective_n
        iacts[r] = ess.integrated_autocorrelation_time
        ess_test = ess_corrected_chisq_uniform_test(counts, ess.effective_n)
        ess_stats[r] = ess_test.statistic
        ess_pvalues[r] = ess_test.pvalue
    end

    return (; full_freqs, freqs, stats, pvalues, effective_ns, iacts,
            ess_stats, ess_pvalues, trimmed_n)
end

function run_uniform_marginal_validation(; n::Int = 100_000,
                                           k::Int = 8,
                                           replicates::Int = 20,
                                           seed::Int = 20260630,
                                           trim_fraction::Real = 0.1)
    alphabet = Symbol.("s" .* string.(1:k))
    target = fill(1 / k, k)
    rows = NamedTuple[]
    histogram_rows = NamedTuple[]

    for (method, factory) in marginal_generator_factories()
        result = uniform_marginal_replicates(factory, n, alphabet;
                                             replicates, seed, trim_fraction)
        mean_freq = vec(mean(result.freqs; dims = 1))
        sd_freq = vec(std(result.freqs; dims = 1))
        full_deviations = abs.(result.full_freqs .- reshape(target, 1, :))
        deviations = abs.(result.freqs .- reshape(target, 1, :))
        pvalues = result.pvalues
        push!(rows, (;
            method,
            n,
            trimmed_n = result.trimmed_n,
            k,
            replicates,
            trim_fraction = Float64(trim_fraction),
            chi2_mean = mean(result.stats),
            chi2_max = maximum(result.stats),
            pvalue_min = minimum(pvalues),
            pvalue_median = median(pvalues),
            reject_005 = count(<(0.05), pvalues),
            effective_n_mean = mean(result.effective_ns),
            effective_n_min = minimum(result.effective_ns),
            iact_mean = mean(result.iacts),
            chi2_ess_mean = mean(result.ess_stats),
            pvalue_ess_min = minimum(result.ess_pvalues),
            pvalue_ess_median = median(result.ess_pvalues),
            reject_ess_005 = count(<(0.05), result.ess_pvalues),
            mean_full_total_variation =
                mean([total_variation(result.full_freqs[r, :], target)
                      for r in 1:replicates]),
            full_max_abs_error = maximum(full_deviations),
            mean_total_variation = mean([total_variation(result.freqs[r, :], target)
                                         for r in 1:replicates]),
            max_abs_error = maximum(deviations),
        ))
        for i in 1:k
            push!(histogram_rows, (;
                method,
                symbol = string(alphabet[i]),
                target = target[i],
                mean_frequency = mean_freq[i],
                sd_frequency = sd_freq[i],
            ))
        end
    end

    return rows, histogram_rows
end

function csv_escape(value)
    s = string(value)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

function write_csv(path::AbstractString, rows)
    isempty(rows) && throw(ArgumentError("cannot write CSV with no rows"))
    mkpath(dirname(path))
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(names, ","))
        for row in rows
            println(io, join((csv_escape(getproperty(row, name)) for name in names), ","))
        end
    end
    return path
end

function svg_escape(s)
    replace(replace(replace(string(s), "&" => "&amp;"), "<" => "&lt;"), ">" => "&gt;")
end

function write_uniform_histogram_svg(path::AbstractString, histogram_rows;
                                     k::Int, n::Int, trimmed_n::Int,
                                     replicates::Int)
    methods = unique(row.method for row in histogram_rows)
    rows_by_method = Dict(method => [row for row in histogram_rows if row.method == method]
                          for method in methods)
    panel_w, panel_h = 360, 250
    cols = 3
    rows = cld(length(methods), cols)
    width = cols * panel_w
    height = rows * panel_h + 90
    max_y = maximum(row.mean_frequency + row.sd_frequency for row in histogram_rows)
    max_y = max(max_y, 1 / k) * 1.18

    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<text x="$(width / 2)" y="32" text-anchor="middle" font-family="sans-serif" font-size="22">Uniform marginal validation, K=$k</text>""")
        println(io, """<text x="$(width / 2)" y="58" text-anchor="middle" font-family="sans-serif" font-size="14">n=$n, trimmed n=$trimmed_n after dropping first and last 10%, replicates=$replicates</text>""")

        for (midx, method) in enumerate(methods)
            c = (midx - 1) % cols
            r = (midx - 1) ÷ cols
            x0 = c * panel_w + 24
            y0 = r * panel_h + 82
            plot_left = x0 + 50
            plot_right = x0 + panel_w - 22
            plot_top = y0 + 38
            plot_bottom = y0 + panel_h - 48
            plot_w = plot_right - plot_left
            plot_h = plot_bottom - plot_top
            sy(v) = plot_bottom - v / max_y * plot_h

            println(io, """<text x="$(x0 + panel_w / 2)" y="$(y0 + 18)" text-anchor="middle" font-family="sans-serif" font-size="15">$(svg_escape(method))</text>""")
            println(io, """<line x1="$plot_left" y1="$plot_bottom" x2="$plot_right" y2="$plot_bottom" stroke="black"/>""")
            println(io, """<line x1="$plot_left" y1="$plot_top" x2="$plot_left" y2="$plot_bottom" stroke="black"/>""")
            target_y = round(sy(1 / k); digits = 2)
            println(io, """<line x1="$plot_left" y1="$target_y" x2="$plot_right" y2="$target_y" stroke="#d62728" stroke-width="2" stroke-dasharray="6 4"/>""")
            println(io, """<text x="$(plot_right - 4)" y="$(target_y - 5)" text-anchor="end" font-family="sans-serif" font-size="10" fill="#d62728">target</text>""")

            bar_gap = 5
            bar_w = (plot_w - (k + 1) * bar_gap) / k
            for (i, row) in enumerate(rows_by_method[method])
                x = plot_left + bar_gap + (i - 1) * (bar_w + bar_gap)
                y = sy(row.mean_frequency)
                h = plot_bottom - y
                println(io, """<rect x="$(round(x; digits=2))" y="$(round(y; digits=2))" width="$(round(bar_w; digits=2))" height="$(round(h; digits=2))" fill="#1f77b4"/>""")
                err_top = sy(row.mean_frequency + row.sd_frequency)
                err_bottom = sy(max(row.mean_frequency - row.sd_frequency, 0.0))
                cx = x + bar_w / 2
                println(io, """<line x1="$(round(cx; digits=2))" y1="$(round(err_top; digits=2))" x2="$(round(cx; digits=2))" y2="$(round(err_bottom; digits=2))" stroke="#333333" stroke-width="1"/>""")
                println(io, """<text x="$(round(cx; digits=2))" y="$(plot_bottom + 16)" text-anchor="middle" font-family="sans-serif" font-size="9">$(i)</text>""")
            end
            println(io, """<text x="$(plot_left - 8)" y="$(sy(1 / k) + 4)" text-anchor="end" font-family="sans-serif" font-size="10">$(round(1 / k; digits=3))</text>""")
            println(io, """<text x="$(x0 + panel_w / 2)" y="$(y0 + panel_h - 15)" text-anchor="middle" font-family="sans-serif" font-size="11">symbol index</text>""")
            println(io, """<text x="$(x0 + 14)" y="$(y0 + panel_h / 2)" transform="rotate(-90 $(x0 + 14) $(y0 + panel_h / 2))" text-anchor="middle" font-family="sans-serif" font-size="11">frequency</text>""")
        end
        println(io, """</svg>""")
    end
    svg_to_pdf_if_available(path)
    return path
end

function svg_to_pdf_if_available(svg_path::AbstractString)
    converter = Sys.which("rsvg-convert")
    pdf_path = replace(svg_path, r"\.svg$" => ".pdf")
    if converter === nothing
        inkscape = Sys.which("inkscape")
        inkscape === nothing && return nothing
        try
            run(`$inkscape $svg_path --export-type=pdf --export-filename=$pdf_path`)
        catch err
            @warn "Could not convert SVG validation plot to PDF with Inkscape" svg_path exception = err
            return nothing
        end
        return pdf_path
    end
    try
        run(`$converter -f pdf -o $pdf_path $svg_path`)
    catch err
        @warn "Could not convert SVG validation plot to PDF" svg_path exception = err
        return nothing
    end
    return pdf_path
end

function write_uniform_marginal_artifacts(; outdir::AbstractString = DEFAULT_MARGINAL_OUTDIR,
                                            kwargs...)
    rows, histogram_rows = run_uniform_marginal_validation(; kwargs...)
    mkpath(outdir)
    summary_path = write_csv(joinpath(outdir, "uniform_marginal_k8_summary.csv"), rows)
    histogram_table_path =
        write_csv(joinpath(outdir, "uniform_marginal_k8_histogram_data.csv"),
                  histogram_rows)
    first_row = first(rows)
    plot_path = write_uniform_histogram_svg(
        joinpath(outdir, "uniform_marginal_histograms_k8.svg"),
        histogram_rows;
        k = first_row.k,
        n = first_row.n,
        trimmed_n = first_row.trimmed_n,
        replicates = first_row.replicates,
    )
    provenance_path = joinpath(outdir, "README.md")
    open(provenance_path, "w") do io
        println(io, "# Uniform Marginal Validation")
        println(io)
        println(io, "Generated by `validation/marginal_control.jl`.")
        println(io)
        println(io, "- Target marginal: uniform categorical distribution.")
        println(io, "- Alphabet size: `k = $(first_row.k)`.")
        println(io, "- Sequence length: `n = $(first_row.n)`.")
        println(io, "- Trimmed length: `$(first_row.trimmed_n)` after dropping the first and last 10%.")
        println(io, "- Replicates per method: `$(first_row.replicates)`.")
        println(io)
        println(io, "The chi-squared p-values use the iid multinomial reference distribution.")
        println(io, "They are useful frequency diagnostics, but they are not exact tests for")
        println(io, "dependent LRD sequences. LRD can slow convergence of empirical")
        println(io, "marginals, so this reference is often conservative.")
        println(io)
        println(io, "The `effective_n_*`, `chi2_ess_*`, `pvalue_ess_*`, and")
        println(io, "`reject_ess_005` columns apply an approximate effective-sample-size")
        println(io, "correction. The correction estimates an integrated autocorrelation time")
        println(io, "for each centered one-hot symbol indicator, uses the smallest symbol")
        println(io, "effective sample size, and scales the chi-squared statistic by")
        println(io, "`effective_n / trimmed_n`. It is a dependence-aware diagnostic, not an")
        println(io, "exact categorical LRD test.")
        println(io)
        println(io, "The summary CSV reports both full-sequence marginal errors and")
        println(io, "trimmed-window marginal errors. The chi-squared diagnostics and")
        println(io, "histogram data use the trimmed window.")
        println(io)
        println(io, "More formal marginal tests should use a calibrated dependence-aware")
        println(io, "null, such as block/subsampling calibration or a parametric Monte Carlo")
        println(io, "envelope generated from the same configured method.")
    end
    return (; summary_path, histogram_table_path, plot_path, provenance_path)
end

function _iid_wavelet_markov(p, alphabet)
    k = length(alphabet)
    P = repeat(reshape(p, 1, k), k, 1)
    WaveletMarkov(0.8, alphabet, [P, P]; regime_weights = [0.5, 0.5])
end

function _marginal_lamp_transition(p)
    return lamp_repeat_transition(
        p; repeat_probability = MARGINAL_VALIDATION_LAMP_REPEAT_PROBABILITY)
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
    println()
    artifacts = write_uniform_marginal_artifacts()
    println("Wrote uniform marginal artifacts:")
    println("  ", artifacts.summary_path)
    println("  ", artifacts.histogram_table_path)
    println("  ", artifacts.plot_path)
end
