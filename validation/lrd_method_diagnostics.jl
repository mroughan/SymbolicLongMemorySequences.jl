import IncCSV
using S5
using StableRNGs
using Statistics
using LinearAlgebra: I

include(joinpath(@__DIR__, "lrd_symbol_diagnostics.jl"))

const DEFAULT_N = 100_000
const DEFAULT_REPLICATES = 30
const DEFAULT_SEED = 20260605
const DEFAULT_ALPHABET = ['A', 'B', 'C', 'D', 'E']

function generator_factories(alphabet)
    k = length(alphabet)
    marginal = fill(1 / k, k)

    regime_matrices = observable_regime_matrices(k; dominance = 0.72, persistence = 0.35)
    Q = fill(1.0 / (k - 1), k, k)
    for i in 1:k
        Q[i, i] = 0.0
    end

    return [
        ("PB1_SpectralFGN", () -> SpectralFGN(0.8, alphabet, marginal)),
        ("PB2_LGCM", () -> LGCM(0.8, alphabet, marginal; calibration_iters = 8)),
        ("PB3_WaveletMarkov_Haar", () ->
            WaveletMarkov(0.8, alphabet, regime_matrices; driver = :haar)),
        ("PB3_WaveletMarkov_Spectral", () ->
            WaveletMarkov(0.8, alphabet, regime_matrices; driver = :spectral)),
        ("MB1a_LAMP", () -> LAMP(0.4, alphabet, marginal; d = 200, epsilon = 0.02)),
        ("MB1b_DyadicLAMP", () -> DyadicLAMP(0.4, alphabet, marginal;
                                             d = 100_000, epsilon = 0.02)),
        ("MB2_OnOffMarkov", () -> OnOffMarkov(1.4, alphabet, regime_matrices, Q; L_min = 50.0)),
        ("MB3_FSS", () -> FSS(1.4, alphabet; rates = ones(k))),
        ("MB4_HawkesSymbol", () -> HawkesSymbol(0.4, alphabet;
            baseline = fill(1.0, k),
            excitation = 6.0 .* Matrix{Float64}(I, k, k),
            d = 20_000)),
    ]
end

function stationary_transition_matrix(p::AbstractVector{<:Real}, persistence::Real)
    p = Float64.(p)
    ρ = Float64(persistence)
    k = length(p)
    P = Matrix{Float64}(undef, k, k)
    for i in 1:k, j in 1:k
        P[i, j] = (i == j ? ρ : 0.0) + (1 - ρ) * p[j]
    end
    return P
end

function observable_regime_matrices(k::Int; dominance::Real, persistence::Real)
    matrices = Matrix{Float64}[]
    for r in 1:k
        p = fill((1 - dominance) / (k - 1), k)
        p[r] = dominance
        push!(matrices, stationary_transition_matrix(p, persistence))
    end
    return matrices
end

function logbin(x::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                bins::Int = 220, positive_y::Bool = true)
    pairs = [(Float64(xi), Float64(yi)) for (xi, yi) in zip(x, y)
             if xi > 0 && (!positive_y || yi > 0) && isfinite(xi) && isfinite(yi)]
    isempty(pairs) && return Float64[], Float64[]

    lx = log10.(first.(pairs))
    lo, hi = minimum(lx), maximum(lx)
    edges = range(lo, hi; length = bins + 1)
    bx = Float64[]
    by = Float64[]

    for b in 1:bins
        left, right = edges[b], edges[b + 1]
        vals = if b == bins
            [(xv, yv) for (xv, yv) in pairs if left ≤ log10(xv) ≤ right]
        else
            [(xv, yv) for (xv, yv) in pairs if left ≤ log10(xv) < right]
        end
        isempty(vals) && continue
        push!(bx, 10.0^mean(log10.(first.(vals))))
        push!(by, mean(last.(vals)))
    end
    order = sortperm(bx)
    bx = bx[order]
    by = by[order]
    return bx, by
end

function append_rows!(rows, method, xname, yname, x, y)
    for (xi, yi) in zip(x, y)
        push!(rows, (;
            method,
            Symbol(xname) => xi,
            Symbol(yname) => yi,
        ))
    end
end

function svg_escape(s)
    replace(replace(replace(string(s), "&" => "&amp;"), "<" => "&lt;"), ">" => "&gt;")
end

function diagnostic_lag_limit(n::Int; overlap_fraction::Real = 0.9)
    0 < overlap_fraction < 1 ||
        throw(ArgumentError("overlap_fraction must be in (0, 1)"))
    n ≥ 2 || throw(ArgumentError("n must be at least 2"))
    return max(1, floor(Int, n * (1 - Float64(overlap_fraction))))
end

function intrinsic_lag_limit(g)
    return nothing
end

intrinsic_lag_limit(g::LAMP) = g.d
intrinsic_lag_limit(g::DyadicLAMP) = g.d
intrinsic_lag_limit(g::HawkesSymbol) = g.d

function acf_limit_annotations(g, n::Int)
    finite_limit = diagnostic_lag_limit(n)
    annotations = [(;
        x = finite_limit,
        label = "finite-sample limit n/10",
        color = "#d62728",
    )]
    intrinsic = intrinsic_lag_limit(g)
    if intrinsic !== nothing && intrinsic > 0 && intrinsic != finite_limit
        push!(annotations, (;
            x = intrinsic,
            label = "generator memory limit",
            color = "#9467bd",
        ))
    end
    return annotations
end

function spectrum_limit_annotations(g, n::Int)
    annotations = NamedTuple[]
    for ann in acf_limit_annotations(g, n)
        push!(annotations, (;
            x = 1 / ann.x,
            label = ann.label,
            color = ann.color,
        ))
    end
    return annotations
end

nominal_acf_decay_exponent(g::Union{SpectralFGN,LGCM,WaveletMarkov}) = 2 - 2g.H
nominal_acf_decay_exponent(g::Union{LAMP,DyadicLAMP}) = g.beta
nominal_acf_decay_exponent(g::Union{OnOffMarkov,FSS}) = g.alpha - 1
nominal_acf_decay_exponent(g::HawkesSymbol) = g.beta

function nominal_reference_line(x::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                                exponent::Real, label::AbstractString)
    pts = [(Float64(xi), Float64(yi)) for (xi, yi) in zip(x, y)
           if xi > 0 && yi > 0 && isfinite(xi) && isfinite(yi)]
    isempty(pts) && return NamedTuple[]
    sort!(pts; by = first)
    x0, y0 = first(pts)
    exp = Float64(exponent)
    ref = [(xi, y0 * (xi / x0)^exp) for (xi, _) in pts]
    return [(;
        x = first.(ref),
        y = last.(ref),
        label,
        color = "#555555",
        dash = "4 4",
    )]
end

function acf_power_law_reference(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, g)
    beta = nominal_acf_decay_exponent(g)
    return nominal_reference_line(x, y;
        exponent = -beta,
        label = "nominal ACF beta=$(round(beta; digits = 3))")
end

function spectrum_power_law_reference(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, g)
    beta = nominal_acf_decay_exponent(g)
    return nominal_reference_line(x, y;
        exponent = beta - 1,
        label = "nominal spectrum beta=$(round(beta; digits = 3))")
end

function write_loglog_svg(path, title, xlabel, ylabel, x, y;
                          vertical_lines = NamedTuple[],
                          reference_lines = NamedTuple[])
    pts = [(Float64(xi), Float64(yi)) for (xi, yi) in zip(x, y)
           if xi > 0 && yi > 0 && isfinite(xi) && isfinite(yi)]
    sort!(pts; by = first)
    isempty(pts) && throw(ArgumentError("cannot write log-log SVG with no positive points"))
    ref_pts = Tuple{Float64, Float64}[]
    for ref in reference_lines
        append!(ref_pts, [(Float64(xi), Float64(yi)) for (xi, yi) in zip(ref.x, ref.y)
                          if xi > 0 && yi > 0 && isfinite(xi) && isfinite(yi)])
    end
    width, height = 900, 650
    ml, mr, mt, mb = 90, 30, 55, 80
    all_pts = isempty(ref_pts) ? pts : vcat(pts, ref_pts)
    xmin, xmax = extrema(first.(all_pts))
    ymin, ymax = extrema(last.(all_pts))
    lxmin, lxmax = log10(xmin), log10(xmax)
    lymin, lymax = log10(ymin), log10(ymax)
    lymin == lymax && (lymax += 1)

    sx(v) = ml + (log10(v) - lxmin) / (lxmax - lxmin) * (width - ml - mr)
    sy(v) = height - mb - (log10(v) - lymin) / (lymax - lymin) * (height - mt - mb)

    poly = join(["$(round(sx(xi); digits=2)),$(round(sy(yi); digits=2))"
                 for (xi, yi) in pts], " ")

    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<text x="$(width/2)" y="28" text-anchor="middle" font-family="sans-serif" font-size="20">$(svg_escape(title))</text>""")
        println(io, """<line x1="$ml" y1="$(height-mb)" x2="$(width-mr)" y2="$(height-mb)" stroke="black"/>""")
        println(io, """<line x1="$ml" y1="$mt" x2="$ml" y2="$(height-mb)" stroke="black"/>""")

        for e in floor(Int, lxmin):ceil(Int, lxmax)
            val = 10.0^e
            xmin ≤ val ≤ xmax || continue
            xpix = sx(val)
            println(io, """<line x1="$xpix" y1="$(height-mb)" x2="$xpix" y2="$(height-mb+5)" stroke="black"/>""")
            println(io, """<text x="$xpix" y="$(height-mb+25)" text-anchor="middle" font-family="sans-serif" font-size="12">10^$e</text>""")
        end
        for e in floor(Int, lymin):ceil(Int, lymax)
            val = 10.0^e
            ymin ≤ val ≤ ymax || continue
            ypix = sy(val)
            println(io, """<line x1="$(ml-5)" y1="$ypix" x2="$ml" y2="$ypix" stroke="black"/>""")
            println(io, """<text x="$(ml-10)" y="$(ypix+4)" text-anchor="end" font-family="sans-serif" font-size="12">10^$e</text>""")
        end

        legend_y = mt + 18
        for ref in reference_lines
            rpts = [(Float64(xi), Float64(yi)) for (xi, yi) in zip(ref.x, ref.y)
                    if xi > 0 && yi > 0 && isfinite(xi) && isfinite(yi)]
            isempty(rpts) && continue
            sort!(rpts; by = first)
            rpoly = join(["$(round(sx(xi); digits=2)),$(round(sy(yi); digits=2))"
                          for (xi, yi) in rpts], " ")
            color = get(ref, :color, "#555555")
            dash = get(ref, :dash, "4 4")
            label = svg_escape(get(ref, :label, "power-law reference"))
            println(io, """<polyline fill="none" stroke="$color" stroke-width="2" stroke-dasharray="$dash" points="$rpoly"/>""")
            println(io, """<line x1="$(width-300)" y1="$legend_y" x2="$(width-260)" y2="$legend_y" stroke="$color" stroke-width="2" stroke-dasharray="$dash"/>""")
            println(io, """<text x="$(width-252)" y="$(legend_y+4)" font-family="sans-serif" font-size="12">$(label)</text>""")
            legend_y += 18
        end

        println(io, """<polyline fill="none" stroke="#1f77b4" stroke-width="2" points="$poly"/>""")
        for line in vertical_lines
            xv = Float64(line.x)
            xmin ≤ xv ≤ xmax || continue
            xpix = round(sx(xv); digits = 2)
            color = get(line, :color, "#d62728")
            label = svg_escape(get(line, :label, "diagnostic limit"))
            println(io, """<line x1="$xpix" y1="$mt" x2="$xpix" y2="$(height-mb)" stroke="$color" stroke-width="2" stroke-dasharray="7 5"/>""")
            println(io, """<line x1="$(width-300)" y1="$legend_y" x2="$(width-260)" y2="$legend_y" stroke="$color" stroke-width="2" stroke-dasharray="7 5"/>""")
            println(io, """<text x="$(width-252)" y="$(legend_y+4)" font-family="sans-serif" font-size="12">$(label)</text>""")
            legend_y += 18
        end
        println(io, """<text x="$(width/2)" y="$(height-25)" text-anchor="middle" font-family="sans-serif" font-size="16">$(svg_escape(xlabel))</text>""")
        println(io, """<text x="25" y="$(height/2)" transform="rotate(-90 25 $(height/2))" text-anchor="middle" font-family="sans-serif" font-size="16">$(svg_escape(ylabel))</text>""")
        println(io, """</svg>""")
    end
end

function write_diagnostic_inc(path, rows, title, columns)
    metadata = Dict(
        "title" => title,
        "package" => "S5",
        "created_by" => "validation/lrd_method_diagnostics.jl",
        "columns" => columns,
    )
    IncCSV.writeinc(path, rows; metadata)
    return path
end

function run_lrd_diagnostics(; n::Int = DEFAULT_N,
                               replicates::Int = DEFAULT_REPLICATES,
                               seed::Int = DEFAULT_SEED,
                               outdir::AbstractString = "validation/results/lrd_diagnostics",
                               save_sequences::Bool = true)
    alphabet = DEFAULT_ALPHABET
    mkpath(outdir)
    seqdir = joinpath(outdir, "sequences")
    plotdir = joinpath(outdir, "plots")
    mkpath(plotdir)
    save_sequences && mkpath(seqdir)

    acf_rows = NamedTuple[]
    pxx_rows = NamedTuple[]
    acf_plot_rows = NamedTuple[]
    pxx_plot_rows = NamedTuple[]

    for (method_index, (method, factory)) in enumerate(generator_factories(alphabet))
        println("Running $method")
        maxlag = n ÷ 2
        avg_acf = zeros(Float64, maxlag)
        avg_power = zeros(Float64, n ÷ 2)
        freqs = nothing
        diagnostic_gen = nothing

        for r in 1:replicates
            rng = StableRNG(seed + 10_000 * method_index + r)
            gen = factory()
            r == 1 && (diagnostic_gen = gen)
            seq = generate(gen, n; rng)
            if save_sequences
                seqpath = joinpath(seqdir, "$(method)_$(lpad(r, 2, '0')).inc")
                save_sequence(seqpath, seq, gen)
            end
            acf, f, pxx = indicator_diagnostics(seq, alphabet; maxlag)
            avg_acf .+= acf
            avg_power .+= pxx
            freqs = f
            println("  replicate $r / $replicates")
        end

        avg_acf ./= replicates
        avg_power ./= replicates
        lags = collect(1:maxlag)
        b_lags, b_acf = logbin(lags, avg_acf; positive_y = true)
        b_freqs, b_power = logbin(freqs, avg_power; positive_y = true)

        append_rows!(acf_rows, method, "lag", "autocorrelation", lags, avg_acf)
        append_rows!(pxx_rows, method, "frequency", "power", freqs, avg_power)
        append_rows!(acf_plot_rows, method, "lag", "autocorrelation", b_lags, b_acf)
        append_rows!(pxx_plot_rows, method, "frequency", "power", b_freqs, b_power)

        write_loglog_svg(joinpath(plotdir, "$(method)_autocorrelation.svg"),
                         "$method average autocorrelation",
                         "lag", "positive average autocorrelation",
                         b_lags, b_acf;
                         vertical_lines = acf_limit_annotations(diagnostic_gen, n),
                         reference_lines = acf_power_law_reference(b_lags, b_acf,
                                                                   diagnostic_gen))
        write_loglog_svg(joinpath(plotdir, "$(method)_power_spectrum.svg"),
                         "$method average power spectrum",
                         "frequency", "average periodogram",
                         b_freqs, b_power;
                         vertical_lines = spectrum_limit_annotations(diagnostic_gen, n),
                         reference_lines = spectrum_power_law_reference(b_freqs,
                                                                        b_power,
                                                                        diagnostic_gen))
    end

    write_diagnostic_inc(
        joinpath(outdir, "average_autocorrelation.inc"),
        acf_rows,
        "Average signed one-hot symbol autocorrelation by generator",
        Dict(
            "method" => "S5 generator name",
            "lag" => "integer lag",
            "autocorrelation" => "signed mean one-hot autocorrelation averaged across symbols and replicates",
        ),
    )
    write_diagnostic_inc(
        joinpath(outdir, "average_power_spectrum.inc"),
        pxx_rows,
        "Average one-hot symbol power spectrum by generator",
        Dict(
            "method" => "S5 generator name",
            "frequency" => "Fourier frequency after log-binning",
            "power" => "mean periodogram averaged across symbols and replicates",
        ),
    )
    write_diagnostic_inc(
        joinpath(outdir, "plot_autocorrelation_logbins.inc"),
        acf_plot_rows,
        "Log-binned positive average autocorrelation used for SVG plots",
        Dict(
            "method" => "S5 generator name",
            "lag" => "geometric mean lag in log bin",
            "autocorrelation" => "positive binned mean autocorrelation",
        ),
    )
    write_diagnostic_inc(
        joinpath(outdir, "plot_power_spectrum_logbins.inc"),
        pxx_plot_rows,
        "Log-binned average power spectrum used for SVG plots",
        Dict(
            "method" => "S5 generator name",
            "frequency" => "geometric mean Fourier frequency in log bin",
            "power" => "binned mean periodogram",
        ),
    )

    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    outdir = run_lrd_diagnostics()
    println("Wrote diagnostics to $outdir")
end
