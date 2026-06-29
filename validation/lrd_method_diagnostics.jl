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
        ("PB4_IntermittentMapSymbols", () ->
            IntermittentMapSymbols(1.6, alphabet, marginal; burnin = 1000)),
        ("MB1a_LAMP", () -> LAMP(0.4, alphabet, marginal; d = 200, epsilon = 0.02)),
        ("MB1b_DyadicLAMP", () -> DyadicLAMP(0.4, alphabet, marginal;
                                             d = 100_000, epsilon = 0.02)),
        ("MB1c_CalibratedAdditiveMarkov", () ->
            CalibratedAdditiveMarkov(0.4, alphabet, marginal; d = 200,
                                     strength = 0.75)),
        ("MB2_OnOffMarkov", () -> OnOffMarkov(1.4, alphabet, regime_matrices, Q; L_min = 50.0)),
        ("MB3_FSS", () -> FSS(1.4, alphabet; rates = ones(k))),
        ("MB4_HawkesSymbol", () -> HawkesSymbol(0.4, alphabet;
            baseline = fill(1.0, k),
            excitation = 6.0 .* Matrix{Float64}(I, k, k),
            d = 20_000)),
        ("MB5_DuplicationMutation", () ->
            DuplicationMutation(1.4, alphabet, marginal;
                                mutation_probability = 0.02,
                                seed_length = 128,
                                max_block_length = 20_000)),
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
intrinsic_lag_limit(g::CalibratedAdditiveMarkov) = g.d
intrinsic_lag_limit(g::HawkesSymbol) = g.d
intrinsic_lag_limit(g::DuplicationMutation) = g.max_block_length

function asymptotic_lag_threshold(g; slope_fraction::Real = 0.9)
    return nothing
end

asymptotic_lag_threshold(g::OnOffMarkov; slope_fraction::Real = 0.9) =
    max(1, ceil(Int, g.L_min))

function asymptotic_lag_threshold(g::HawkesSymbol; slope_fraction::Real = 0.9)
    0 < slope_fraction < 1 ||
        throw(ArgumentError("slope_fraction must be in (0, 1)"))

    # For an offset kernel (lag + c)^(-beta), the log-log slope magnitude is
    # beta * lag / (lag + c). Mark the lag where it reaches the requested
    # fraction of its asymptotic value beta.
    return max(1, ceil(Int, slope_fraction * g.c / (1 - slope_fraction)))
end

function acf_limit_annotations(g, n::Int)
    finite_limit = diagnostic_lag_limit(n)
    annotations = [(;
        x = finite_limit,
        label = "finite-sample limit n/10",
        color = "#d62728",
    )]
    threshold = asymptotic_lag_threshold(g)
    if threshold !== nothing && threshold > 1 && threshold != finite_limit
        push!(annotations, (;
            x = threshold,
            label = "approx. power-law onset",
            color = "#2ca02c",
        ))
    end
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
nominal_acf_decay_exponent(g::IntermittentMapSymbols) = 2 - g.z
nominal_acf_decay_exponent(g::Union{LAMP,DyadicLAMP,CalibratedAdditiveMarkov}) = g.beta
nominal_acf_decay_exponent(g::Union{OnOffMarkov,FSS}) = g.alpha - 1
nominal_acf_decay_exponent(g::HawkesSymbol) = g.beta
nominal_acf_decay_exponent(g::DuplicationMutation) = g.alpha - 1

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

function _loglog_panel_data(x, y, reference_lines)
    pts = [(Float64(xi), Float64(yi)) for (xi, yi) in zip(x, y)
           if xi > 0 && yi > 0 && isfinite(xi) && isfinite(yi)]
    sort!(pts; by = first)
    isempty(pts) && throw(ArgumentError("cannot write log-log panel with no positive points"))
    ref_pts = Tuple{Float64, Float64}[]
    for ref in reference_lines
        append!(ref_pts, [(Float64(xi), Float64(yi)) for (xi, yi) in zip(ref.x, ref.y)
                          if xi > 0 && yi > 0 && isfinite(xi) && isfinite(yi)])
    end
    all_pts = isempty(ref_pts) ? pts : vcat(pts, ref_pts)
    xmin, xmax = extrema(first.(all_pts))
    ymin, ymax = extrema(last.(all_pts))
    lxmin, lxmax = log10(xmin), log10(xmax)
    lymin, lymax = log10(ymin), log10(ymax)
    lymin == lymax && (lymax += 1)
    return (; pts, xmin, xmax, ymin, ymax, lxmin, lxmax, lymin, lymax)
end

function _draw_loglog_panel(io, panel, title, xlabel, ylabel, x, y;
                            vertical_lines = NamedTuple[],
                            reference_lines = NamedTuple[])
    data = _loglog_panel_data(x, y, reference_lines)
    x0, y0, width, height = panel.x, panel.y, panel.width, panel.height
    ml, mr, mt, mb = 78, 22, 48, 68
    plot_left = x0 + ml
    plot_right = x0 + width - mr
    plot_top = y0 + mt
    plot_bottom = y0 + height - mb
    plot_width = plot_right - plot_left
    plot_height = plot_bottom - plot_top

    sx(v) = plot_left + (log10(v) - data.lxmin) / (data.lxmax - data.lxmin) * plot_width
    sy(v) = plot_bottom - (log10(v) - data.lymin) / (data.lymax - data.lymin) * plot_height

    poly = join(["$(round(sx(xi); digits=2)),$(round(sy(yi); digits=2))"
                 for (xi, yi) in data.pts], " ")

    println(io, """<text x="$(x0 + width / 2)" y="$(y0 + 25)" text-anchor="middle" font-family="sans-serif" font-size="18">$(svg_escape(title))</text>""")
    println(io, """<line x1="$plot_left" y1="$plot_bottom" x2="$plot_right" y2="$plot_bottom" stroke="black"/>""")
    println(io, """<line x1="$plot_left" y1="$plot_top" x2="$plot_left" y2="$plot_bottom" stroke="black"/>""")

    for e in floor(Int, data.lxmin):ceil(Int, data.lxmax)
        val = 10.0^e
        data.xmin ≤ val ≤ data.xmax || continue
        xpix = sx(val)
        println(io, """<line x1="$xpix" y1="$plot_bottom" x2="$xpix" y2="$(plot_bottom + 5)" stroke="black"/>""")
        println(io, """<text x="$xpix" y="$(plot_bottom + 24)" text-anchor="middle" font-family="sans-serif" font-size="11">10^$e</text>""")
    end
    for e in floor(Int, data.lymin):ceil(Int, data.lymax)
        val = 10.0^e
        data.ymin ≤ val ≤ data.ymax || continue
        ypix = sy(val)
        println(io, """<line x1="$(plot_left - 5)" y1="$ypix" x2="$plot_left" y2="$ypix" stroke="black"/>""")
        println(io, """<text x="$(plot_left - 9)" y="$(ypix + 4)" text-anchor="end" font-family="sans-serif" font-size="11">10^$e</text>""")
    end

    legend_x = plot_right - 235
    legend_y = plot_top + 16
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
        println(io, """<line x1="$legend_x" y1="$legend_y" x2="$(legend_x + 34)" y2="$legend_y" stroke="$color" stroke-width="2" stroke-dasharray="$dash"/>""")
        println(io, """<text x="$(legend_x + 42)" y="$(legend_y + 4)" font-family="sans-serif" font-size="11">$(label)</text>""")
        legend_y += 17
    end

    println(io, """<polyline fill="none" stroke="#1f77b4" stroke-width="2" points="$poly"/>""")
    for line in vertical_lines
        xv = Float64(line.x)
        data.xmin ≤ xv ≤ data.xmax || continue
        xpix = round(sx(xv); digits = 2)
        color = get(line, :color, "#d62728")
        label = svg_escape(get(line, :label, "diagnostic limit"))
        println(io, """<line x1="$xpix" y1="$plot_top" x2="$xpix" y2="$plot_bottom" stroke="$color" stroke-width="2" stroke-dasharray="7 5"/>""")
        println(io, """<line x1="$legend_x" y1="$legend_y" x2="$(legend_x + 34)" y2="$legend_y" stroke="$color" stroke-width="2" stroke-dasharray="7 5"/>""")
        println(io, """<text x="$(legend_x + 42)" y="$(legend_y + 4)" font-family="sans-serif" font-size="11">$(label)</text>""")
        legend_y += 17
    end
    println(io, """<text x="$(x0 + width / 2)" y="$(y0 + height - 20)" text-anchor="middle" font-family="sans-serif" font-size="14">$(svg_escape(xlabel))</text>""")
    println(io, """<text x="$(x0 + 22)" y="$(y0 + height / 2)" transform="rotate(-90 $(x0 + 22) $(y0 + height / 2))" text-anchor="middle" font-family="sans-serif" font-size="14">$(svg_escape(ylabel))</text>""")
end

function write_diagnostic_pair_svg(path, method, acf_x, acf_y, spectrum_x, spectrum_y;
                                   acf_vertical_lines = NamedTuple[],
                                   spectrum_vertical_lines = NamedTuple[],
                                   acf_reference_lines = NamedTuple[],
                                   spectrum_reference_lines = NamedTuple[])
    width, height = 1500, 650
    left = (; x = 20, y = 58, width = 720, height = 570)
    right = (; x = 760, y = 58, width = 720, height = 570)
    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<text x="$(width / 2)" y="30" text-anchor="middle" font-family="sans-serif" font-size="22">$(svg_escape(method)) validation diagnostics</text>""")
        _draw_loglog_panel(io, left, "Autocorrelation", "lag",
                           "positive average autocorrelation", acf_x, acf_y;
                           vertical_lines = acf_vertical_lines,
                           reference_lines = acf_reference_lines)
        _draw_loglog_panel(io, right, "Power spectrum", "frequency",
                           "average periodogram", spectrum_x, spectrum_y;
                           vertical_lines = spectrum_vertical_lines,
                           reference_lines = spectrum_reference_lines)
        println(io, """</svg>""")
    end
    svg_to_pdf_if_available(path)
    return path
end

function write_property_diagnostic_svg(path, method,
                                       latent_acf_x, latent_acf_y,
                                       latent_spectrum_x, latent_spectrum_y,
                                       symbol_acf_x, symbol_acf_y,
                                       symbol_spectrum_x, symbol_spectrum_y;
                                       acf_vertical_lines = NamedTuple[],
                                       spectrum_vertical_lines = NamedTuple[],
                                       acf_reference_lines = NamedTuple[],
                                       spectrum_reference_lines = NamedTuple[])
    width, height = 1500, 1180
    panels = (
        (; x = 20, y = 58, width = 720, height = 520),
        (; x = 760, y = 58, width = 720, height = 520),
        (; x = 20, y = 628, width = 720, height = 520),
        (; x = 760, y = 628, width = 720, height = 520),
    )
    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<text x="$(width / 2)" y="30" text-anchor="middle" font-family="sans-serif" font-size="22">$(svg_escape(method)) latent and symbol validation diagnostics</text>""")
        _draw_loglog_panel(io, panels[1], "Latent autocorrelation", "lag",
                           "positive latent autocorrelation",
                           latent_acf_x, latent_acf_y;
                           vertical_lines = acf_vertical_lines,
                           reference_lines = acf_reference_lines)
        _draw_loglog_panel(io, panels[2], "Latent power spectrum", "frequency",
                           "latent periodogram",
                           latent_spectrum_x, latent_spectrum_y;
                           vertical_lines = spectrum_vertical_lines,
                           reference_lines = spectrum_reference_lines)
        _draw_loglog_panel(io, panels[3], "Symbol autocorrelation", "lag",
                           "positive average one-hot autocorrelation",
                           symbol_acf_x, symbol_acf_y;
                           vertical_lines = acf_vertical_lines,
                           reference_lines = acf_reference_lines)
        _draw_loglog_panel(io, panels[4], "Symbol power spectrum", "frequency",
                           "average one-hot periodogram",
                           symbol_spectrum_x, symbol_spectrum_y;
                           vertical_lines = spectrum_vertical_lines,
                           reference_lines = spectrum_reference_lines)
        println(io, """</svg>""")
    end
    svg_to_pdf_if_available(path)
    return path
end

function svg_to_pdf_if_available(svg_path::AbstractString)
    converter = Sys.which("rsvg-convert")
    converter === nothing && return nothing
    pdf_path = replace(svg_path, r"\.svg$" => ".pdf")
    try
        run(`$converter -f pdf -o $pdf_path $svg_path`)
    catch err
        @warn "Could not convert SVG validation plot to PDF" svg_path exception = err
        return nothing
    end
    return pdf_path
end

latent_diagnostics_available(::Any) = false
latent_diagnostics_available(::Union{SpectralFGN,LGCM,WaveletMarkov,
                                     IntermittentMapSymbols,
                                     PropertyBasedGenerator}) = true

function latent_numeric_diagnostics(latent::AbstractMatrix{<:Real};
                                    maxlag::Int = size(latent, 2) ÷ 2)
    series = Vector{Float64}[]
    for row in axes(latent, 1)
        x = Float64.(latent[row, :])
        x .-= mean(x)
        mean(abs2, x) > 0 && push!(series, x)
    end
    isempty(series) && throw(ArgumentError("all latent series have zero variance"))

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
    for path in readdir(plotdir; join = true)
        if endswith(path, "_autocorrelation.svg") || endswith(path, "_power_spectrum.svg")
            rm(path; force = true)
        end
    end

    acf_rows = NamedTuple[]
    pxx_rows = NamedTuple[]
    latent_acf_rows = NamedTuple[]
    latent_pxx_rows = NamedTuple[]
    acf_plot_rows = NamedTuple[]
    pxx_plot_rows = NamedTuple[]
    latent_acf_plot_rows = NamedTuple[]
    latent_pxx_plot_rows = NamedTuple[]

    for (method_index, (method, factory)) in enumerate(generator_factories(alphabet))
        println("Running $method")
        maxlag = n ÷ 2
        avg_acf = zeros(Float64, maxlag)
        avg_power = zeros(Float64, n ÷ 2)
        avg_latent_acf = zeros(Float64, maxlag)
        avg_latent_power = zeros(Float64, n ÷ 2)
        freqs = nothing
        latent_freqs = nothing
        diagnostic_gen = nothing
        has_latent = false

        for r in 1:replicates
            rng = StableRNG(seed + 10_000 * method_index + r)
            gen = factory()
            r == 1 && (diagnostic_gen = gen)
            if latent_diagnostics_available(gen)
                seq, latent = generate_with_latent(gen, n; rng)
                latent_acf, latent_f, latent_pxx =
                    latent_numeric_diagnostics(latent; maxlag)
                avg_latent_acf .+= latent_acf
                avg_latent_power .+= latent_pxx
                latent_freqs = latent_f
                has_latent = true
            else
                seq = generate(gen, n; rng)
            end
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
        if has_latent
            avg_latent_acf ./= replicates
            avg_latent_power ./= replicates
        end
        lags = collect(1:maxlag)
        b_lags, b_acf = logbin(lags, avg_acf; positive_y = true)
        b_freqs, b_power = logbin(freqs, avg_power; positive_y = true)

        append_rows!(acf_rows, method, "lag", "autocorrelation", lags, avg_acf)
        append_rows!(pxx_rows, method, "frequency", "power", freqs, avg_power)
        append_rows!(acf_plot_rows, method, "lag", "autocorrelation", b_lags, b_acf)
        append_rows!(pxx_plot_rows, method, "frequency", "power", b_freqs, b_power)
        if has_latent
            b_latent_lags, b_latent_acf =
                logbin(lags, avg_latent_acf; positive_y = true)
            b_latent_freqs, b_latent_power =
                logbin(latent_freqs, avg_latent_power; positive_y = true)
            append_rows!(latent_acf_rows, method, "lag", "autocorrelation",
                         lags, avg_latent_acf)
            append_rows!(latent_pxx_rows, method, "frequency", "power",
                         latent_freqs, avg_latent_power)
            append_rows!(latent_acf_plot_rows, method, "lag", "autocorrelation",
                         b_latent_lags, b_latent_acf)
            append_rows!(latent_pxx_plot_rows, method, "frequency", "power",
                         b_latent_freqs, b_latent_power)

            write_property_diagnostic_svg(
                joinpath(plotdir, "$(method)_diagnostics.svg"),
                method,
                b_latent_lags,
                b_latent_acf,
                b_latent_freqs,
                b_latent_power,
                b_lags,
                b_acf,
                b_freqs,
                b_power;
                acf_vertical_lines = acf_limit_annotations(diagnostic_gen, n),
                spectrum_vertical_lines = spectrum_limit_annotations(diagnostic_gen, n),
                acf_reference_lines = acf_power_law_reference(b_lags, b_acf,
                                                              diagnostic_gen),
                spectrum_reference_lines = spectrum_power_law_reference(b_freqs,
                                                                        b_power,
                                                                        diagnostic_gen),
            )
        else
            write_diagnostic_pair_svg(
                joinpath(plotdir, "$(method)_diagnostics.svg"),
                method,
                b_lags,
                b_acf,
                b_freqs,
                b_power;
                acf_vertical_lines = acf_limit_annotations(diagnostic_gen, n),
                spectrum_vertical_lines = spectrum_limit_annotations(diagnostic_gen, n),
                acf_reference_lines = acf_power_law_reference(b_lags, b_acf,
                                                              diagnostic_gen),
                spectrum_reference_lines = spectrum_power_law_reference(b_freqs,
                                                                        b_power,
                                                                        diagnostic_gen),
            )
        end
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
    write_diagnostic_inc(
        joinpath(outdir, "latent_average_autocorrelation.inc"),
        latent_acf_rows,
        "Average latent numerical autocorrelation by property-based generator",
        Dict(
            "method" => "S5 property-based generator name",
            "lag" => "integer lag",
            "autocorrelation" => "mean latent autocorrelation averaged across latent streams and replicates",
        ),
    )
    write_diagnostic_inc(
        joinpath(outdir, "latent_average_power_spectrum.inc"),
        latent_pxx_rows,
        "Average latent numerical power spectrum by property-based generator",
        Dict(
            "method" => "S5 property-based generator name",
            "frequency" => "Fourier frequency",
            "power" => "mean latent periodogram averaged across latent streams and replicates",
        ),
    )
    write_diagnostic_inc(
        joinpath(outdir, "latent_plot_autocorrelation_logbins.inc"),
        latent_acf_plot_rows,
        "Log-binned latent autocorrelation used for property-based SVG plots",
        Dict(
            "method" => "S5 property-based generator name",
            "lag" => "geometric mean lag in log bin",
            "autocorrelation" => "positive binned mean latent autocorrelation",
        ),
    )
    write_diagnostic_inc(
        joinpath(outdir, "latent_plot_power_spectrum_logbins.inc"),
        latent_pxx_plot_rows,
        "Log-binned latent power spectrum used for property-based SVG plots",
        Dict(
            "method" => "S5 property-based generator name",
            "frequency" => "geometric mean Fourier frequency in log bin",
            "power" => "binned mean latent periodogram",
        ),
    )

    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    outdir = run_lrd_diagnostics()
    println("Wrote diagnostics to $outdir")
end
