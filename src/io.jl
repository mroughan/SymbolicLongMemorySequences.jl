"""
    save_sequence(filepath, seq, gen; created = string(today())) -> filepath

Write a generated symbol sequence to an INC file (IncCSV.jl format) with full
provenance metadata.

The INC file contains:
- A metadata block recording the S5.jl package version, the generator type and
  all its parameters, and the creation date.
- A two-column CSV body: `index` (1-based integer) and `symbol` (string).

# Arguments
- `filepath::AbstractString`: output path (`.inc` extension recommended).
- `seq::AbstractVector`: symbol sequence as returned by [`generate`](@ref).
- `gen::LRDGenerator`: the generator instance used to produce `seq`.

# Keyword Arguments
- `created::String`: creation date (default: today's date in ISO 8601 format).

# Returns
`filepath`, to allow chaining.

# Examples
```julia
julia> g   = SpectralFGN(0.8, [:a, :b, :c])
julia> seq = generate(g, 1000)
julia> save_sequence("output.inc", seq, g)
"output.inc"
```
"""
function save_sequence(filepath::AbstractString, seq::AbstractVector,
                       gen::LRDGenerator;
                       created::String = string(today()))
    rows = [(; index = i, symbol = string(seq[i])) for i in eachindex(seq)]
    meta = _build_metadata(gen, length(seq), created)
    IncCSV.writeinc(filepath, rows; metadata = meta)
    return filepath
end

# ── Per-generator metadata builders ───────────────────────────────────────────

function _build_metadata(gen::SpectralFGN, n::Int, created::String)
    k   = length(gen.alphabet)
    alp = join(string.(gen.alphabet), ",")
    mar = join(string.(round.(gen.marginal; digits = 8)), ",")
    Dict(
        "title"     => "S5.jl synthetic LRD symbol sequence",
        "package"   => "S5",
        "version"   => string(pkgversion(@__MODULE__)),
        "created"   => created,
        "n"         => n,
        "generator" => "SpectralFGN",
        "method"    => "PB1",
        "generator_params" => Dict(
            "H"             => string(gen.H),
            "alphabet_size" => k,
            "alphabet"      => alp,
            "marginal"      => mar,
        ),
        "columns" => Dict(
            "index"  => "time index (1-based)",
            "symbol" => "generated symbol from the alphabet",
        ),
    )
end

function _build_metadata(gen::LGCM, n::Int, created::String)
    k   = length(gen.alphabet)
    alp = join(string.(gen.alphabet), ",")
    mar = join(string.(round.(gen.marginal; digits = 8)), ",")
    Dict(
        "title"     => "S5.jl synthetic LRD symbol sequence",
        "package"   => "S5",
        "version"   => string(pkgversion(@__MODULE__)),
        "created"   => created,
        "n"         => n,
        "generator" => "LGCM",
        "method"    => "PB2",
        "generator_params" => Dict(
            "H"                 => string(gen.H),
            "alphabet_size"     => k,
            "alphabet"          => alp,
            "marginal"          => mar,
            "calibration_iters" => gen.calibration_iters,
            "calibration_rate"  => string(gen.calibration_rate),
        ),
        "columns" => Dict(
            "index"  => "time index (1-based)",
            "symbol" => "generated symbol from the alphabet",
        ),
    )
end

function _build_metadata(gen::LAMP, n::Int, created::String)
    k   = length(gen.alphabet)
    alp = join(string.(gen.alphabet), ",")
    mar = join(string.(round.(gen.marginal; digits = 8)), ",")
    Dict(
        "title"     => "S5.jl synthetic LRD symbol sequence",
        "package"   => "S5",
        "version"   => string(pkgversion(@__MODULE__)),
        "created"   => created,
        "n"         => n,
        "generator" => "LAMP",
        "method"    => "MB1",
        "generator_params" => Dict(
            "beta"          => string(gen.beta),
            "d"             => gen.d,
            "alphabet_size" => k,
            "alphabet"      => alp,
            "marginal"      => mar,
        ),
        "columns" => Dict(
            "index"  => "time index (1-based)",
            "symbol" => "generated symbol from the alphabet",
        ),
    )
end

function _build_metadata(gen::OnOffMarkov, n::Int, created::String)
    k   = length(gen.alphabet)
    R   = length(gen.transition_matrices)
    alp = join(string.(gen.alphabet), ",")
    H   = (3 - gen.alpha) / 2
    Dict(
        "title"     => "S5.jl synthetic LRD symbol sequence",
        "package"   => "S5",
        "version"   => string(pkgversion(@__MODULE__)),
        "created"   => created,
        "n"         => n,
        "generator" => "OnOffMarkov",
        "method"    => "MB2",
        "generator_params" => Dict(
            "alpha"         => string(gen.alpha),
            "H_nominal"     => string(round(H; digits = 8)),
            "L_min"         => string(gen.L_min),
            "alphabet_size" => k,
            "alphabet"      => alp,
            "n_regimes"     => R,
        ),
        "columns" => Dict(
            "index"  => "time index (1-based)",
            "symbol" => "generated symbol from the alphabet",
        ),
    )
end

function _build_metadata(gen::FSS, n::Int, created::String)
    k    = length(gen.alphabet)
    alp  = join(string.(gen.alphabet), ",")
    rats = join(string.(round.(gen.rates; digits = 8)), ",")
    H    = (3 - gen.alpha) / 2
    Dict(
        "title"     => "S5.jl synthetic LRD symbol sequence",
        "package"   => "S5",
        "version"   => string(pkgversion(@__MODULE__)),
        "created"   => created,
        "n"         => n,
        "generator" => "FSS",
        "method"    => "MB3",
        "generator_params" => Dict(
            "alpha"         => string(gen.alpha),
            "H_nominal"     => string(round(H; digits = 8)),
            "x_min"         => string(gen.x_min),
            "alphabet_size" => k,
            "alphabet"      => alp,
            "rates"         => rats,
        ),
        "columns" => Dict(
            "index"  => "time index (1-based)",
            "symbol" => "generated symbol from the alphabet",
        ),
    )
end
