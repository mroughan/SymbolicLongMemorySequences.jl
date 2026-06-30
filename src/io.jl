"""
    save_sequence(filepath, seq, gen; created = string(today())) -> filepath

Write a generated symbol sequence to an INC file (IncCSV.jl format) with full
provenance metadata.

The INC file contains:
- A metadata block recording the SymbolicLongMemorySequences.jl package version, the generator type and
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
        "title"     => "SymbolicLongMemorySequences.jl synthetic LRD symbol sequence",
        "package"   => "SymbolicLongMemorySequences",
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
        "title"     => "SymbolicLongMemorySequences.jl synthetic LRD symbol sequence",
        "package"   => "SymbolicLongMemorySequences",
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

function _build_metadata(gen::WaveletMarkov, n::Int, created::String)
    k   = length(gen.alphabet)
    R   = length(gen.transition_matrices)
    alp = join(string.(gen.alphabet), ",")
    rw  = join(string.(round.(gen.regime_weights; digits = 8)), ",")
    Dict(
        "title"     => "SymbolicLongMemorySequences.jl synthetic LRD symbol sequence",
        "package"   => "SymbolicLongMemorySequences",
        "version"   => string(pkgversion(@__MODULE__)),
        "created"   => created,
        "n"         => n,
        "generator" => "WaveletMarkov",
        "method"    => "PB3",
        "generator_params" => Dict(
            "H"              => string(gen.H),
            "alphabet_size"  => k,
            "alphabet"       => alp,
            "n_regimes"      => R,
            "regime_weights" => rw,
            "cascade_depth"  => gen.cascade_depth,
            "driver"         => string(gen.driver),
        ),
        "columns" => Dict(
            "index"  => "time index (1-based)",
            "symbol" => "generated symbol from the alphabet",
        ),
    )
end

function _build_metadata(gen::IntermittentMapSymbols, n::Int, created::String)
    k   = length(gen.alphabet)
    alp = join(string.(gen.alphabet), ",")
    mar = join(string.(round.(gen.marginal; digits = 8)), ",")
    Dict(
        "title"     => "SymbolicLongMemorySequences.jl synthetic LRD symbol sequence",
        "package"   => "SymbolicLongMemorySequences",
        "version"   => string(pkgversion(@__MODULE__)),
        "created"   => created,
        "n"         => n,
        "generator" => "IntermittentMapSymbols",
        "method"    => "PB4",
        "generator_params" => Dict(
            "z"             => string(gen.z),
            "burnin"        => gen.burnin,
            "alphabet_size" => k,
            "alphabet"      => alp,
            "marginal"      => mar,
            "latent_driver" => "Pomeau-Manneville-style intermittent map",
        ),
        "columns" => Dict(
            "index"  => "time index (1-based)",
            "symbol" => "generated symbol from the alphabet",
        ),
    )
end

function _build_metadata(gen::PropertyBasedGenerator, n::Int, created::String)
    sym = gen.symbolizer
    alp = join(string.(_symbolizer_alphabet(sym)), ",")
    Dict(
        "title"     => "SymbolicLongMemorySequences.jl synthetic LRD symbol sequence",
        "package"   => "SymbolicLongMemorySequences",
        "version"   => string(pkgversion(@__MODULE__)),
        "created"   => created,
        "n"         => n,
        "generator" => "PropertyBasedGenerator",
        "method"    => "PB-composed",
        "generator_params" => Dict(
            "source"        => string(nameof(typeof(gen.source))),
            "symbolizer"    => string(nameof(typeof(sym))),
            "latent_width"  => latent_width(sym),
            "alphabet_size" => length(_symbolizer_alphabet(sym)),
            "alphabet"      => alp,
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
    P = join((join(string.(round.(gen.transition_matrix[i, :]; digits = 8)), ",")
              for i in axes(gen.transition_matrix, 1)), ";")
    Dict(
        "title"     => "SymbolicLongMemorySequences.jl synthetic LRD symbol sequence",
        "package"   => "SymbolicLongMemorySequences",
        "version"   => string(pkgversion(@__MODULE__)),
        "created"   => created,
        "n"         => n,
        "generator" => "LAMP",
        "method"    => "MB1a",
        "generator_params" => Dict(
            "beta"          => string(gen.beta),
            "d"             => gen.d,
            "alphabet_size" => k,
            "alphabet"      => alp,
            "marginal"      => mar,
            "transition_matrix" => P,
        ),
        "columns" => Dict(
            "index"  => "time index (1-based)",
            "symbol" => "generated symbol from the alphabet",
        ),
    )
end

function _build_metadata(gen::CalibratedAdditiveMarkov, n::Int, created::String)
    k   = length(gen.alphabet)
    alp = join(string.(gen.alphabet), ",")
    mar = join(string.(round.(gen.marginal; digits = 8)), ",")
    Dict(
        "title"     => "SymbolicLongMemorySequences.jl synthetic LRD symbol sequence",
        "package"   => "SymbolicLongMemorySequences",
        "version"   => string(pkgversion(@__MODULE__)),
        "created"   => created,
        "n"         => n,
        "generator" => "CalibratedAdditiveMarkov",
        "method"    => "MB1c",
        "generator_params" => Dict(
            "beta"          => string(gen.beta),
            "d"             => gen.d,
            "strength"      => string(gen.strength),
            "alphabet_size" => k,
            "alphabet"      => alp,
            "marginal"      => mar,
            "memory_function" => "centered additive power law",
        ),
        "columns" => Dict(
            "index"  => "time index (1-based)",
            "symbol" => "generated symbol from the alphabet",
        ),
    )
end

function _build_metadata(gen::DyadicLAMP, n::Int, created::String)
    k   = length(gen.alphabet)
    alp = join(string.(gen.alphabet), ",")
    mar = join(string.(round.(gen.marginal; digits = 8)), ",")
    P = join((join(string.(round.(gen.transition_matrix[i, :]; digits = 8)), ",")
              for i in axes(gen.transition_matrix, 1)), ";")
    Dict(
        "title"     => "SymbolicLongMemorySequences.jl synthetic LRD symbol sequence",
        "package"   => "SymbolicLongMemorySequences",
        "version"   => string(pkgversion(@__MODULE__)),
        "created"   => created,
        "n"         => n,
        "generator" => "DyadicLAMP",
        "method"    => "MB1b",
        "generator_params" => Dict(
            "beta"          => string(gen.beta),
            "d"             => gen.d,
            "alphabet_size" => k,
            "alphabet"      => alp,
            "marginal"      => mar,
            "transition_matrix" => P,
            "history_representation" => "dyadic buckets",
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
        "title"     => "SymbolicLongMemorySequences.jl synthetic LRD symbol sequence",
        "package"   => "SymbolicLongMemorySequences",
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
        "title"     => "SymbolicLongMemorySequences.jl synthetic LRD symbol sequence",
        "package"   => "SymbolicLongMemorySequences",
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

function _build_metadata(gen::HawkesSymbol, n::Int, created::String)
    k = length(gen.alphabet)
    alp = join(string.(gen.alphabet), ",")
    base = join(string.(round.(gen.baseline; digits = 8)), ",")
    E = join((join(string.(round.(gen.excitation[i, :]; digits = 8)), ",")
              for i in axes(gen.excitation, 1)), ";")
    Dict(
        "title"     => "SymbolicLongMemorySequences.jl synthetic LRD symbol sequence",
        "package"   => "SymbolicLongMemorySequences",
        "version"   => string(pkgversion(@__MODULE__)),
        "created"   => created,
        "n"         => n,
        "generator" => "HawkesSymbol",
        "method"    => "MB4",
        "generator_params" => Dict(
            "beta"          => string(gen.beta),
            "d"             => gen.d,
            "c"             => string(gen.c),
            "alphabet_size" => k,
            "alphabet"      => alp,
            "baseline"      => base,
            "excitation"    => E,
            "kernel"        => "normalized power law",
            "time_model"    => "discrete",
        ),
        "columns" => Dict(
            "index"  => "time index (1-based)",
            "symbol" => "generated symbol from the alphabet",
        ),
    )
end

function _build_metadata(gen::DuplicationMutation, n::Int, created::String)
    k = length(gen.alphabet)
    alp = join(string.(gen.alphabet), ",")
    mar = join(string.(round.(gen.marginal; digits = 8)), ",")
    Dict(
        "title"     => "SymbolicLongMemorySequences.jl synthetic LRD symbol sequence",
        "package"   => "SymbolicLongMemorySequences",
        "version"   => string(pkgversion(@__MODULE__)),
        "created"   => created,
        "n"         => n,
        "generator" => "DuplicationMutation",
        "method"    => "MB5",
        "generator_params" => Dict(
            "alpha"                => string(gen.alpha),
            "mutation_probability" => string(gen.mutation_probability),
            "seed_length"          => gen.seed_length,
            "max_block_length"     => gen.max_block_length,
            "alphabet_size"        => k,
            "alphabet"             => alp,
            "marginal"             => mar,
            "growth_model"         => "power-law lag copy-and-mutate growth",
        ),
        "columns" => Dict(
            "index"  => "time index (1-based)",
            "symbol" => "generated symbol from the alphabet",
        ),
    )
end
