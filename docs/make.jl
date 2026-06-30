using Documenter
using SymbolicLongMemorySequences

makedocs(
    sitename = "SymbolicLongMemorySequences.jl",
    modules  = [SymbolicLongMemorySequences],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = nothing,
    ),
    pages = [
        "Home"                      => "index.md",
        "API"                       => "api.md",
        "Validation and Benchmarks" => "validation_benchmarks.md",
        "Reference"                 => "reference.md",
    ],
    checkdocs = :exports,
    warnonly  = false,
)

deploydocs(
    repo = "github.com/mroughan/SymbolicLongMemorySequences.jl.git",
    devbranch = "main",
)
