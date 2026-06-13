using Documenter
using S5

makedocs(
    sitename = "S5.jl",
    modules  = [S5],
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
    repo = "github.com/mroughan/S5.jl.git",
    devbranch = "main",
)
