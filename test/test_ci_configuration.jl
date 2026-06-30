@testset "CI configuration" begin
    root = normpath(joinpath(@__DIR__, ".."))
    workflows = joinpath(root, ".github", "workflows")

    expected = ("ci.yml", "aqua.yml", "jet.yml", "codecov.yml",
                "documentation.yml")
    @test all(isfile(joinpath(workflows, file)) for file in expected)

    for file in ("aqua.yml", "jet.yml")
        workflow = read(joinpath(workflows, file), String)
        @test contains(workflow, "version: '1'")
        @test !contains(workflow, "matrix:")
    end

    coverage = read(joinpath(workflows, "codecov.yml"), String)
    @test contains(coverage, "julia-processcoverage@v1")
    @test contains(coverage, "codecov/codecov-action@v7")

    documentation = read(joinpath(workflows, "documentation.yml"), String)
    @test contains(documentation, "docs/make.jl")
    @test contains(documentation, "Pkg.develop([")
    @test contains(documentation, "PackageSpec(path=pwd())")
    @test contains(documentation,
                   "PackageSpec(url=\"https://github.com/mroughan/IncCSV.jl\")")
    @test !contains(documentation, "Pkg.add(url=")

    makefile = read(joinpath(root, "docs", "make.jl"), String)
    @test contains(makefile, "deploydocs(")
    @test contains(makefile, "devbranch = \"main\"")

    readme = read(joinpath(root, "README.md"), String)
    for badge in ("actions/workflows/ci.yml/badge.svg",
                  "actions/workflows/aqua.yml/badge.svg",
                  "actions/workflows/jet.yml/badge.svg",
                  "codecov.io/gh/mroughan/SymbolicLongMemorySequences.jl/branch/main/graph/badge.svg",
                  "actions/workflows/documentation.yml/badge.svg")
        @test contains(readme, badge)
    end
end
