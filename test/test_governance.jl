@testset "Development governance" begin
    root = normpath(joinpath(@__DIR__, ".."))

    architecture = read(joinpath(root, "ARCHITECTURE.md"), String)
    for heading in ("## Goals", "## Non-Goals", "## Philosophy",
                    "## Development Pathway", "## Testing And Validation")
        @test contains(architecture, heading)
    end

    agents = read(joinpath(root, "AGENTS.md"), String)
    for required_path in ("../GUARDRAILS.md", "ARCHITECTURE.md", "CHANGELOG.md",
                          "validation/")
        @test contains(agents, required_path)
    end

    references = read(joinpath(root, "references", "README.md"), String)
    @test contains(references, "Downloaded on:")
    @test contains(references, "canonical online sources remain authoritative")
    for filename in ("julia-style-guide.html", "julia-performance-tips.html",
                     "semantic-versioning.html", "keep-a-changelog.html",
                     "aqua-readme.md", "jet-readme.md",
                     "documenter-hosting.md", "codecov-action-readme.md",
                     "longmemory-api-notes.md")
        @test isfile(joinpath(root, "references", filename))
    end
end
