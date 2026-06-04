target_marginal(g::Union{SpectralFGN,LGCM,LAMP}) = copy(Float64.(g.marginal))
function target_marginal(g::OnOffMarkov)
    regime_π = stationary_distribution(g.switching_matrix)
    p = zeros(Float64, length(g.alphabet))
    for (r, P) in enumerate(g.transition_matrices)
        p .+= regime_π[r] .* stationary_distribution(P)
    end
    p ./= sum(p)
    return p
end
target_marginal(g::FSS) = Float64.(g.rates) ./ sum(g.rates)
