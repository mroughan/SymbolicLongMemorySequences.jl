target_marginal(g::Union{SpectralFGN,LAMP}) = copy(Float64.(g.marginal))
target_marginal(g::FSS) = Float64.(g.rates) ./ sum(g.rates)
