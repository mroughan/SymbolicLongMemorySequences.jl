target_marginal(g::Union{SpectralFGN,LGCM,LAMP,DyadicLAMP}) = copy(Float64.(g.marginal))
function target_marginal(g::WaveletMarkov)
    p = zeros(Float64, length(g.alphabet))
    for (r, P) in enumerate(g.transition_matrices)
        p .+= g.regime_weights[r] .* stationary_distribution(P)
    end
    p ./= sum(p)
    return p
end
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
target_marginal(g::HawkesSymbol) = Float64.(g.baseline) ./ sum(g.baseline)

control_capabilities(::SpectralFGN) =
    ControlCapabilities(:exact, :finite_sample, :induced, :induced, :approximate)
control_capabilities(::LGCM) =
    ControlCapabilities(:exact, :empirical, :induced, :induced, :latent_approximate)
control_capabilities(::WaveletMarkov) =
    ControlCapabilities(:exact, :implied, :per_regime, :induced, :latent_approximate)
control_capabilities(::LAMP) =
    ControlCapabilities(:exact, :innovation_target, :induced, :induced, :finite_history)
control_capabilities(::DyadicLAMP) =
    ControlCapabilities(:exact, :innovation_target, :induced, :induced, :finite_history)
control_capabilities(::OnOffMarkov) =
    ControlCapabilities(:exact, :implied, :per_regime, :induced, :nominal)
control_capabilities(::FSS) =
    ControlCapabilities(:exact, :asymptotic, :induced, :induced, :nominal)
control_capabilities(::HawkesSymbol) =
    ControlCapabilities(:exact, :implied, :induced, :induced, :finite_history)
