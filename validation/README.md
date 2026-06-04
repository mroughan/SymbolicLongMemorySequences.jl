# S5.jl Validation Studies

This folder contains reproducible simulation studies for generator controllability.
These are not LRD estimators. They test whether a generator respects user-facing
controls such as alphabet membership and target marginal probabilities.

The scripts use `StableRNGs.StableRNG` so results are reproducible across Julia
sessions and package updates.

## Marginal Control

Run from the package root:

```julia
julia --project=. validation/marginal_control.jl
```

The script prints aggregate total-variation and maximum absolute marginal errors for
`SpectralFGN`, `LGCM`, `LAMP`, `OnOffMarkov`, and `FSS` across a small grid of
sequence lengths, alphabet sizes, and marginal distributions.

`LGCM` is more expensive than the other methods because it calibrates latent
offsets over an `n × k` matrix. Increase `replicates`, `ns`, or
`calibration_iters` manually when running longer studies.

Keep generated data out of the repository unless a result table is intentionally being
tracked.
