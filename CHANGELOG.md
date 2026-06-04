# Changelog

All notable changes to S5.jl are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [0.1.0] — 2026-06-04

### Added
- `LRDGenerator` abstract type and `generate(g, n; rng)` common interface.
- **PB1 — `SpectralFGN`**: property-based generator using approximate spectral
  (FFT) synthesis of fractional Gaussian noise followed by quantization to symbols.
  Hurst parameter `H ∈ (0.5, 1.0)`. O(n log n) time.
- **PB2 — `LGCM`**: property-based Latent Gaussian Categorical Model using one
  latent fGn stream per symbol and calibrated argmax offsets for marginal control.
  Hurst parameter `H ∈ (0.5, 1.0)`.
- **MB1 — `LAMP`**: model-based Linear-Additive Markov Process with power-law
  history weights `wₖ ∝ k^{-(1+β)}`. ACF decay exponent `β ∈ (0, 1)`. O(n·d) time.
- **MB2 — `OnOffMarkov`**: model-based heavy-tailed regime-switching Markov chain
  with Pareto sojourns and per-regime transition matrices.
- **MB3 — `FSS`**: model-based Fractal Symbol Sequence via independent
  Pareto-distributed renewal processes, one per symbol, merged in event-time order.
  Tail index `α ∈ (1, 2)`, `H = (3−α)/2`. O(n·k) time.
- `save_sequence(filepath, seq, gen)`: writes generated sequences to INC format
  (IncCSV.jl) with full provenance metadata (package version, generator type and
  parameters, creation date).
- `quantize_to_symbols`: utility for mapping real-valued sequences to symbol
  alphabets via sample-quantile thresholding.
- Full test suite covering output correctness, argument validation, statistical
  properties, and INC round-trip.
- Documenter.jl documentation with API reference.
- GitHub Actions CI: tests on Julia 1.10 and latest stable.

### Changed
- `quantize_to_symbols` now uses deterministic rank binning with integer counts
  from `bin_counts`, giving finite-sample marginals as close as possible to the
  requested target.
- `LAMP` now supports an `epsilon` marginal innovation term to improve
  finite-sample marginal control and avoid finite-history absorption.
- `FSS` Pareto renewal draws now use Distributions.jl.
- Tests and validation studies use StableRNGs.jl for reproducible simulation.

### Validation
- Marginal/local-structure helpers: `target_marginal`, `empirical_marginal`,
  `empirical_bigram`, `empirical_trigram`, `total_variation`, and
  `rowwise_total_variation`.
- Markov helpers: `validate_transition_matrix` and `stationary_distribution`.
- `validation/` folder with a reproducible Monte Carlo marginal-control study.
