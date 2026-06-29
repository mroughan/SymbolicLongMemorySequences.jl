# Changelog

All notable changes to S5.jl are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added
- `ARCHITECTURE.md` defining project goals, non-goals, boundaries, stable
  contracts, and the development pathway.
- `AGENTS.md` translating the architecture and repository guardrails into a
  required workflow for AI coding agents.
- Governance tests that preserve the core development-pathway sections and
  required agent references.
- Provenance-labeled local snapshots of the external development guidelines
  cited by the architecture.
- Separate GitHub Actions workflows for package tests, Aqua quality checks, JET
  static analysis, Codecov coverage uploads, and Documenter builds and deployment.
- Latest-Julia-only Aqua and JET quality scripts under `quality/`.
- `MarkovSpec`, a validated reusable first-order local-structure specification,
  with convenience constructors for `WaveletMarkov` and `OnOffMarkov`.
- `ControlCapabilities` and `control_capabilities(g)` for programmatic,
  strength-aware generator control contracts.
- A reproducible local-structure validation study and focused tests for
  per-regime Markov specifications and aggregate bigram control.
- README status badges for package tests, Aqua, JET, Codecov, and Documenter.
- `LocalStructureSpec` and `local_structure_order` as the extension path for
  future higher-order local-structure specifications.
- `benchmark/benchmarks.jl` and `benchmark/Project.toml` for BenchmarkTools.jl
  performance studies across implemented generators.
- Retained benchmark outputs under `benchmark/results/`, including CSV timings,
  histogram-style relative-time SVGs, sequence-length scaling SVGs, and
  `benchmark/RESULTS.md`.
- A rare `S5_BENCHMARK_SCALING=true` benchmark mode for sequence-length scaling
  studies over `n = 100` through `1_000_000`, using repeated independently
  seeded syntheses per BenchmarkTools trial and deferring the `k = 64` case.
- `generate_with_latent(g, n; rng)` for property-based generators, returning
  both the symbolic sequence and the numerical latent process used before
  symbolization.
- `VALIDATION_POLICY.md` describing fast tests, manual validation studies, and
  opt-in large validation and benchmark runs.
- Tests protecting the validation policy and benchmark infrastructure.
- Reusable LRD symbol-diagnostic helpers for centered one-hot transformations,
  autocorrelation, autocovariance, and periodogram conventions.
- A `validation/longmemory_comparison.jl` script comparing S5's formalized
  diagnostic helpers with LongMemory.jl's `autocovariance`, `autocorrelation`,
  and `periodogram` functions.
- Dashed interpretation-limit markers on LRD validation SVG plots, including the
  finite-sample `n / 10` lag scale and explicit generator memory limits where
  available.
- `lamp_repeat_transition` and transition-matrix support for `LAMP`, allowing
  repeat-biased identity/dyad transition patterns over history symbols.
- `DyadicLAMP` as MB1b, a scalable dyadic-bucket approximation to LAMP for large
  effective history depths.
- A selectable PB3 latent-driver path for `WaveletMarkov`, with `driver =
  :spectral` using spectral fGn rank-binning and `driver = :haar` retaining the
  original cascade for validation comparison.

### Changed
- Added explicit compatibility bounds for standard-library dependencies so Aqua
  can enforce complete dependency compatibility metadata.
- Fixed fresh-checkout Documenter CI setup by developing the unregistered S5 and
  IncCSV packages together before instantiating the documentation environment.
- Removed obsolete proposal-specific references from package documentation and
  docstrings.
- Documented that MB1/LAMP remains a finite-history approximation and that PB3's
  current Haar-style driver needs further validation or replacement.
- Changed `LAMP` generation so `d > n` uses observed history only and assigns
  missing pre-history weight to the target marginal instead of sampling an
  artificial random pre-history.
- Reclassified exact `LAMP` provenance as MB1a and `DyadicLAMP` provenance as
  MB1b.
- Refreshed LRD validation tables and SVG plots for the MB1a/MB1b split, and
  aligned benchmark labels with the same naming.
- Clarified README complexity notation for `n`, `d`, `k`, and `I` in the
  implemented-methods summary.
- Changed the default `WaveletMarkov` PB3 latent driver from the original
  Haar-style cascade to the spectral fGn rank-binning path, and recorded the
  driver choice in INC provenance metadata.
- Added nominal exact power-law reference overlays to autocorrelation and
  power-spectrum validation SVG plots.
- Added approximate short-range/asymptotic-onset markers to validation plots for
  generators with defensible onset scales, and documented MB4's weak
  power-spectrum behavior under the current finite Hawkes-style construction.
- Removed the now-redundant status column from README method summary tables.
- `HawkesSymbol` as MB4, a finite-history discrete-time Hawkes-style symbolic
  generator with power-law self/cross-excitation over observed history.
- Focused MB4 tests, INC provenance metadata, benchmark coverage, validation
  diagnostics, and documentation.
- Added `IntermittentMapSymbols` (PB4), `CalibratedAdditiveMarkov` (MB1c),
  and `DuplicationMutation` (MB5), with provenance metadata, tests,
  benchmarks, validation coverage, and documentation.
- Updated `ARCHITECTURE.md` to include PB4, MB1a/MB1b/MB1c, MB4, and MB5
  in the stable generator taxonomy.
- Regenerated LRD diagnostic SVG plots and log-binned plotting tables for the
  expanded PB1-PB4 and MB1a-MB5 method set.
- Revised MB5 `DuplicationMutation` to use a power-law copy-distance kernel
  instead of uniformly selected source blocks with power-law block lengths,
  after validation showed the earlier version produced nearly flat
  autocorrelation.
- Added the uniform factory API `method_ids`, `method_info`, and
  `make_generator` so users can list methods, inspect standard parameters, and
  construct standard generator cases through one entry point.
- Added `ParameterInfo` and `method_parameters` so factory methods can report
  accepted keyword inputs, defaults, domains, and short descriptions.
- Added a composable property-based API with `LatentSource`, `Symbolizer`, and
  `PropertyBasedGenerator`, plus spectral, Haar, intermittent-map, quantile,
  argmax, and Markov-regime building blocks for mix-and-match PB studies.
- Replaced separate LRD validation autocorrelation and power-spectrum SVGs with
  paired diagnostics plots that show both panels side by side for each method.
- Extended property-based validation plots and tables to include the numerical
  latent generator diagnostics alongside the final symbolic diagnostics.
- Expanded README, Documenter pages, validation notes, benchmark notes, and
  public docstrings with updated factory, PB4, MB1c, MB4, and MB5 examples.
- Reorganized Documenter documentation into `Home`, `API`, `Validation and
  Benchmarks`, and `Reference` pages, with the raw docstring reference moved to
  the dedicated `Reference` page.
- Removed `paper/` and `background/` from public version control while keeping
  local files in the working tree, and ignored both directories going forward.

## [0.1.0] — 2026-06-04

### Added
- `LRDGenerator` abstract type and `generate(g, n; rng)` common interface.
- **PB1 — `SpectralFGN`**: property-based generator using approximate spectral
  (FFT) synthesis of fractional Gaussian noise followed by quantization to symbols.
  Hurst parameter `H ∈ (0.5, 1.0)`. O(n log n) time.
- **PB2 — `LGCM`**: property-based Latent Gaussian Categorical Model using one
  latent fGn stream per symbol and calibrated argmax offsets for marginal control.
  Hurst parameter `H ∈ (0.5, 1.0)`.
- **PB3 — `WaveletMarkov`**: property-based multiscale Haar-like latent driver
  selecting Markov regimes with per-regime transition matrices.
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
- Documentation motivations for estimator testing, excess entropy and entropy-rate
  experiments, non-language LLM-style sequence training, context-length diagnostics,
  anomaly/change detection, and privacy-preserving simulation.
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
