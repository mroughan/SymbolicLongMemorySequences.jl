# S5.jl Validation Studies

This folder contains reproducible simulation studies for generator controllability.
These are not LRD estimators. They test whether a generator respects user-facing
controls such as alphabet membership and target marginal probabilities.

The scripts use `StableRNGs.StableRNG` so results are reproducible across Julia
sessions and package updates.

See `../VALIDATION_POLICY.md` for the project validation tiers. Fast contract tests
belong in `test/`; broader empirical evidence belongs here and should be run
manually or behind explicit `S5_` environment flags.

## Marginal Control

Run from the package root:

```julia
julia --project=. validation/marginal_control.jl
```

The script prints aggregate total-variation and maximum absolute marginal errors for
`SpectralFGN`, `LGCM`, `WaveletMarkov`, `LAMP`, `OnOffMarkov`, and `FSS` across
a small grid of sequence lengths, alphabet sizes, and marginal distributions.

`LGCM` is more expensive than the other methods because it calibrates latent
offsets over an `n × k` matrix. Increase `replicates`, `ns`, or
`calibration_iters` manually when running longer studies.

Keep generated data out of the repository unless a result table is intentionally being
tracked.

Large validation grids should be opt-in. Prefer keyword arguments in the script API
and environment variables such as `S5_VALIDATION_LARGE=true` or
`S5_VALIDATION_REPLICATES=<integer>` when a script grows beyond the default manual
run.

## Local Structure Control

Run from the package root:

```julia
julia --project=. validation/local_structure.jl
```

The script measures row-wise and aggregate transition-matrix total variation for
`WaveletMarkov` and `OnOffMarkov` across iid, persistent, and cyclic first-order
Markov specifications. Each method uses identical `MarkovSpec` values in every
regime, so the common transition matrix is an unambiguous aggregate bigram target.
This study does not claim that mixtures of different regime specifications have a
simple aggregate target.

## LRD Method Diagnostics

Run from the package root:

```julia
julia --project=. validation/lrd_method_diagnostics.jl
```

The diagnostic transformations and numerical conventions are formalized in
`lrd_symbol_diagnostics.jl`. Symbol sequences are transformed into one centered
one-hot numeric series per alphabet symbol:

```julia
x_t = 1{X_t = symbol} - mean(1{X_t = symbol})
```

Centering removes the symbol marginal so autocorrelation and spectrum summaries
focus on dependence rather than the zero-frequency mass. Zero-variance indicator
series are skipped because autocorrelation is undefined for constant series.

This generates 30 sequences of length 100,000 for each implemented method using the
alphabet `{A,B,C,D,E}`. Sequences are saved as INC files under
`validation/results/lrd_diagnostics/sequences/`. That sequence directory is ignored
by Git because these large generated files are reproducible and change whenever the
diagnostic configuration changes.

For PB3 (`WaveletMarkov`) and MB2 (`OnOffMarkov`), the diagnostic uses five
observable regimes whose stationary distributions are each biased toward one
alphabet symbol (`dominance = 0.72`) with moderate within-regime Markov
persistence (`persistence = 0.35`). The regimes are balanced overall, but the
symbol-level one-hot ACF and spectrum can see the long-memory regime process.
MB2 uses `L_min = 50.0` so the heavy-tailed sojourn mechanism is visible at
`n = 100_000`. If all regimes share the same stationary marginal, these
diagnostics can look short-memory even when the latent regime process has
long-range structure.

The script computes the one-hot symbol autocorrelation and power spectrum for each
sequence, averages across symbols and replicates, and writes:

- `average_autocorrelation.inc`: full signed averaged autocorrelation by lag;
- `average_power_spectrum.inc`: full averaged periodogram by Fourier frequency;
- `plot_autocorrelation_logbins.inc`: positive log-binned values used for plotting;
- `plot_power_spectrum_logbins.inc`: log-binned spectrum values used for plotting.

The log-binning code keeps the final bin closed only at its upper edge and sorts
the binned x-values before writing tables and SVG polylines, so plot x-values are
strictly increasing within each method. Log-log SVG plots are written under
`validation/results/lrd_diagnostics/plots/`.

## LongMemory.jl Comparison

S5.jl does not depend on estimator packages at runtime, but validation can compare
the formalized diagnostic transformations with
[`LongMemory.jl`](https://github.com/everval/LongMemory.jl). Instantiate the
validation environment, then run:

```julia
julia --project=validation -e 'using Pkg; Pkg.instantiate()'
julia --project=validation validation/longmemory_comparison.jl
```

The comparison script adapts LongMemory.jl conventions explicitly:

- `autocovariance(x, k)` and `autocorrelation(x, k)` return lags `0:k-1`;
- S5 plots use lags `1:maxlag`, so the lag-zero autocorrelation is dropped;
- LongMemory.jl periodograms report angular frequencies and include zero;
- S5 plots use cycles per observation and drop zero frequency.

The script reports maximum absolute differences between S5's local
LongMemory-compatible helpers and LongMemory.jl's exported `autocovariance`,
`autocorrelation`, and `periodogram` functions on the centered one-hot series.

## Benchmarks

Benchmarks are separate from validation studies and live under `../benchmark/`.
Run the default suite from the package root:

```julia
julia --project=benchmark benchmark/benchmarks.jl
```

Run the larger opt-in suite with:

```julia
S5_BENCHMARK_LARGE=true julia --project=benchmark benchmark/benchmarks.jl
```
