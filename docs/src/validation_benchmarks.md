# Validation and Benchmarks

S5.jl separates package correctness from empirical research evidence. Fast
tests protect public contracts; validation scripts provide reproducible
scientific diagnostics; benchmarks provide machine-specific performance
evidence.

See also `VALIDATION_POLICY.md`, `validation/README.md`, and
`benchmark/README.md` in the repository.

## Fast Tests

Run the package test suite from the repository root:

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```

These tests are the main development path. They cover deterministic contracts
such as constructor validation, output length, element type, alphabet membership,
reproducible RNG use, factory discovery, and INC provenance metadata. They also
include small seeded statistical checks where exact assertions are not
meaningful.

## Manual Validation Studies

Manual validation scripts live in `validation/`. They use reproducible RNGs and
may write aggregate result tables or SVG/PDF plots.

```julia
julia --project=. validation/marginal_control.jl
julia --project=. validation/local_structure.jl
julia --project=. validation/lrd_method_diagnostics.jl
julia --project=validation validation/longmemory_comparison.jl
```

Validation scripts may use `make_generator` for smoke tests of standard cases,
but studies that support method-specific claims should prefer explicit
constructors so the scientific parameters and assumptions remain visible.

## Marginal Control

`validation/marginal_control.jl` compares empirical frequencies against declared
targets across sequence lengths, alphabet sizes, and marginal distributions. It
reports aggregate total-variation and maximum absolute marginal errors.

This study is designed to test user-facing marginal controls, not to estimate a
Hurst parameter.

## Local Structure Control

`validation/local_structure.jl` measures row-wise and aggregate
transition-matrix total variation for `WaveletMarkov` and `OnOffMarkov` under
validated `MarkovSpec` values.

The study uses identical `MarkovSpec` values in each regime when it needs an
unambiguous aggregate bigram target. Mixtures of different regimes may be useful
generators, but they do not generally have a simple aggregate Markov target.

## LRD Method Diagnostics

`validation/lrd_method_diagnostics.jl` creates one-hot symbol diagnostics for all
implemented methods. Property-based methods additionally report diagnostics for
the numerical latent process returned by `generate_with_latent`. Each symbol
sequence is transformed into centered indicator series before autocorrelation,
autocovariance, or periodogram calculations:

```julia
x_t = 1{X_t = symbol} - mean(1{X_t = symbol})
```

Centering removes the symbol marginal so the summaries focus on dependence
rather than zero-frequency mass. Zero-variance indicator series are skipped.

The script writes:

- `average_autocorrelation.inc`;
- `average_power_spectrum.inc`;
- `latent_average_autocorrelation.inc`;
- `latent_average_power_spectrum.inc`;
- `plot_autocorrelation_logbins.inc`;
- `plot_power_spectrum_logbins.inc`;
- `latent_plot_autocorrelation_logbins.inc`;
- `latent_plot_power_spectrum_logbins.inc`;
- paired SVG plots under `validation/results/lrd_diagnostics/plots/`, with
  autocorrelation on the left and power spectrum on the right. If
  `rsvg-convert` is available, PDF copies are written for paper inclusion.

For property-based methods the SVG/PDF plots have four panels: latent
autocorrelation, latent power spectrum, symbolic autocorrelation, and symbolic
power spectrum. This makes the transformation loss or distortion visible.

Autocorrelation plots include a vertical dashed finite-sample interpretation
limit at lag `n / 10`. Methods with explicit finite memory, such as `LAMP`,
`DyadicLAMP`, `CalibratedAdditiveMarkov`, `HawkesSymbol`, and
`DuplicationMutation`, also mark the generator cutoff. Power-spectrum plots show
the reciprocal scales. Methods with a defensible asymptotic-onset scale also
mark an approximate power-law onset: `OnOffMarkov` uses its Pareto scale `L_min`,
and `HawkesSymbol` uses the lag where the local log-log slope of
`(lag + c)^(-beta)` reaches 90% of its asymptotic slope. Gray dashed reference
lines show the nominal power-law slope implied by each generator's configured
decay parameter; they are visual guides, not fitted curves and not proofs of LRD
behavior.

The MB4 (`HawkesSymbol`) plots are a current cautionary case. The finite
discrete-time, probability-normalized implementation can produce short-range
burstiness while its centered one-hot power spectrum remains close to white
noise. Treat MB4 as implemented but not yet a strong symbolic LRD baseline.

## LongMemory.jl Comparison

`validation/longmemory_comparison.jl` checks S5's formalized diagnostic helpers
against LongMemory.jl's `autocovariance`, `autocorrelation`, and `periodogram`.
The comparison documents the adaptations needed for symbolic data and plotting:

- centered one-hot symbol series are used as numeric inputs;
- lag-zero autocorrelation is dropped for S5 plots;
- angular frequencies are converted to cycles per observation;
- zero frequency is dropped for log-log spectral plots.

Run it in the validation environment:

```julia
julia --project=validation -e 'using Pkg; Pkg.instantiate()'
julia --project=validation validation/longmemory_comparison.jl
```

## Large Validation Flags

Long-running studies should be opt-in. Preferred environment variables use the
`S5_` prefix:

- `S5_VALIDATION_LARGE=true`;
- `S5_VALIDATION_REPLICATES=<integer>`;
- `S5_VALIDATION_N=<integer>` or a script-specific size variable.

Retained outputs should record enough provenance to identify generator settings,
sequence length, replicate count, random seed, package version when available,
and creation date.

## Benchmarks

Benchmarks live in `benchmark/` and use their own `Project.toml` with
BenchmarkTools.jl:

```julia
julia --project=benchmark benchmark/benchmarks.jl
```

The default suite covers all implemented generators across moderate sequence
lengths and alphabet sizes. Larger runs are opt-in:

```julia
S5_BENCHMARK_LARGE=true julia --project=benchmark benchmark/benchmarks.jl
```

Rare scaling runs are also opt-in:

```julia
S5_BENCHMARK_SCALING=true julia --project=benchmark benchmark/benchmarks.jl
```

The scaling suite uses `n = 100, 1_000, 10_000, 100_000, 1_000_000`, defers
`k = 64`, and runs `k = 2, 8`. Each BenchmarkTools trial synthesizes 10
independently seeded sequences by default, and retained times are reported as
per-synthesis averages.

Additional knobs:

- `S5_BENCHMARK_SAMPLES=<integer>`;
- `S5_BENCHMARK_SECONDS=<seconds>`;
- `S5_BENCHMARK_SYNTH_REPEATS=<integer>`;
- `S5_BENCHMARK_WRITE_RESULTS=false` to suppress retained result artifacts.

By default the benchmark script writes:

- `benchmark/RESULTS.md`, a machine-specific summary;
- `benchmark/results/benchmarks.csv`, with one row per method, `k`, and `n`;
- histogram-style relative-time SVGs for the largest benchmarked `n`;
- log-log scaling SVGs for each alphabet size.

Benchmark labels include complexity-relevant settings such as `k`, `d`, stream
count for `FSS`, copy-distance settings for `DuplicationMutation`, and
FFT/rank-binning behavior for property-based methods. Interpret benchmark
results as machine- and Julia-version-specific performance evidence, not as
platform-independent speed guarantees.

The retained scaling run shows the expected separation between generator
families: direct sequential methods are fastest, FFT/rank-binning methods scale
well with sequence length, `LGCM` grows with alphabet size, and explicit
history-based methods pay for their configured memory depth.
