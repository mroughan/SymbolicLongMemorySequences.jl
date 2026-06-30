# SymbolicLongMemorySequences.jl Validation Studies

This folder contains reproducible simulation studies for generator controllability.
These are not LRD estimators. They test whether a generator respects user-facing
controls such as alphabet membership and target marginal probabilities.
Standard validation factories should use the explicit scientific constructors
when a study needs exact parameter visibility; `make_generator` is available for
smoke tests and examples where standard cases are enough.

The scripts use `StableRNGs.StableRNG` so results are reproducible across Julia
sessions and package updates.

See `../VALIDATION_POLICY.md` for the project validation tiers. Fast contract tests
belong in `test/`; broader empirical evidence belongs here and should be run
manually or behind explicit `SLMS_` environment flags.

## Marginal Control

Run from the package root:

```julia
julia --project=. validation/marginal_control.jl
```

The script prints aggregate total-variation and maximum absolute marginal errors for
all implemented generators across a small grid of sequence lengths, alphabet sizes,
and marginal distributions.

The same script also writes a focused uniform-marginal validation case for
paper and report use:

- `validation/results/marginal_control/uniform_marginal_k8_summary.csv`;
- `validation/results/marginal_control/uniform_marginal_k8_histogram_data.csv`;
- `validation/results/marginal_control/uniform_marginal_histograms_k8.svg`;
- `validation/results/marginal_control/uniform_marginal_histograms_k8.pdf`
  when a local converter such as `rsvg-convert` or Inkscape is available.

This case uses `k = 8`, `n = 100_000`, 20 replicates, and a uniform categorical
target. The first and last 10% of each generated sequence are dropped before
frequencies are computed, reducing sensitivity to initialization and finite-end
effects. The CSV reports chi-squared frequency diagnostics against the intended
uniform marginal. Because the generated sequences are dependent, the iid
multinomial chi-squared p-values are only a reference diagnostic, not exact
hypothesis-test p-values. LRD can make empirical marginals converge more slowly
than the iid multinomial model predicts, so this reference is often too
conservative while still being informative.
The summary also includes an approximate effective-sample-size correction. For
each replicate, centered one-hot indicators estimate an integrated
autocorrelation time for every symbol; the smallest symbol ESS is used as a
conservative `effective_n`, and the chi-squared statistic is scaled by
`effective_n / trimmed_n`. The adjusted p-values remain diagnostics, not exact
LRD categorical tests, but they give a more realistic frequency check when
dependence is strong.
The CSV also reports full-sequence total-variation and maximum absolute errors
beside the trimmed-window errors. This distinction matters for rank-binned and
empirically calibrated property-based methods: `SpectralFGN`,
`IntermittentMapSymbols`, and usually `LGCM` can match the full generated sample
very closely while still showing interior-window deviations after trimming.

For MB1a (`LAMP`) and MB1b (`DyadicLAMP`), this marginal-control study uses
`lamp_repeat_transition(p; repeat_probability = 0.4)`. This keeps a
repeat-biased LAMP mechanism while making the transition matrix ergodic with
the requested marginal as its stationary distribution. A pure identity
transition is a useful stress case, but it can lock in early finite-sample
imbalances and is not the standard marginal-control validation case.

A more formal marginal test should calibrate the null distribution under
dependence. Good candidates are block/subsampling tests or a parametric Monte
Carlo envelope generated from the same configured generator. Those are better
suited to manual validation reports than to the fast package test suite.

`LGCM` is more expensive than the other methods because it calibrates latent
offsets over an `n × k` matrix. Increase `replicates`, `ns`, or
`calibration_iters` manually when running longer studies.

Keep generated data out of the repository unless a result table is intentionally being
tracked.

Large validation grids should be opt-in. Prefer keyword arguments in the script API
and environment variables such as `SLMS_VALIDATION_LARGE=true` or
`SLMS_VALIDATION_REPLICATES=<integer>` when a script grows beyond the default manual
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

For larger alphabets, directly testing every entry of a transition matrix becomes
data-hungry because there are `k^2` cells and some rows may be visited rarely.
Shorter one-step diagnostics should test interpretable contrasts first:
stationary-weighted row total variation, repeat probability, selected important
rows, grouped symbol transitions, or the action of the matrix on a small number
of contrast vectors. Full row-by-row tests are still appropriate for small `k`
or explicitly large validation runs.

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
PB3 is plotted as separate Haar and spectral-driver variants so the latent-driver
choice can be compared directly. MB2 uses `L_min = 50.0` so the heavy-tailed
sojourn mechanism is visible at `n = 100_000`. If all regimes share the same
stationary marginal, these diagnostics can look short-memory even when the
latent regime process has long-range structure.
MB4 (`HawkesSymbol`) uses identity excitation so recent symbols increase their
own future intensity through the configured finite power-law history kernel.
MB1c (`CalibratedAdditiveMarkov`) uses a centered additive memory function.
MB5 (`DuplicationMutation`) uses copy-mutate growth with a truncated power-law
copy-distance kernel. PB4 (`IntermittentMapSymbols`) uses an
intermittent latent driver followed by rank binning.

The PB3 split is diagnostic evidence rather than a completed calibration claim.
The legacy Haar cascade tends to retain a high, flat autocorrelation shelf, while
the spectral driver can be damped strongly by the rank-binning and Markov
emission layers. Further PB3 work should validate the latent driver, regime
indicators, and emitted symbols separately.

The current MB4 validation should also be read cautiously. `HawkesSymbol`
produces short-range burstiness, but the centered one-hot power spectrum can look
close to white noise under the default finite discrete-time construction. A
better MB4 path likely needs a near-critical or event-count formulation rather
than only increasing the identity-excitation strength in the present
probability-normalized sampler.

The script computes the one-hot symbol autocorrelation and power spectrum for each
sequence, averages across symbols and replicates, and writes:

- `average_autocorrelation.inc`: full signed averaged autocorrelation by lag;
- `average_power_spectrum.inc`: full averaged periodogram by Fourier frequency;
- `plot_autocorrelation_logbins.inc`: positive log-binned values used for plotting;
- `plot_power_spectrum_logbins.inc`: log-binned spectrum values used for plotting.

For property-based methods, the script calls `generate_with_latent` and also
writes:

- `latent_average_autocorrelation.inc`: full latent numerical autocorrelation;
- `latent_average_power_spectrum.inc`: full latent numerical periodogram;
- `latent_plot_autocorrelation_logbins.inc`: log-binned latent ACF values;
- `latent_plot_power_spectrum_logbins.inc`: log-binned latent spectrum values.

The log-binning code keeps the final bin closed only at its upper edge and sorts
the binned x-values before writing tables and SVG polylines, so plot x-values are
strictly increasing within each method. Paired log-log SVG plots are written
under `validation/results/lrd_diagnostics/plots/`, with autocorrelation on the
left and power spectrum on the right. Property-based plots have two rows: the
latent numerical process on top and the final symbolic one-hot diagnostics on
the bottom. If `rsvg-convert` is installed, matching PDF files are produced for
paper inclusion.

The SVG plots include vertical dashed interpretation limits. The red line marks
the finite-sample lag limit `n / 10`, chosen so autocorrelation estimates are not
interpreted deep into the range where too few overlapping pairs remain. On power
spectrum plots this same scale is shown as frequency `10 / n`. Methods with an
explicit internal memory limit may add a second dashed line; for example, `LAMP`,
`DyadicLAMP`, `CalibratedAdditiveMarkov`, and `HawkesSymbol` mark their
configured history depth `d`, while `DuplicationMutation` marks its configured
maximum copy-distance window.
Where a generator has a defensible asymptotic-onset scale, the plots also mark an
approximate power-law onset. For `OnOffMarkov` this is the Pareto scale `L_min`.
For `HawkesSymbol`, whose kernel is `(lag + c)^(-beta)`, the onset is the lag at
which the kernel's local log-log slope reaches 90% of the asymptotic slope:
`ceil(0.9c / 0.1)`. This is only a visual guide; the transition is gradual and
does not imply that earlier lags are irrelevant.
Autocorrelation and power-spectrum plots also include gray dashed nominal
power-law reference lines, anchored at the first positive plotted value. The ACF
reference has slope `lag^(-beta)`, while the spectral-density reference has the
corresponding low-frequency slope `frequency^(beta - 1)`. The references use each
generator's nominal decay exponent and are visual guidance, not fitted curves.

These limits are visual guidance only. They do not prove LRD behavior, and they
make visible when poor diagnostics are caused by a generator's own truncation
rather than by estimator plotting alone.

## LongMemory.jl Comparison

SymbolicLongMemorySequences.jl does not depend on estimator packages at runtime, but validation can compare
the formalized diagnostic transformations with
[`LongMemory.jl`](https://github.com/everval/LongMemory.jl). Instantiate the
validation environment, then run:

```julia
julia --project=validation -e 'using Pkg; Pkg.instantiate()'
julia --project=validation validation/longmemory_comparison.jl
```

The comparison script adapts LongMemory.jl conventions explicitly:

- `autocovariance(x, k)` and `autocorrelation(x, k)` return lags `0:k-1`;
- SymbolicLongMemorySequences plots use lags `1:maxlag`, so the lag-zero autocorrelation is dropped;
- LongMemory.jl periodograms report angular frequencies and include zero;
- SymbolicLongMemorySequences plots use cycles per observation and drop zero frequency.

The script reports maximum absolute differences between SymbolicLongMemorySequences's local
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
SLMS_BENCHMARK_LARGE=true julia --project=benchmark benchmark/benchmarks.jl
```

Run the rare sequence-length scaling suite with:

```julia
SLMS_BENCHMARK_SCALING=true julia --project=benchmark benchmark/benchmarks.jl
```

The benchmark script writes `../benchmark/RESULTS.md`, CSV results, and SVG
plots by default. The relative-time plots are histogram-style bars normalized to
the fastest method for each alphabet size at the largest benchmarked sequence
length. The scaling plots show generation time against sequence length on
log-log axes. The retained scaling run uses 10 independently seeded syntheses
inside each BenchmarkTools trial and reports per-synthesis average times.
