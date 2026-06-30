# SymbolicLongMemorySequences.jl Validation Policy

SymbolicLongMemorySequences.jl separates package correctness from empirical research evidence.

Fast tests protect public contracts: constructor validation, reproducible RNG use,
output length and element type, alphabet membership, provenance metadata, and small
bounded checks for declared controls. These tests run through:

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```

Longer validation studies live under `validation/`. They may generate larger
sequences, run more replicates, write result tables, or create plots. These studies
support empirical claims, but they are not substitutes for unit tests and are not
LRD estimators.

## Validation Tiers

### Fast Path

The main development path is the package test suite. Tests should finish quickly
enough for normal development and CI. Use exact assertions for API contracts and
small, seeded statistical checks only when the behavior cannot be tested exactly.

### Manual Validation

Manual validation scripts are reproducible Julia programs in `validation/`.
Run them from the package root:

```julia
julia --project=. validation/marginal_control.jl
julia --project=. validation/local_structure.jl
julia --project=. validation/lrd_method_diagnostics.jl
julia --project=validation validation/longmemory_comparison.jl
```

Scripts should use `StableRNGs.StableRNG`, state their grid of sample sizes and
replicates, and write only aggregate outputs unless generated sequences are needed
for a documented diagnostic. Large generated sequences should stay out of version
control.

Marginal validation should include at least one simple, interpretable target.
The default paper-facing case is a uniform categorical distribution with
`k = 8`, dropping the first and last 10% of each sequence before computing
frequencies. The retained outputs should include a histogram plot and a table of
chi-squared frequency diagnostics. The raw chi-squared reference distribution is
an iid multinomial approximation; for dependent LRD sequences it should be
described as a diagnostic, not as an exact hypothesis test. Because LRD can slow
convergence of empirical averages, the iid reference is usually too strict but
still useful as an indicator of possible marginal-control problems. Retained
uniform-marginal outputs should also include the implemented
effective-sample-size correction, which estimates the integrated autocorrelation
time of centered one-hot symbol indicators and scales the chi-squared statistic
by the conservative ratio `effective_n / trimmed_n`.

When a validation claim depends on formal marginal testing, prefer a
dependence-aware calibration instead of a naive chi-squared p-value. The ESS
correction is an approximate diagnostic; stronger extensions include block or
subsampling tests, or a parametric Monte Carlo envelope generated from the same
configured generator. These belong in manual validation, not the fast test path.

Validation scripts may use `make_generator` for smoke tests of standard cases,
but studies that support method-specific claims should prefer explicit
constructors so the scientific parameters and assumptions remain visible in the
script.

### Large Validation Flags

Long-running studies should expose keyword arguments and environment-variable
flags rather than becoming part of the default test path. Preferred flag names use
the `SLMS_` prefix, for example:

- `SLMS_VALIDATION_LARGE=true` for larger validation grids;
- `SLMS_VALIDATION_REPLICATES=<integer>` for replicate counts;
- `SLMS_VALIDATION_N=<integer>` or script-specific size variables when useful.

If a script writes retained output, it should record enough provenance to identify
the generator settings, sequence length, replicate count, random seed, package
version when available, and creation date.

One-step Markov validation should start with aggregate contrasts when `k` is
large. Full matrix tests have `k^2` cells and can require much longer sequences
than marginal checks. Useful shorter diagnostics include repeat probability,
stationary-weighted row total variation, selected row contrasts, and grouped
symbol transitions. Full row-by-row tests should be reserved for smaller
alphabets or explicitly large validation runs.

Diagnostic transformations must be code, not just plotting lore. For LRD visual
diagnostics, the symbolic sequence is transformed into centered one-hot numeric
series before calling autocorrelation, autocovariance, or periodogram routines.
Any adaptation to an external package, such as dropping LongMemory.jl's lag-zero
autocorrelation or converting angular frequency to cycles per observation, should
live in a named helper and be documented beside the validation script.

Visual diagnostics should mark interpretation limits. For autocorrelation plots,
the default finite-sample limit is lag `n / 10`, because estimates at much larger
lags rely on a shrinking number of overlapping pairs. Spectrum plots should mark
the reciprocal scale. If a generator has an explicit internal cutoff, such as
`LAMP.d`, the plot should mark that as a separate generator limit.
When a generator has a defensible short-range-to-tail onset scale, plots may also
mark an approximate power-law onset. This is a visual guide, not a fitted
transition point. For example, a Pareto sojourn model can use its scale
parameter, while an offset kernel `(lag + c)^(-beta)` can mark the lag where its
local log-log slope reaches a chosen fraction of the asymptotic slope.

## Benchmarking

Benchmarks are performance evidence, not correctness tests. They live under
`benchmark/` and use their own `Project.toml` with `BenchmarkTools.jl`.
Retained benchmark results should include a machine-specific summary, a
machine-readable table, and plots that make relative timing and scaling with
sequence length visible. Rare scaling studies should be opt-in, may use much
larger `n`, and should average each reported timing over multiple independently
seeded syntheses rather than relying on a single generated sequence.

Run the default benchmark suite from the package root:

```julia
julia --project=benchmark benchmark/benchmarks.jl
```

The default suite should remain moderate. Larger benchmark runs are opt-in:

```julia
SLMS_BENCHMARK_LARGE=true julia --project=benchmark benchmark/benchmarks.jl
SLMS_BENCHMARK_SCALING=true julia --project=benchmark benchmark/benchmarks.jl
```

Benchmark results should be interpreted as machine- and Julia-version-specific.
Use them to catch performance regressions and document scaling behavior, not to
make platform-independent speed guarantees.

## Future Trigram Validation

SymbolicLongMemorySequences.jl already provides `empirical_trigram` for diagnostics, but it does not expose
a concrete trigram-control specification. Future trigram-control work should add a
higher-order local-structure specification, focused constructor validation, fast
unit tests for the specification contract, and manual validation studies for
empirical controllability.
