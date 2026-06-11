# S5.jl Validation Policy

S5.jl separates package correctness from empirical research evidence.

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

### Large Validation Flags

Long-running studies should expose keyword arguments and environment-variable
flags rather than becoming part of the default test path. Preferred flag names use
the `S5_` prefix, for example:

- `S5_VALIDATION_LARGE=true` for larger validation grids;
- `S5_VALIDATION_REPLICATES=<integer>` for replicate counts;
- `S5_VALIDATION_N=<integer>` or script-specific size variables when useful.

If a script writes retained output, it should record enough provenance to identify
the generator settings, sequence length, replicate count, random seed, package
version when available, and creation date.

Diagnostic transformations must be code, not just plotting lore. For LRD visual
diagnostics, the symbolic sequence is transformed into centered one-hot numeric
series before calling autocorrelation, autocovariance, or periodogram routines.
Any adaptation to an external package, such as dropping LongMemory.jl's lag-zero
autocorrelation or converting angular frequency to cycles per observation, should
live in a named helper and be documented beside the validation script.

## Benchmarking

Benchmarks are performance evidence, not correctness tests. They live under
`benchmark/` and use their own `Project.toml` with `BenchmarkTools.jl`.

Run the default benchmark suite from the package root:

```julia
julia --project=benchmark benchmark/benchmarks.jl
```

The default suite should remain moderate. Larger benchmark runs are opt-in:

```julia
S5_BENCHMARK_LARGE=true julia --project=benchmark benchmark/benchmarks.jl
```

Benchmark results should be interpreted as machine- and Julia-version-specific.
Use them to catch performance regressions and document scaling behavior, not to
make platform-independent speed guarantees.

## Future Trigram Validation

S5.jl already provides `empirical_trigram` for diagnostics, but it does not expose
a concrete trigram-control specification. Future trigram-control work should add a
higher-order local-structure specification, focused constructor validation, fast
unit tests for the specification contract, and manual validation studies for
empirical controllability.
