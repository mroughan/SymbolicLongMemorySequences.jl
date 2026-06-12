# S5.jl Benchmarks

This folder contains performance benchmarks for the implemented S5.jl generators.
Benchmarks are machine-specific evidence, not correctness tests.

Run from the package root:

```julia
julia --project=benchmark benchmark/benchmarks.jl
```

The default suite covers all implemented generators across moderate sequence
lengths and alphabet sizes. Larger runs are opt-in:

```julia
S5_BENCHMARK_LARGE=true julia --project=benchmark benchmark/benchmarks.jl
```

Additional knobs:

- `S5_BENCHMARK_SAMPLES=<integer>` controls BenchmarkTools sample count.
- `S5_BENCHMARK_SECONDS=<seconds>` controls the per-benchmark time budget.

The benchmark labels include complexity-relevant settings such as `k`, `d` for
`LAMP` and `HawkesSymbol`, stream count for `FSS`, and FFT length behavior for
`SpectralFGN`.
