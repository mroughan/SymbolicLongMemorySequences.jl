# SymbolicLongMemorySequences.jl Benchmarks

This folder contains performance benchmarks for the implemented SymbolicLongMemorySequences.jl generators.
Benchmarks are machine-specific evidence, not correctness tests.

Run from the package root:

```julia
julia --project=benchmark benchmark/benchmarks.jl
```

The benchmark script writes retained artifacts by default:

- `benchmark/RESULTS.md`: a machine-specific results summary;
- `benchmark/results/benchmarks.csv`: one row per method, alphabet size, and
  sequence length;
- `benchmark/results/relative_times_k*_n*.svg`: histogram-style relative-time
  bar charts at the largest benchmarked sequence length;
- `benchmark/results/scaling_k*.svg`: log-log scaling plots of generation time
  against sequence length.

The default suite covers all implemented generators across moderate sequence
lengths and alphabet sizes. Larger runs are opt-in:

```julia
SLMS_BENCHMARK_LARGE=true julia --project=benchmark benchmark/benchmarks.jl
```

Rare sequence-length scaling studies are opt-in and intentionally heavier:

```julia
SLMS_BENCHMARK_SCALING=true julia --project=benchmark benchmark/benchmarks.jl
```

The scaling suite uses `n = 100, 1_000, 10_000, 100_000, 1_000_000`, defers the
`k = 64` case, and runs `k = 2, 8`. Each BenchmarkTools trial synthesizes 10
independently seeded sequences by default, and retained times are reported as
per-synthesis averages from those trial timings.

Additional knobs:

- `SLMS_BENCHMARK_SAMPLES=<integer>` controls BenchmarkTools sample count.
- `SLMS_BENCHMARK_SECONDS=<seconds>` controls the per-benchmark time budget.
- `SLMS_BENCHMARK_SYNTH_REPEATS=<integer>` controls the number of independently
  seeded syntheses inside each BenchmarkTools trial.
- `SLMS_BENCHMARK_WRITE_RESULTS=false` disables writing `RESULTS.md`, CSV, and
  SVG plot artifacts.

The benchmark labels include complexity-relevant settings such as `k`, `d` for
history-based generators, stream count for `FSS`, copy-distance settings for
`DuplicationMutation`, and FFT/rank-binning behavior for property-based methods.
The factory API itself is intentionally thin and is covered by package tests;
benchmarks measure generator hot paths after construction.

The retained scaling results show a clear split: direct sequential generators
such as `OnOffMarkov`, `FSS`, and `DuplicationMutation` are fastest in the
extended grid, FFT/rank-binning property-based methods scale well with sequence
length, `LGCM` grows with alphabet size, and explicit history methods such as
`LAMP`, `HawkesSymbol`, and `DyadicLAMP` pay for their configured memory depth.
