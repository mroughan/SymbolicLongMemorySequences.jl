# LongMemory.jl API Notes

Downloaded/inspected on: 2026-06-11

Source repository: <https://github.com/everval/LongMemory.jl>

Commit inspected: `4f95f98f05e243096bf7b57a92cae687d75252f7`

LongMemory.jl version observed in `Project.toml`: `1.0.1`

Relevant source files:

- `src/ClassicEstimators.jl`
- `src/LogPeriodEstimators.jl`
- `test/runtests.jl`

Observed API conventions used by SymbolicLongMemorySequences validation:

- `autocovariance(x::Array, k::Int)` returns lags `0:k-1` as a `k × 1` array.
- `autocovariance` subtracts the sample mean and divides each lag sum by
  `T = length(x)`, not by `T - lag`.
- `autocorrelation(x::Array, k::Int; flag = false)` normalizes by the lag-zero
  autocovariance and returns lags `0:k-1`.
- `periodogram(x::Array)` returns `(I_w, w)`, where `w` is angular frequency and
  includes zero frequency.

SymbolicLongMemorySequences validation adapts these conventions by converting symbolic sequences to
centered one-hot numeric series, dropping the lag-zero autocorrelation for plots
over lags `1:maxlag`, converting angular frequency to cycles per observation,
and dropping zero frequency before log-log spectrum plots.
