# S5.jl

**Self-Similar Symbols Sequence Synthesis**

S5.jl generates Long-Range Dependent (LRD) sequences of categorical
(non-numerical) symbols. It produces controllable synthetic data for LRD
estimator tests, information-theoretic experiments, and LLM-style neural
sequence models trained on challenging non-language symbolic streams.

This package implements Task 1 of the ARC Discovery Grant project *"Analysis
and Synthesis of Long-Range Structure in Non-Numerical Time Series"*
(Roughan & Willinger, 2023).

## Quick start

```julia
using S5, StableRNGs

rng = StableRNG(42)

# Property-based: spectral fGn + quantization
g1 = SpectralFGN(0.8, [:a, :b, :c], [0.2, 0.3, 0.5])
seq1 = generate(g1, 10_000; rng)

# Property-based: latent Gaussian categorical model
g2 = LGCM(0.8, [:a, :b, :c], [0.2, 0.3, 0.5])
seq2 = generate(g2, 10_000; rng)

# Property-based: multiscale driver + Markov regimes
P1 = [0.9 0.1; 0.2 0.8]
P2 = [0.3 0.7; 0.6 0.4]
g3 = WaveletMarkov(0.8, [:a, :b], [P1, P2])
seq3 = generate(g3, 10_000; rng)

# Model-based: Linear-Additive Markov Process
g4 = LAMP(0.5, [:a, :b, :c], [0.2, 0.3, 0.5]; d = 500, epsilon = 0.02)
seq4 = generate(g4, 10_000; rng)

# Model-based: heavy-tailed regime-switching Markov chain
Q = [0.2 0.8; 0.8 0.2]
g5 = OnOffMarkov(1.5, [:a, :b], [P1, P2], Q)
seq5 = generate(g5, 10_000; rng)

# Model-based: Fractal Symbol Sequence  (H = (3-α)/2 = 0.75)
g6 = FSS(1.5, [:a, :b, :c]; rates = [2.0, 3.0, 5.0])
seq6 = generate(g6, 10_000; rng)

empirical_marginal(seq1, g1.alphabet)

# Save to INC format with full provenance metadata
save_sequence("seq_pb1.inc", seq1, g1)
```

## Generators

| ID  | Type          | LRD mechanism                          | Short-range control        | Complexity      |
|-----|---------------|----------------------------------------|----------------------------|-----------------|
| PB1 | `SpectralFGN` | Spectral $1/f^\alpha$ shaping          | Poor (set by quantization) | $O(n \log n)$  |
| PB2 | `LGCM`        | Latent fGn streams + argmax            | Offset-calibrated marginals | $O(n \cdot k \cdot I)$ |
| PB3 | `WaveletMarkov` | Multiscale driver + Markov regimes   | Per-regime Markov matrices | $O(n \log n + n \cdot k)$ |
| MB1 | `LAMP`        | Power-law history weights              | Weight tensor              | $O(n \cdot d)$ |
| MB2 | `OnOffMarkov` | Heavy-tailed regime sojourns           | Per-regime Markov matrices | $O(n \cdot k)$ |
| MB3 | `FSS`         | Pareto renewal process per symbol      | Poor (independent streams) | $O(n \cdot k)$ |

### LRD parameters

Each generator exposes either a direct Hurst parameter $H \in (1/2, 1)$ or a
tail/decay parameter with a nominal relationship to $H$:

| Type          | Input parameter  | Relationship to $H$    |
|---------------|------------------|------------------------|
| `SpectralFGN` | `H`              | direct                 |
| `LGCM`        | `H`              | direct                 |
| `WaveletMarkov` | `H`            | latent-driver target   |
| `LAMP`        | `beta` ($\beta$) | $H = (2 - \beta) / 2$ |
| `OnOffMarkov` | `alpha` ($\alpha$) | nominal $H = (3 - \alpha) / 2$ |
| `FSS`         | `alpha` ($\alpha$) | $H = (3 - \alpha) / 2$ |

## Controllability

All implemented generators accept an explicit ordered alphabet. Duplicate alphabet
entries are rejected because they make empirical frequency tables ambiguous.

Use `control_capabilities(g)` to inspect the strength and scope of the controls a
generator claims. It distinguishes direct finite-sample or empirical controls from
implied, asymptotic, induced, latent, and nominal behavior.

| Type | Marginal control | Local structure control |
|------|------------------|-------------------------|
| `SpectralFGN` | direct `marginal`; rank binning gives integer counts as close as possible to target | none |
| `LGCM` | direct `marginal`; calibrated latent offsets | none |
| `WaveletMarkov` | aggregate stationary marginal implied by regimes | per-regime bigram matrices |
| `LAMP` | direct `marginal` mixed through `epsilon`; larger `epsilon` improves marginal control | history-weighted dependence, not arbitrary bigrams |
| `OnOffMarkov` | aggregate stationary marginal implied by regimes | per-regime bigram matrices |
| `FSS` | asymptotic `rates / sum(rates)` | none |

First-order local structure can be represented by a validated `MarkovSpec`:

```julia
spec = MarkovSpec([:a, :b], [0.9 0.1; 0.2 0.8])
g1 = WaveletMarkov(0.8, [spec, spec])
g2 = OnOffMarkov(1.5, [spec, spec], [0.2 0.8; 0.8 0.2])

control_capabilities(g1).bigram
```

`WaveletMarkov` and `OnOffMarkov` provide `:per_regime` bigram control. Their
aggregate observed bigrams depend on regime differences and switching behavior. If
all regimes use the same `MarkovSpec`, that common transition matrix is also the
unambiguous aggregate target.

`MarkovSpec` is the current concrete `LocalStructureSpec`, and
`local_structure_order(spec)` returns its order. This is the intended extension
point for future sparse higher-order controls. S5.jl currently provides
`empirical_trigram` for diagnostics, but it does not expose a trigram-control
specification.

For `WaveletMarkov` and `OnOffMarkov`, one-hot symbol diagnostics need regimes
with different stationary symbol distributions. If each regime has the same
stationary marginal, the latent regime process can be long-memory while the
symbol-level ACF and spectrum look nearly short-memory.

Reproducible controllability studies live in `validation/`, for example:

```julia
julia --project=. validation/marginal_control.jl
julia --project=. validation/local_structure.jl
julia --project=. validation/lrd_method_diagnostics.jl
julia --project=validation validation/longmemory_comparison.jl
```

The LRD diagnostic script writes generated sequences and summary tables as INC
files under `validation/results/lrd_diagnostics/`, and writes log-log SVG plots of
log-binned one-hot autocorrelation and power-spectrum summaries.
The symbolic-to-numeric transformation is formalized in
`validation/lrd_symbol_diagnostics.jl`: each symbol is converted to a centered
one-hot indicator series before autocorrelation, autocovariance, or periodogram
calculations. `validation/longmemory_comparison.jl` compares those helpers with
LongMemory.jl's `autocovariance`, `autocorrelation`, and `periodogram`, including
the documented lag-zero and angular-frequency adaptations.

See `VALIDATION_POLICY.md` for the validation tiers. Fast tests run through the
package test suite. Longer validation studies and larger benchmark runs are manual
or flag-controlled.

## Benchmarking

Benchmarking uses a separate environment under `benchmark/`:

```julia
julia --project=benchmark benchmark/benchmarks.jl
S5_BENCHMARK_LARGE=true julia --project=benchmark benchmark/benchmarks.jl
```

The default suite covers all implemented generators across moderate sequence
lengths and alphabet sizes. The large suite adds longer runs and should be treated
as machine-specific performance evidence rather than a correctness test.

## Motivations

S5.jl is intended for several related research uses:

- testing NN-LRD estimators on sequences with known generator settings;
- probing excess entropy rate, entropy-rate convergence, and other
  information-theoretic summaries of long-memory symbolic data;
- creating non-language but language-like long-context sequences for training and
  stress-testing LLM-style neural networks;
- separating local bigram/trigram competence from genuine long-context modelling;
- studying anomaly detection, change detection, and synthetic privacy-preserving
  categorical data with realistic burstiness.

## INC output format

All generated sequences should be saved with `save_sequence`, which writes the
data and full provenance metadata (package version, generator type and
parameters, creation date) in the
[INC format](https://github.com/mroughan/IncCSV.jl).

```julia
save_sequence("output.inc", seq, g)
```

The resulting file has this structure:

```
---
title = S5.jl synthetic LRD symbol sequence
package = S5
version = "0.1.0"
generator = SpectralFGN
method = PB1
n = 10000
created = "2026-06-04"
[generator_params]
H = "0.8"
alphabet_size = 3
alphabet = "a,b,c"
marginal = "0.33333333,0.33333333,0.33333334"
[columns]
index = time index (1-based)
symbol = generated symbol from the alphabet
---
index,symbol
1,b
2,a
…
```

## Background

Long-Range Dependence (LRD) means that the large-scale statistical structure
of a sequence is as important as its short-range structure. Formally, a
sequence is LRD if its autocovariance function (ACF) decays as a power law:

$$\gamma_k \sim c_\gamma |k|^{-\beta}, \quad \beta \in (0, 1),$$

with Hurst parameter $H = (2-\beta)/2 \in (1/2, 1)$.

LRD is ubiquitous in human-generated data (natural language, Internet traffic,
genomics, social media), yet almost all LRD synthesis tools target numerical
data. S5.jl fills that gap for **symbol sequences** — data taking values in a
finite, unordered alphabet.

## References

- Roughan, M. & Willinger, W. (2023). ARC Discovery Grant proposal.
- Paxson, V. (1997). Fast, approximate synthesis of fractional Gaussian noise. *CCR* 27.
- Kumar, R., et al. (2017). Linear additive Markov processes. *WWW '17*.
- Lowen, S. B. & Teich, M. C. (1995). Fractal stochastic point processes. *Fractals* 3(1).
- Pipiras, V. & Taqqu, M. S. (2017). *Long-Range Dependence and Self-Similarity*. Cambridge UP.
