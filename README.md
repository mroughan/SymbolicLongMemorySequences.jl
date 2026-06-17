# S5.jl

[![Package tests](https://github.com/mroughan/S5.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/mroughan/S5.jl/actions/workflows/ci.yml)
[![Aqua](https://github.com/mroughan/S5.jl/actions/workflows/aqua.yml/badge.svg)](https://github.com/mroughan/S5.jl/actions/workflows/aqua.yml)
[![JET](https://github.com/mroughan/S5.jl/actions/workflows/jet.yml/badge.svg)](https://github.com/mroughan/S5.jl/actions/workflows/jet.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://mroughan.github.io/S5.jl/dev/)
[![Codecov](https://codecov.io/gh/mroughan/S5.jl/branch/main/graph/badge.svg?token=WYXfD9ij0s)](https://codecov.io/gh/mroughan/S5.jl)
[![Documenter](https://github.com/mroughan/S5.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/mroughan/S5.jl/actions/workflows/documentation.yml)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mroughan.github.io/S5.jl/dev)
[![Julia](https://img.shields.io/badge/julia-1.10%2B-blue.svg)](https://julialang.org/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Self-Similar Symbols Sequence Synthesis**

Long-memory is a feature of many natural sequences, closely relate to statistical self-similarity. In the past, numerical sequences starting with the size of Nile river floods have been the main topic of study. But there are many domains where there are (hypothesised) long-memory sequences that are non-numerical, symbolic sequences. For instance: human writing, DNA, and so forth.

S5.jl is a Julia package for generating Long-Range Dependent (LRD) sequences of
non-numerical (categorical/symbolic) data. It provides controllable synthetic
sequences for testing LRD estimators, probing information-theoretic quantities,
and training or stress-testing sequence models on non-language data with
language-like long context.

---

## Motivations

S5.jl is useful beyond basic generator benchmarking:

- **Estimator testing:** generate labelled NN-LRD (NonNumerical LRD) sequences where alphabet, marginal
  distribution, local structure, and nominal LRD mechanism are known.

- **Information-theoretic experiments:** create controlled cases for ideas such as
  excess entropy rate, entropy-rate convergence, and the gap between local and
  long-range predictability in LRD processes.

- **Non-language sequence modelling:** train or stress-test LLM-style neural sequence
  models on symbolic data that is not text, such as event logs, vulnerability classes,
  genomic symbols, workflow traces, or synthetic protocol states.

- **Context-length diagnostics:** test whether models exploit genuinely long context
  rather than only short-range bigram/trigram cues.

- **Anomaly and change-detection studies:** create controlled shifts in marginal,
  local Markov structure, regime persistence, or long-range behaviour.

- **Privacy-preserving simulation:** produce synthetic categorical sequences with
  realistic burstiness without copying a sensitive corpus.

These applications require the abililty to generate 

+ very long sequences (million or billions of tokens seems a reasonable starting point); with

+ control over short-term behaviour (marginals) as well as enforcing long-memory. 
 
---

## Background

Long-Range Dependence (LRD) means that the large-scale statistical structure of a
sequence is as important as its short-range structure. Formally, a sequence is LRD if
its autocovariance function (ACF) decays as a power law,

$$\gamma_k \sim c_\gamma |k|^{-\beta}, \quad \beta \in (0, 1),$$

so that the sum of the tail ACF diverges. LRD is characterised by closely related
parameters: the ACF decay exponent $\beta$, the spectral exponent $\alpha = 1 - \beta$,
and the Hurst parameter $H = (2 - \beta)/2$ with $H \in (1/2, 1)$.

LRD is ubiquitous in human-generated data (text, Internet traffic, genomics, social
media), yet almost all synthesis tools target numerical data. S5.jl fills that gap for
**symbol sequences** — data that takes values in a finite, unordered alphabet such as
`{Orange, Apple,  Pear, ...}` or `{G, A, C, T}`.

---

## Implemented Methods

Methods are broadly classified into **property-based** and **model-based**. The former largely aim to synthesize LRD by starting with a numerical sequence with known LRD properties and then crafting a symbolic sequence by transforming the numerical data. The latter start from a model that has properties such as hierarchical structure or power-law distributed times to drive the sequence generation directly.

Complexity notation used below:

- `n`: generated sequence length;
- `d`: configured history depth or effective memory cutoff;
- `k`: alphabet size;
- `I`: number of calibration iterations.

Standard cases can be constructed through the uniform factory API:

```julia
using S5, StableRNGs

alphabet = [:a, :b, :c]
g = make_generator(:PB1, alphabet; H = 0.8, marginal = [0.2, 0.3, 0.5])
seq = generate(g, 10_000; rng = StableRNG(42))

method_ids()
method_info(:PB1).defaults
```

`make_generator(id, alphabet; kwargs...)` accepts IDs such as `:PB1`, `"MB1c"`,
or type-name aliases such as `:SpectralFGN`. It is intended for common starting
points; the explicit constructors below remain the full-control API.

### Property-Based Methods

These generate one or more underlying numerical LRD processes and then map them
to symbols. In code, this can now be expressed directly as a composition of a
`LatentSource` and a `Symbolizer`:

```julia
source = SpectralFGNSource(0.8)
symbolizer = QuantileSymbolizer([:a, :b, :c], [0.2, 0.3, 0.5])
g = PropertyBasedGenerator(source, symbolizer)
seq = generate(g, 10_000; rng = StableRNG(42))
```

The named PB methods remain the standard, documented cases. The composable API
makes the latent-source/symbolization split explicit, while still checking
compatibility at construction time. For example, a quantile symbolizer needs one
latent series, an argmax symbolizer needs one latent series per symbol, and an
intermittent-map source is currently single-stream only.

| Named method | Latent source | Symbolizer |
|--------------|---------------|------------|
| PB1 | `SpectralFGNSource` | `QuantileSymbolizer` |
| PB2 | `SpectralFGNSource` | `ArgmaxSymbolizer` |
| PB3 | `SpectralFGNSource` or `HaarLRDSource` | `MarkovRegimeSymbolizer` |
| PB4 | `IntermittentMapSource` | `QuantileSymbolizer` |

The LRD property is inherited from the numerical layer, then altered by the
symbolization step. Validation therefore reports the behavior of the full
composition, not just the latent process.

| ID | Name | LRD mechanism | Short-range control | Complexity | Novel? |
|----|------|---------------|---------------------|------------|--------|
| PB1 | Spectral fGn + quantization | Spectral $1/f^\alpha$ shaping | Poor (set by quantization) | $O(n \log n)$ | No |
| PB2 | Latent Gaussian categorical (LGCM) | fGn streams + argmax | Via calibrated offsets | $O(n k I)$ | No |
| PB3 | Spectral/wavelet driver + Markov state machine | Latent LRD driver rank-binned into regimes | Markov transition matrices | $O(n \log n + n k)$ | Partial |
| PB4 | Intermittent map + quantization | Latent intermittent dynamics | Poor (set by quantization) | $O(n \log n)$ | No |

**PB1 — Spectral fGn + quantization.**
Fractional Gaussian noise with Hurst parameter $H$ is synthesized using the fast,
approximate spectral method of Paxson (1997). The real-valued output is sorted into
$k$ rank bins; each bin maps to one symbol, with integer bin counts chosen to match a
target marginal distribution as closely as possible for the finite sample. This is the
simplest approach and serves as the primary validation baseline.

**PB2 — Latent Gaussian categorical model (LGCM).**
A vector of $k$ latent fGn streams is generated, one stream per symbol. At each time
step the symbol is the argmax of the latent vector plus calibrated per-symbol offsets.
The offsets shift marginal probabilities while the latent streams carry the LRD
structure. This is a practical finite-sample approximation to the latent Gaussian
categorical model of Gal, Chen & Ghahramani (ICML 2015).

**PB3 — Latent LRD driver with a Markov state machine.**
A latent long-memory driver controls which Markov regime is active at each step.
Each regime has its own symbol transition matrix, so local bigram structure can be
prescribed while the latent driver injects persistence across scales. The default
`driver = :spectral` uses approximate spectral fGn synthesis followed by
rank-binning into regimes; `driver = :haar` retains the original Haar-like cascade
as a comparison path. This is a practical implementation of the wavelet/state-machine
idea in Roughan, Veitch & Abry (2000), with a fully calibrated wavelet variant left
as a research extension.

**PB4 — Intermittent map + quantization.**
A Pomeau-Manneville-style intermittent map generates a latent real-valued driver
with long laminar episodes. The driver is rank-binned to symbols, so finite-sample
marginal counts are controlled in the same spirit as PB1. This keeps the method in
the property-based family: long-range structure lives in a latent dynamical system,
not in an explicit symbolic transition rule.

---

### Model-Based Methods

These produce LRD through the stochastic model itself rather than via mapping.

| ID | Name | LRD mechanism | Short-range control | Complexity | Novel? |
|----|------|---------------|---------------------|------------|--------|
| MB1a | Linear-Additive Markov Process (LAMP) | Exact power-law history weights | History-weighted transition matrix | $O(n \cdot \min(d,n))$ | No |
| MB1b | Dyadic-bucket LAMP | Dyadic approximation to power-law history | History-weighted transition matrix | $O(n k \log n \log \min(d,n))$ | Partial |
| MB1c | Calibrated additive Markov chain | Centered power-law memory function | Symbol recurrence through additive memory | $O(n \cdot \min(d,n))$ | No |
| MB2 | Heavy-tailed On/Off doubly-stochastic Markov chain | Pareto regime sojourn times | Per-regime Markov chains | $O(n \cdot k)$ | No |
| MB3 | Fractal Symbol Sequence (FSS) via FRP/FSNP | Fractal point process inter-arrivals | Poor (independent streams) | $O(n \cdot k)$ | **Yes** |
| MB4 | Hawkes-style symbolic process | Power-law self/cross-excitation over history | Excitation matrix | $O(n \cdot k \cdot \min(d,n))$ | No |
| MB5 | Duplication-mutation growth | Power-law lag copy/mutate growth | Poor (copy structure induced by growth) | $O(n \log d + d)$ | No |

**MB1a/MB1b/MB1c — Additive history models.**
Transition probabilities are a weighted sum over transition-matrix rows selected by
the observed history,

$$q(s) = (1-\epsilon)\sum_{k=1}^{d} w_k \cdot P[X_{t-k}, s] + \epsilon p(s),$$

with weights $w_k \propto k^{-(1+\beta)}$ targeting a power-law decay up to the
finite observed history range. If `d` exceeds the sequence length, only observed
history contributes and the missing pre-history mass is assigned to the target
marginal. The default transition matrix is identity, recovering copy-from-history
behavior. A simple repeat-biased choice is
`lamp_repeat_transition(marginal; repeat_probability)`, an identity/dyad mixture
whose dyad rows equal the requested marginal. The small innovation term
$\epsilon p(s)$ keeps the requested marginal distribution active after
initialization and prevents finite-history absorption.
The Custom Decay Language Model (CDLM) of Singh, Greenberg & Klakow (2016) is a close
variant demonstrated on text. For large alphabets the weight tensor may be compressed
via low-rank approximations.

`LAMP` is now treated as **MB1a**, the exact finite-history implementation. It is
useful for testing and moderate sequence lengths, but becomes expensive when
`d >= n`. `DyadicLAMP` is **MB1b**, a scalable approximation that groups history
lags into dyadic age buckets such as `1`, `2:3`, `4:7`, and so on. Each bucket
contributes its total power-law weight times the empirical symbol mix in that age
range. MB1b keeps the same transition-matrix controls while making much larger
effective memory depths feasible.

`CalibratedAdditiveMarkov` is **MB1c**, a centered additive Markov-chain memory
function:

$$q(s) = p(s) + \rho \sum_{k=1}^{d} w_k \left(1[X_{t-k}=s] - p(s)\right),$$

with $w_k \propto k^{-\beta}$ and $\rho \in [0,1]$. It is closer to the additive
Markov-chain memory-function literature than LAMP's transition-row mixture, and
is the clearest path toward future correlation-calibrated symbolic generators.

**MB2 — Heavy-tailed On/Off doubly-stochastic Markov chain.**
The sequence is generated by a Markov chain that alternates between two or more regimes.
Sojourn times in each regime follow a Pareto distribution with tail index
$\alpha \in (1, 2)$, so the variance of the symbol count function grows super-linearly,
with nominal $H = (3 - \alpha)/2$. Within each regime a standard (SRD) Markov chain
governs symbol emissions, giving direct control of local statistics. Analogous to
Fractal Shot Noise Processes adapted to symbol sequences (Ryu & Lowen 1998; Garrett
& Willinger 1994).

**MB3 — Fractal Symbol Sequence (FSS) via FRP/FSNP.**
Each symbol $s_i$ is assigned an independent Fractal Renewal Process (FRP) or Fractal
Shot Noise Process (FSNP) governing the times at which that symbol is emitted. The
final sequence merges all symbol streams, with the earliest pending event at each step
determining the output. LRD arises in each symbol's count process through heavy-tailed
inter-arrival times. The known "missing scales" pitfall of naive FRP construction
(Roughan, Yates & Veitch 1999) is addressed by using FSNP or a corrected FRP with a
verified scale range.

**MB4 — Hawkes-style symbolic process.**
`HawkesSymbol` is a discrete-time, finite-history analogue of Hawkes-process word
occurrence models. Each symbol has a positive baseline intensity, and each recent
symbol adds a power-law weighted row of an excitation matrix to the current
intensity vector. Identity-like excitation creates bursty repeat behavior; off-
diagonal excitation can encode cross-symbol triggering. The model is motivated by
Ogura, Hanada, Amano & Kondo (2022), who model long-range dynamic correlations of
word occurrences in written text with Hawkes processes. In S5.jl this is a
symbol-sequence generator rather than a fitted continuous-time text model.

**MB5 — Duplication-mutation symbolic growth.**
`DuplicationMutation` starts from an iid seed and grows the sequence by repeatedly
copying from a power-law-distributed previous lag and mutating copied symbols.
The method is motivated by expansion-modification and duplication-mutation models
for DNA-like symbolic sequences; it is a finite copy/mutate simulator, not a
biological genome model. The lag-copy mechanism is deliberate: validation showed
that power-law block lengths with uniformly chosen source blocks mostly created
local patches rather than the intended decaying autocorrelation curve.

---

## Controllability

All implemented generators accept an explicit ordered alphabet and reject duplicate
alphabet entries. `target_marginal(g)` reports the marginal distribution the generator
claims to target; `empirical_marginal(seq, alphabet)` and `empirical_bigram(seq,
alphabet)` provide lightweight checks for simulated data.

| Type | Alphabet | Marginal control | Bigram/trigram control |
|------|----------|------------------|------------------------|
| `SpectralFGN` | explicit `alphabet` | direct `marginal`; rank binning gives near-exact finite-sample counts | no direct control |
| `LGCM` | explicit `alphabet` | direct `marginal`; calibrated latent offsets | no direct control |
| `WaveletMarkov` | explicit `alphabet` | aggregate stationary marginal implied by regimes | direct per-regime bigram matrices |
| `IntermittentMapSymbols` | explicit `alphabet` | direct `marginal`; rank binning gives near-exact finite-sample counts | no direct control |
| `LAMP` | explicit `alphabet` | direct `marginal` mixed through `epsilon`; larger `epsilon` improves marginal control but weakens history dependence | exact history-weighted transition matrix |
| `DyadicLAMP` | explicit `alphabet` | direct `marginal` mixed through `epsilon`; larger `epsilon` improves marginal control but weakens history dependence | dyadic-bucket approximation to history-weighted transition matrix |
| `CalibratedAdditiveMarkov` | explicit `alphabet` | centered additive memory around `marginal`; `strength = 0` is iid | additive memory induces recurrence, not arbitrary bigrams |
| `OnOffMarkov` | explicit `alphabet` | aggregate stationary marginal implied by regimes | direct per-regime bigram matrices |
| `FSS` | explicit `alphabet` | `rates / sum(rates)` asymptotically | no direct control |
| `HawkesSymbol` | explicit `alphabet` | baseline distribution reported, but output marginal is implied by excitation and finite history | excitation matrix induces bursty local and long-context structure |
| `DuplicationMutation` | explicit `alphabet` | seed/mutation replacement distribution; output marginal is shaped by copied history | copy/mutate growth induces local and long-context structure |

For regime-driven methods (`WaveletMarkov` and `OnOffMarkov`), symbol-level ACF
and spectrum diagnostics only see the LRD regime process when regimes have
different observable stationary symbol distributions. Regimes with identical
stationary marginals can carry latent long memory while looking nearly
short-memory to one-hot symbol diagnostics.

Reproducible simulation studies live in `validation/`. For example:

```julia
julia --project=. validation/marginal_control.jl
julia --project=. validation/lrd_method_diagnostics.jl
julia --project=validation validation/longmemory_comparison.jl
```

These studies test controllability of simulated data; LRD-parameter estimation is
intended for a future separate estimator package.

The LRD diagnostic transformation is formalized in
`validation/lrd_symbol_diagnostics.jl`: symbols are converted to centered one-hot
numeric series before autocorrelation, autocovariance, and periodogram
calculations. The LongMemory.jl comparison script documents and tests the needed
adaptations to LongMemory.jl's lag and frequency conventions.
Autocorrelation validation plots include dashed vertical interpretation limits:
a finite-sample `n / 10` lag limit, and explicit generator limits where they
exist, such as `LAMP.d`. Spectrum plots show the same scales as reciprocal
frequencies. Where the model has a defensible asymptotic-onset scale, plots also
mark an approximate power-law onset. For example, `OnOffMarkov` uses its Pareto
scale `L_min`, while `HawkesSymbol` uses the lag where the offset kernel
`(lag + c)^(-beta)` reaches 90% of its asymptotic log-log slope. Autocorrelation
plots also include a gray dashed nominal power-law reference line with slope
`lag^(-beta)`, anchored to the first positive plotted autocorrelation value.
Power-spectrum plots include the corresponding gray dashed low-frequency
reference with slope `frequency^(beta - 1)`.

Current MB4 (`HawkesSymbol`) diagnostics should be read cautiously: the finite
discrete-time implementation can produce short-range burstiness while its
centered one-hot power spectrum remains close to white noise. Improving MB4
probably requires a more faithful near-critical/event-count Hawkes construction,
not just stronger identity excitation.

Validation policy is documented in `VALIDATION_POLICY.md`. The fast package test
suite remains the main development pathway, while larger empirical studies are run
manually or behind explicit flags such as `S5_VALIDATION_LARGE=true`.

Benchmarks live in `benchmark/` and use a separate `Project.toml` with
BenchmarkTools.jl:

```julia
julia --project=benchmark benchmark/benchmarks.jl
S5_BENCHMARK_LARGE=true julia --project=benchmark benchmark/benchmarks.jl
```

First-order local-structure controls use `MarkovSpec`. Trigram diagnostics are
available through `empirical_trigram`, but a concrete trigram-control
specification is intentionally left for future work. The code now exposes
`LocalStructureSpec` and `local_structure_order` as the extension path for that
higher-order API.



---

## References

- Paxson, V. (1997). Fast, approximate synthesis of fractional Gaussian noise. *CCR* 27.
- Dieker, T. (2004). *Simulation of fractional Brownian motion*. PhD thesis, U. Twente.
- Roughan, M., Veitch, D., & Abry, P. (2000). Real-time estimation of LRD parameters. *IEEE/ACM ToN* 8(4).
- Gal, Y., Chen, Y., & Ghahramani, Z. (2015). Latent Gaussian processes for distribution estimation of multivariate categorical data. *ICML*.
- Li, W. (1991). Expansion-modification systems: a model for spatial 1/f spectra. *Physical Review A* 43.
- Li, W., Marr, T. G., & Kaneko, K. (1994). Understanding long-range correlations in DNA sequences. *Physica D* 75.
- Melnyk, S. S., Usatenko, O. V., & Yampol'skii, V. A. (2006). Memory functions of the additive Markov chains. *Physica A* 361.
- Mayzelis, Z. A., Apostolov, S. S., Melnyk, S. S., Usatenko, O. V., & Yampol'skii, V. A. (2006). Additive N-step Markov chains as prototype model of symbolic stochastic dynamical systems with long-range correlations.
- Kumar, R., Raghu, M., Sarlós, T., & Tomkins, A. (2017). Linear additive Markov processes. *WWW '17*.
- Singh, M., Greenberg, C., & Klakow, D. (2016). The custom decay language model. *TSD*.
- Ryu, B., & Lowen, S. (1998). Point process models for self-similar network traffic. *Stochastic Models* 14(3).
- Roughan, M., Yates, J., & Veitch, D. (1999). The mystery of the missing scales. *Heavy Tails Workshop*.
- Provata, A., & Beck, C. (2012). Coupled intermittent maps modelling the statistics of genomic sequences: a network approach. arXiv:1205.2249.
- Pipiras, V., & Taqqu, M. S. (2017). *Long-Range Dependence and Self-Similarity*. Cambridge UP.

---

## AI Disclosure

This package was developed with assistance from **Claude Sonnet 4.6** (Anthropic), and Codex; 
AI coding assistants. The design goals, overall architecture, methods design were human, but 
Claude and Codex contributed to the design of the package architecture,
the coding itself, and some of the write-up of synthesis methods. 
