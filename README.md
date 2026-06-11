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

### Property-Based Methods

These generate an underlying numerical LRD process and then map it to symbols.
The LRD property is inherited from the numerical layer.

| ID | Name | Status | LRD mechanism | Short-range control | Complexity | Novel? |
|----|------|--------|---------------|---------------------|------------|--------|
| PB1 | Spectral fGn + quantization | Implemented | Spectral $1/f^\alpha$ shaping | Poor (set by quantization) | $O(n \log n)$ | No |
| PB2 | Latent Gaussian categorical (LGCM) | Implemented | fGn streams + argmax | Via calibrated offsets | $O(n k I)$ | No |
| PB3 | Wavelet-cascade + Markov state machine | Implemented | Multiscale Haar-like cascade | Markov transition matrices | $O(n \log n + n k)$ | Partial |

**PB1 — Spectral fGn + quantization.**
Fractional Gaussian noise with Hurst parameter $H$ is synthesised using the fast,
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

**PB3 — Wavelet-cascade driving a Markov state machine.**
A latent multiscale Haar-like driver controls which Markov regime is active at each
step. Each regime has its own symbol transition matrix, so local bigram structure can
be prescribed while the latent cascade injects persistence across scales. This is a
practical implementation of the wavelet/state-machine idea in Roughan, Veitch & Abry
(2001), with the full calibrated wavelet variant left as a research extension.

---

### Model-Based Methods

These produce LRD through the stochastic model itself rather than via mapping.

| ID | Name | Status | LRD mechanism | Short-range control | Complexity | Novel? |
|----|------|--------|---------------|---------------------|------------|--------|
| MB1 | Linear-Additive Markov Process (LAMP) | Implemented | Power-law history weights | Weight tensor | $O(n \cdot d)$ | No |
| MB2 | Heavy-tailed On/Off doubly-stochastic Markov chain | Implemented | Pareto regime sojourn times | Per-regime Markov chains | $O(n \cdot k)$ | No |
| MB3 | Fractal Symbol Sequence (FSS) via FRP/FSNP | Implemented | Fractal point process inter-arrivals | Poor (independent streams) | $O(n \cdot k)$ | **Yes** |

**MB1 — Linear-Additive Markov Process (LAMP).**
Transition probabilities are a weighted sum over the history,

$$q(s) = (1-\epsilon)\sum_{k=1}^{d} w_k \cdot \mathbf{1}[X_{t-k} = s] + \epsilon p(s),$$

with weights $w_k \propto k^{-(1+\beta)}$ enforcing a power-law ACF decay directly.
The small innovation term $\epsilon p(s)$ keeps the requested marginal distribution
active after initialisation and prevents finite-history absorption.
The Custom Decay Language Model (CDLM) of Singh, Greenberg & Klakow (2016) is a close
variant demonstrated on text. For large alphabets the weight tensor may be compressed
via low-rank approximations.

**MB2 — Heavy-tailed On/Off doubly-stochastic Markov chain.**
The sequence is generated by a Markov chain that alternates between two or more regimes.
Sojourn times in each regime follow a Pareto distribution with tail index
$\alpha \in (1, 2)$, so the variance of the symbol count function grows super-linearly
(LRD by Definition 3 of the proposal), and $H = (3 - \alpha)/2$. Within each regime a
standard (SRD) Markov chain governs symbol emissions, giving direct control of local
statistics. Analogous to Fractal Shot Noise Processes adapted to symbol sequences
(Ryu & Lowen 1998; Garrett & Willinger 1994).

**MB3 — Fractal Symbol Sequence (FSS) via FRP/FSNP.**
Each symbol $s_i$ is assigned an independent Fractal Renewal Process (FRP) or Fractal
Shot Noise Process (FSNP) governing the times at which that symbol is emitted. The
final sequence merges all symbol streams, with the earliest pending event at each step
determining the output. LRD arises in each symbol's count process through heavy-tailed
inter-arrival times (Definition 4 of the proposal). The known "missing scales" pitfall
of naive FRP construction (Roughan, Yates & Veitch 1999) is addressed by using FSNP
or a corrected FRP with a verified scale range. This is the novel method proposed
specifically within the grant project.

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
| `LAMP` | explicit `alphabet` | direct `marginal` mixed through `epsilon`; larger `epsilon` improves marginal control but weakens history dependence | no arbitrary target table |
| `OnOffMarkov` | explicit `alphabet` | aggregate stationary marginal implied by regimes | direct per-regime bigram matrices |
| `FSS` | explicit `alphabet` | `rates / sum(rates)` asymptotically | no direct control |

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
- Roughan, M., Veitch, D., & Abry, P. (2001). Real-time estimation of LRD parameters. *IEEE/ACM ToN* 8(4).
- Gal, Y., Chen, Y., & Ghahramani, Z. (2015). Latent Gaussian processes for distribution estimation of multivariate categorical data. *ICML*.
- Kumar, R., Raghu, M., Sarlos, T., & Tomkins, A. (2017). Linear additive Markov processes. *WWW '17*.
- Singh, M., Greenberg, C., & Klakow, D. (2016). The custom decay language model. *TSD*.
- Ryu, B., & Lowen, S. (1998). Point process models for self-similar network traffic. *Stochastic Models* 14(3).
- Roughan, M., Yates, J., & Veitch, D. (1999). The mystery of the missing scales. *Heavy Tails Workshop*.
- Pipiras, V., & Taqqu, M. S. (2017). *Long-Range Dependence and Self-Similarity*. Cambridge UP.

---

## AI Disclosure

This package was developed with assistance from **Claude Sonnet 4.6** (Anthropic), and Codex; 
AI coding assistants. The design goals, overall architecture, methods design were human, but 
Claude and Codex contributed to the design of the package architecture,
the coding itself, and some of the write-up of synthesis methods. 
