# S5.jl

**Self-Similar Symbols Sequence Synthesis**

S5.jl generates Long-Range Dependent (LRD) sequences of categorical
(non-numerical) symbols. Its primary use is producing ground-truth test data
for LRD estimation algorithms applied to symbol sequences such as natural
language, genomic data, or event logs.

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

# Model-based: Linear-Additive Markov Process
g3 = LAMP(0.5, [:a, :b, :c], [0.2, 0.3, 0.5]; d = 500, epsilon = 0.02)
seq3 = generate(g3, 10_000; rng)

# Model-based: heavy-tailed regime-switching Markov chain
P1 = [0.9 0.1; 0.2 0.8]
P2 = [0.3 0.7; 0.6 0.4]
Q = [0.2 0.8; 0.8 0.2]
g4 = OnOffMarkov(1.5, [:a, :b], [P1, P2], Q)
seq4 = generate(g4, 10_000; rng)

# Model-based: Fractal Symbol Sequence  (H = (3-α)/2 = 0.75)
g5 = FSS(1.5, [:a, :b, :c]; rates = [2.0, 3.0, 5.0])
seq5 = generate(g5, 10_000; rng)

empirical_marginal(seq1, g1.alphabet)

# Save to INC format with full provenance metadata
save_sequence("seq_pb1.inc", seq1, g1)
```

## Generators

| ID  | Type          | LRD mechanism                          | Short-range control        | Complexity      |
|-----|---------------|----------------------------------------|----------------------------|-----------------|
| PB1 | `SpectralFGN` | Spectral $1/f^\alpha$ shaping          | Poor (set by quantization) | $O(n \log n)$  |
| PB2 | `LGCM`        | Latent fGn streams + argmax            | Offset-calibrated marginals | $O(n \cdot k \cdot I)$ |
| MB1 | `LAMP`        | Power-law history weights              | Weight tensor              | $O(n \cdot d)$ |
| MB2 | `OnOffMarkov` | Heavy-tailed regime sojourns           | Per-regime Markov matrices | $O(n \cdot k)$ |
| MB3 | `FSS`         | Pareto renewal process per symbol      | Poor (independent streams) | $O(n \cdot k)$ |

### LRD parameters

All three generators expose a Hurst parameter $H \in (1/2, 1)$:

| Type          | Input parameter  | Relationship to $H$    |
|---------------|------------------|------------------------|
| `SpectralFGN` | `H`              | direct                 |
| `LGCM`        | `H`              | direct                 |
| `LAMP`        | `beta` ($\beta$) | $H = (2 - \beta) / 2$ |
| `OnOffMarkov` | `alpha` ($\alpha$) | nominal $H = (3 - \alpha) / 2$ |
| `FSS`         | `alpha` ($\alpha$) | $H = (3 - \alpha) / 2$ |

## Controllability

All implemented generators accept an explicit ordered alphabet. Duplicate alphabet
entries are rejected because they make empirical frequency tables ambiguous.

| Type | Marginal control | Local structure control |
|------|------------------|-------------------------|
| `SpectralFGN` | direct `marginal`; rank binning gives integer counts as close as possible to target | none |
| `LGCM` | direct `marginal`; calibrated latent offsets | none |
| `LAMP` | direct `marginal` mixed through `epsilon`; larger `epsilon` improves marginal control | history-weighted dependence, not arbitrary bigrams |
| `OnOffMarkov` | aggregate stationary marginal implied by regimes | per-regime bigram matrices |
| `FSS` | asymptotic `rates / sum(rates)` | none |

Reproducible controllability studies live in `validation/`, for example:

```julia
julia --project=. validation/marginal_control.jl
```

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
