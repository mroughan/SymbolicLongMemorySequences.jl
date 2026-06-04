# S5.jl

**Self-Similar Symbols Sequence Synthesis**

S5.jl is a Julia package for generating Long-Range Dependent (LRD) sequences of
non-numerical (categorical/symbolic) data. It provides ground-truth synthetic sequences
for use in testing and benchmarking LRD estimation algorithms applied to non-numerical
time series.

This package implements Task 1 of the ARC Discovery Grant project *"Analysis and
Synthesis of Long-Range Structure in Non-Numerical Time Series"* (Roughan & Willinger).

---

## Background

Long-Range Dependence (LRD) means that the large-scale statistical structure of a
sequence is as important as its short-range structure. Formally, a sequence is LRD if
its autocovariance function (ACF) decays as a power law,

$$\gamma_k \sim c_\gamma |k|^{-\beta}, \quad \beta \in (0, 1),$$

so that the sum of the tail ACF diverges. LRD is characterised by two closely related
parameters: the ACF decay exponent $\beta$, the spectral exponent $\alpha = 1 - \beta$,
and the Hurst parameter $H = (2 - \beta)/2 \in (1/2, 1)$.

LRD is ubiquitous in human-generated data (text, Internet traffic, genomics, social
media), yet almost all synthesis tools target numerical data. S5.jl fills that gap for
**symbol sequences** — data that takes values in a finite, unordered alphabet such as
`{a, b, c, ...}`.

---

## Implemented Methods

Five synthesis methods are currently implemented. Additional methods from the
grant proposal remain on the roadmap.

### Property-Based Methods

These generate an underlying numerical LRD process and then map it to symbols.
The LRD property is inherited from the numerical layer.

| ID | Name | Status | LRD mechanism | Short-range control | Complexity | Novel? |
|----|------|--------|---------------|---------------------|------------|--------|
| PB1 | Spectral fGn + quantization | Implemented | Spectral $1/f^\alpha$ shaping | Poor (set by quantization) | $O(n \log n)$ | No |
| PB2 | Latent Gaussian categorical (LGCM) | Implemented | fGn streams + argmax | Via calibrated offsets | $O(n k I)$ | No |
| PB3 | Wavelet-cascade + Markov state machine | Planned | Wavelet coefficient cascade | Markov transition matrices | $O(n)$ | Partial |

**PB1 — Spectral fGn + quantization.**
Fractional Gaussian noise with Hurst parameter $H$ is synthesised using the fast
spectral method of Paxson (1997) or circulant embedding (Dieker 2004). The real-valued
output is sorted into $k$ rank bins; each bin maps to one symbol, with integer bin
counts chosen to match a target marginal distribution as closely as possible for the
finite sample. This is the simplest approach and serves as the primary validation
baseline.

**PB2 — Latent Gaussian categorical model (LGCM).**
A vector of $k$ latent fGn streams is generated, one stream per symbol. At each time
step the symbol is the argmax of the latent vector plus calibrated per-symbol offsets.
The offsets shift marginal probabilities while the latent streams carry the LRD
structure. This is a practical finite-sample approximation to the latent Gaussian
categorical model of Gal, Chen & Ghahramani (ICML 2015).

**PB3 — Wavelet-cascade driving a Markov state machine.**
A latent LRD intensity signal is generated via the wavelet synthesis method of Roughan,
Veitch & Abry (2001). This signal controls which row of a symbol transition matrix is
active at each step, so the state machine enforces prescribed bigram/trigram statistics
while the wavelet cascade injects LRD at all scales. Provides simultaneous control of
both $H$ and short-range structure at $O(n)$ cost.

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
| `LAMP` | explicit `alphabet` | direct `marginal` mixed through `epsilon`; larger `epsilon` improves marginal control but weakens history dependence | no arbitrary target table |
| `OnOffMarkov` | explicit `alphabet` | aggregate stationary marginal implied by regimes | direct per-regime bigram matrices |
| `FSS` | explicit `alphabet` | `rates / sum(rates)` asymptotically | no direct control |

Reproducible simulation studies live in `validation/`. For example:

```julia
julia --project=. validation/marginal_control.jl
```

These studies test controllability of simulated data; LRD-parameter estimation is
intended for a future separate estimator package.

---

## Implementation Priority

Methods are implemented in the following order:

1. **PB1** — fastest baseline; easiest to validate against known $H$ via spectral estimators.
2. **MB1** — cleanest finite-state theory; already demonstrated on text corpora.
3. **MB3** — novel FSS contribution that distinguishes this project scientifically.
4. **PB2 and MB2** — additional marginal and local-structure controls.
5. **PB3** — follows once wavelet-regime coupling is specified.

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

This package was developed with assistance from **Claude Sonnet 4.6** (Anthropic), an
AI coding assistant. Claude contributed to the design of the package architecture,
the selection and write-up of synthesis methods, and the generation of documentation
and code. All scientific content is grounded in the references cited above and in the
ARC Discovery Grant proposal by Roughan & Willinger (2023). Human review and direction
was provided throughout by the project investigators.
