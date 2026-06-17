# S5.jl TODO

Estimator development will live elsewhere. S5.jl should therefore focus on:

- generating symbol sequences with explicit provenance;
- exposing clear controls for alphabet, marginal distribution, and local structure;
- measuring how well each generator achieves the controls it claims to support;
- providing reproducible simulation studies that external estimator packages can reuse.

---

## Current Status

Implemented:

- [x] Common generator interface: `generate(g, n; rng)`.
- [x] `SpectralFGN` (PB1): spectral fGn plus quantization.
- [x] `LGCM` (PB2): latent Gaussian categorical model.
- [x] `WaveletMarkov` (PB3): multiscale latent driver with Markov regimes.
- [x] `LAMP` (MB1a): exact finite-history Linear-Additive Markov Process.
- [x] `DyadicLAMP` (MB1b): scalable dyadic-bucket LAMP approximation.
- [x] `OnOffMarkov` (MB2): heavy-tailed regime-switching Markov chain.
- [x] `FSS` (MB3): Fractal Symbol Sequence via independent Pareto renewal streams.
- [x] `HawkesSymbol` (MB4): finite-history Hawkes-style self-exciting symbol process.
- [x] `IntermittentMapSymbols` (PB4): intermittent-map latent driver plus quantization.
- [x] `CalibratedAdditiveMarkov` (MB1c): centered additive Markov-chain memory
      function.
- [x] `DuplicationMutation` (MB5): power-law lag copy/mutate symbolic growth
      generator.
- [x] INC output with provenance metadata via `save_sequence`.
- [x] Basic unit tests for construction, output length/type, alphabet membership, and
      simple marginal checks.
- [x] Documenter docs and changelog entry for v0.1.0.
- [x] Stable reproducibility support via StableRNGs in tests and validation studies.
- [x] Random-variable support via Distributions.jl for Pareto renewal sampling.

Not yet implemented:

- [x] Benchmark suite scaffolded under `benchmark/` with a separate
      BenchmarkTools.jl environment and an opt-in large-run flag.
- [x] Initial simulation studies of marginal controllability.
- [x] Expanded simulation studies of local-structure controllability.
- [x] Formalized LRD visual diagnostic transformations and added a
      LongMemory.jl comparison validation script.

Deferred to a future estimation package:

- Spectral, wavelet, Whittle, recurrence-time, Hill, and count-variance estimators of
  LRD parameters.
- Strong validation that generated sequences recover a target Hurst parameter under a
  particular estimator.

S5.jl may include light sanity checks that guard against obvious regressions, but it
should not grow into the estimator package.

---

## Priority 1: Control Contracts

Define, document, and test what each generator can naturally control.

### Alphabet Control

All generators should accept an ordered `alphabet` collection and emit elements from
that collection with the same element type.

- [x] `SpectralFGN(H, alphabet, marginal = uniform)`.
- [x] `LGCM(H, alphabet, marginal = uniform)`.
- [x] `WaveletMarkov(H, alphabet, transition_matrices; regime_weights = uniform)`.
- [x] `LAMP(beta, alphabet, marginal = uniform; d = 1000)`.
- [x] `DyadicLAMP(beta, alphabet, marginal = uniform; d = 1_000_000)`.
- [x] `OnOffMarkov(alpha, alphabet, transition_matrices, switching_matrix; L_min = 1.0)`.
- [x] `FSS(alpha, alphabet; rates = ones(k), x_min = 1.0)`.
- [x] `HawkesSymbol(beta, alphabet; baseline = ones(k), excitation = I, d = 1000)`.
- [x] Add explicit tests for non-symbol alphabets:
      `Char`, `String`, `Int`, `Symbol`, and small custom immutable values if useful.
- [x] Decide whether duplicate alphabet entries should be rejected.
- [x] Add a shared alphabet-validation helper if the constructors start duplicating
      checks.

### Marginal Control

Users should be able to specify target marginal probabilities over the alphabet when a
method has a meaningful marginal-control knob.

Current behavior:

- `SpectralFGN`: accepts `marginal`; rank binning gives integer finite-sample counts
  as close as possible to the requested marginal.
- `LGCM`: accepts `marginal`; latent offsets are calibrated on the generated sample
  to approximate the target marginal.
- `WaveletMarkov`: aggregate marginal is implied by regime weights and per-regime
  stationary distributions.
- `LAMP`: accepts `marginal` and mixes it into exact history-based probabilities
  through `epsilon`. Larger `epsilon` improves finite-sample marginal control but
  weakens history dependence.
- `DyadicLAMP`: accepts the same `marginal`, `epsilon`, and transition matrix
  controls as `LAMP`, but approximates long histories with dyadic age buckets.
- `OnOffMarkov`: aggregate marginal is implied by regime occupancy and per-regime
  stationary distributions.
- `FSS`: accepts `rates`; target marginal is `rates / sum(rates)` asymptotically.
- `HawkesSymbol`: accepts `baseline`; `target_marginal` reports the normalized
  baseline, but realized marginals are implied by excitation and finite history.

Tasks:

- [x] Add `target_marginal(g)` and/or `marginal_control(g)` helpers so tests and docs
      can ask a generator what it claims to control.
- [x] Add Monte Carlo marginal-control tests across many independent sequences:
      compare mean empirical frequencies, standard deviations, and worst-case errors
      to the declared target.
- [ ] Test marginal control across alphabet sizes, e.g. `k in (2, 4, 16, 64)`.
- [x] Test skewed marginals, including Zipf-like distributions, because this is central
      for text and word-token sequences.
- [x] For `SpectralFGN`, test that empirical marginals are close even per sequence.
- [x] For `FSS`, test convergence of frequencies as `n` and replicate count increase.
- [x] For `LAMP`, investigate whether the current process needs an innovation or
      teleportation term, e.g.
      `q = (1 - epsilon) * history_weights + epsilon * marginal`,
      to make marginal control robust and prevent finite-history absorption.

Suggested simulation grid:

- `n in (1_000, 10_000, 100_000)`;
- `replicates in (50, 200)` depending on runtime;
- `k in (2, 8, 32)`;
- marginals: uniform, moderately skewed, Zipf-like.

Store aggregate results, not every generated sequence.

### Local Structure Control

Users may want to specify short-range structure such as bigram or trigram
probabilities. Not every generator can support this naturally, so S5.jl should expose
capabilities rather than pretending all methods can do everything.

Current capability:

- `SpectralFGN`: no direct bigram/trigram control. Local structure is induced by the
  latent Gaussian path and quantization thresholds.
- `LGCM`: no direct bigram/trigram control. Local structure is induced by the latent
  fGn streams and argmax mapping.
- `WaveletMarkov`: direct per-regime bigram control through Markov transition
  matrices driven by a multiscale latent regime process.
- `LAMP`: controls dependence through exact history weights and a transition matrix,
  but not arbitrary user-specified bigram/trigram probabilities.
- `DyadicLAMP`: controls dependence through dyadic-bucket approximations to the
  same history-weighted transition mechanism.
- `OnOffMarkov`: direct per-regime bigram control through Markov transition matrices.
- `FSS`: no direct bigram/trigram control because symbol streams are independent.
- `HawkesSymbol`: local and long-context structure are induced by the excitation
  matrix and power-law history kernel, not by arbitrary target bigram/trigram tables.

Tasks:

- [x] Define a short-range specification type, for example:
      `MarkovSpec(alphabet, transition_matrix)` for bigram control.
- [ ] Consider a higher-order form for trigram control, e.g. a sparse mapping from
      `(previous_symbol_1, previous_symbol_2)` to next-symbol probabilities.
      The code now exposes `LocalStructureSpec` and `local_structure_order` so a
      future trigram specification has a clear extension point without implying
      current trigram control.
- [x] Add validation helpers for local structure:
      empirical unigram, bigram, and trigram frequency tables;
      total variation distance from a target table;
      row-wise transition error for Markov matrices.
- [x] Add capability docs: each generator should state whether it supports
      `:marginal`, `:bigram`, `:trigram`, and how strongly.
- [x] Use these specs first in `MB2`, then in `PB3`.

---

## Priority 2: Empirical Controllability Tests

These tests do not estimate LRD. They test whether the generator respects user-facing
controls.

### Unit Tests

- [x] Constructor rejects invalid marginals: wrong length, negative entries, all-zero,
      and non-finite values.
- [x] Constructor rejects invalid rates for `FSS`: wrong length, non-positive, and
      non-finite values.
- [x] Generated sequence uses only the supplied alphabet.
- [x] Generated sequence preserves expected element type.
- [x] Reproducibility with fixed RNG seeds.

### Simulation Tests

Keep these separate from fast unit tests if runtime becomes large.

- [x] `test/marginal_control.jl`: repeated simulations for each implemented generator.
- [x] `test/local_structure.jl`: empirical bigram/trigram tools, initially tested on
      simple iid and Markov baselines.
- [x] Decide whether large simulations belong under `test/` with `@testset`, under
      `validation/`, or under `examples/` as reproducible scripts. The policy is
      documented in `VALIDATION_POLICY.md`: fast contract checks belong in `test/`;
      larger empirical studies belong in `validation/` and run manually or behind
      explicit flags.
- [x] Move visual LRD diagnostic mechanics into reusable validation code, including
      centered one-hot transformations and LongMemory.jl-compatible autocovariance,
      autocorrelation, and periodogram conventions.
- [x] Add dashed interpretation limits to LRD validation plots: finite-sample
      `n / 10` autocorrelation support, reciprocal spectrum scale, and explicit
      generator memory limits such as `LAMP.d`.
- [ ] Write summary tables to `data/validation/results/` only if those outputs are
      intended to be tracked.

Useful metrics:

- empirical marginal vector;
- absolute error and total variation distance from target marginal;
- distribution of errors across replicates;
- confidence intervals for Monte Carlo error;
- bigram/trigram total variation distance where a target exists;
- runtime and allocation counts.

---

## Priority 3: Benchmarks

- [x] Add `benchmark/benchmarks.jl` using BenchmarkTools.jl.
- [x] Measure wall time and allocations for `n in (10^4, 10^5)`, with
      `10^6` available through `S5_BENCHMARK_LARGE=true`.
- [x] Include alphabet sizes `k in (2, 8, 64)` where feasible.
- [x] Include skewed marginal/rate settings.
- [x] Report complexity-relevant parameters:
      `d` for `LAMP` and `HawkesSymbol`, `k` for `FSS`, FFT length for `SpectralFGN`.

---

## Priority 4: Next Generators

### Candidate Future Generators From The Literature

Recent bibliography work suggests three plausible generator directions beyond
the current PB/MB set:

- [x] **Additive Markov chain with prescribed memory function.** Melnyk,
      Usatenko, and Yampol'skii (2006), Mayzelis et al. (2006), and Melnik and
      Usatenko (2014) give a finite-alphabet additive Markov-chain route from
      desired pair correlations to memory functions. This is the strongest
      candidate for a more principled successor or extension to MB1 because it
      could expose explicit target-correlation controls while preserving a
      symbolic generator identity. Implemented initially as
      `CalibratedAdditiveMarkov`, with full correlation calibration left for
      future refinement.
- [x] **Expansion-modification / duplication-mutation symbolic generator.** Li's
      expansion-modification model and later duplication-mutation DNA work
      provide a natural symbolic mechanism for DNA-like sequences. This is a
      good candidate if S5.jl adds a generator whose scientific identity is
      copy-mutate growth rather than latent Gaussian, LAMP, or renewal-process
      synthesis. Implemented initially as `DuplicationMutation`; the current
      version uses a power-law copy-distance kernel after validation showed that
      power-law block lengths with uniform source selection gave nearly flat
      autocorrelation.
- [ ] **Hierarchical text/rare-event generator.** Alvarez-Lacalle et al. (2006),
      Altmann et al. (2009, 2012), Tanaka-Ishii and Bunde (2016), and Ogura et
      al. (2022) suggest text-oriented mechanisms based on document hierarchy,
      rare-word clustering, and bursty point processes. This may belong as a
      validation/application model rather than as a core generator unless the
      controls can be stated cleanly.
- [x] **Intermittent-map symbolic generator.** Provata and Beck (2012) show that
      coupled intermittent maps can create DNA-like symbolic statistics. This is
      less immediate for S5.jl because map parameters do not directly translate
      to alphabet, marginal, bigram, and LRD control contracts, but it is worth
      keeping as a lower-priority research candidate. Implemented initially as
      `IntermittentMapSymbols`, with the API explicitly framed as a latent
      empirical driver.

### MB2: Heavy-Tailed On/Off Doubly-Stochastic Markov Chain

Implemented as `OnOffMarkov`.

- [x] Define constructor accepting `alphabet`, `transition_matrices`,
      `switching_matrix`, `alpha`, and `L_min`.
- [x] Validate Markov matrices: square, row-stochastic, non-negative, matching alphabet.
- [x] Generate regime sojourns from a heavy-tailed distribution.
- [x] Emit symbols from the active regime's transition matrix.
- [x] Test marginal and bigram control within regimes and in aggregate.
- [x] Document that aggregate marginals depend on regime occupancy and transition
      stationary distributions.

### PB3: Wavelet-Cascade Driving a Markov State Machine

Implemented as `WaveletMarkov`.

- [x] Reuse transition-matrix interface.
- [x] Generate a latent LRD driver.
- [x] Map driver values to regimes.
- [x] Emit via regime-specific transition matrices.
- [x] Test controllability of bigrams conditional on regime and in aggregate.
- [x] Revisit the latent LRD driver. `WaveletMarkov` now separates the latent
      regime driver from the Markov emission layer, defaults to spectral fGn
      rank-binning, and retains the original Haar cascade as `driver = :haar`
      for comparison.
- [ ] Compare PB3 spectral and Haar drivers in the full LRD diagnostics and
      decide whether the Haar path should remain only as a validation baseline.
- [ ] Investigate a fully calibrated wavelet synthesis driver for PB3.
- [ ] Add a PB3 driver-layer validation that separately plots latent-driver,
      regime-indicator, and emitted-symbol autocorrelation. Current evidence
      suggests the legacy Haar path can be too flat, while the spectral path can
      be damped strongly by rank-binning and Markov emission.

### MB4: Hawkes-Style Symbolic Process

Implemented as `HawkesSymbol`, but current validation is weak.

- [ ] Rework or supplement MB4 with a more faithful near-critical/event-count
      Hawkes construction. The present finite discrete-time sampler uses
      probability-normalized intensities and can show short-range burstiness
      while its centered one-hot power spectrum remains close to white noise.
      Simply increasing identity excitation does not appear to recover the
      nominal low-frequency power-law slope.

### MB1: Linear-Additive Markov Process

Implemented as exact `LAMP` (MB1a) and approximate `DyadicLAMP` (MB1b).

- [x] Reassess whether fixed-depth LAMP should be labeled only as a finite-history
      approximation. With a hard history depth `d`, the generator cannot provide
      true asymptotic LRD beyond its configured memory scale.
- [x] Allow `d > n` in a finite-sequence honest way: only observed history
      contributes, while missing pre-history weight is assigned to the target
      marginal.
- [x] Add a history-weighted transition matrix to LAMP, with
      `lamp_repeat_transition` providing a simple identity/dyad mixture for
      repeat-biased behavior.
- [x] Prototype a scalable full-history or multiscale-history variant whose
      effective memory grows with `n`, with explicit provenance for the finite
      simulation cutoff. `DyadicLAMP` implements the dyadic-bucket path for
      power-law history weights.
- [ ] Validate MB1b against MB1a on small sequences and against LRD diagnostics on
      larger sequences where exact MB1a is too expensive.

### PB2: Latent Gaussian Categorical Model

Implemented as `LGCM`.

- [x] Accept `alphabet` and `marginal`.
- [x] Implement marginal calibration through latent means or thresholds.
- [x] Start with a practical FFT/circulant approximation rather than exact Cholesky for
      large `n`.
- [x] Test marginal control by simulation.

---

## Documentation Tasks

- [x] Update README to describe controllability:
      alphabet, marginal, bigram/trigram, and LRD mechanism separately.
- [x] Add a capability matrix to docs:
      method versus alphabet/marginal/bigram/trigram/LRD parameter controls.
- [x] Explain that strong LRD-parameter validation is deferred to external estimator
      packages.
- [x] Add examples for custom alphabets and non-uniform marginals.
- [x] Add examples of what cannot be controlled by each method.

---

## Release Hygiene

- [x] Keep CHANGELOG.md updated as new controls, tests, and methods are added.
- [x] Ensure CI runs fast unit tests.
- [x] Ensure Documenter CI resolves unregistered dependencies from a clean checkout.
- [x] Decide whether longer simulation studies run manually, nightly, or behind an
      environment variable. `VALIDATION_POLICY.md` keeps the main path fast and
      reserves larger studies for manual or flag-controlled runs.
- [ ] Revisit version number after the TODO and docs match the implemented package.
