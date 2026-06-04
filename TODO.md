# S5.jl TODO

This package is the synthesis side of the ARC Discovery Grant project
"Analysis and Synthesis of Long-Range Structure in Non-Numerical Time Series"
(Roughan & Willinger, 2023).

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
- [x] `LAMP` (MB1): Linear-Additive Markov Process.
- [x] `FSS` (MB3): Fractal Symbol Sequence via independent Pareto renewal streams.
- [x] INC output with provenance metadata via `save_sequence`.
- [x] Basic unit tests for construction, output length/type, alphabet membership, and
      simple marginal checks.
- [x] Documenter docs and changelog entry for v0.1.0.
- [x] Stable reproducibility support via StableRNGs in tests and validation studies.
- [x] Random-variable support via Distributions.jl for Pareto renewal sampling.

Not yet implemented:

- [ ] PB2: Latent Gaussian Categorical Model.
- [ ] PB3: Wavelet-cascade driving a Markov state machine.
- [ ] MB2: Heavy-tailed On/Off doubly-stochastic Markov chain.
- [ ] Benchmark suite.
- [x] Initial simulation studies of marginal controllability.
- [ ] Expanded simulation studies of local-structure controllability.

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
- [x] `LAMP(beta, alphabet, marginal = uniform; d = 1000)`.
- [x] `FSS(alpha, alphabet; rates = ones(k), x_min = 1.0)`.
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
- `LAMP`: accepts `marginal` and mixes it into the history-based probabilities through
  `epsilon`. Larger `epsilon` improves finite-sample marginal control but weakens
  history dependence.
- `FSS`: accepts `rates`; target marginal is `rates / sum(rates)` asymptotically.

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
- `LAMP`: controls dependence through history weights, but not arbitrary user-specified
  bigram/trigram probabilities in the current implementation.
- `FSS`: no direct bigram/trigram control because symbol streams are independent.
- `PB3` and `MB2`: natural places to support Markov transition matrices and therefore
  bigram control.

Tasks:

- [ ] Define a short-range specification type, for example:
      `MarkovSpec(alphabet, transition_matrix)` for bigram control.
- [ ] Consider a higher-order form for trigram control, e.g. a sparse mapping from
      `(previous_symbol_1, previous_symbol_2)` to next-symbol probabilities.
- [x] Add validation helpers for local structure:
      empirical unigram, bigram, and trigram frequency tables;
      total variation distance from a target table;
      row-wise transition error for Markov matrices.
- [ ] Add capability docs: each generator should state whether it supports
      `:marginal`, `:bigram`, `:trigram`, and how strongly.
- [ ] Use these specs first in `MB2`, then in `PB3`.

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
- [ ] `test/local_structure.jl`: empirical bigram/trigram tools, initially tested on
      simple iid and Markov baselines.
- [ ] Decide whether large simulations belong under `test/` with `@testset`, under
      `validation/`, or under `examples/` as reproducible scripts.
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

- [ ] Add `benchmark/benchmarks.jl` using BenchmarkTools.jl.
- [ ] Measure wall time and allocations for `n in (10^4, 10^5, 10^6)`.
- [ ] Include alphabet sizes `k in (2, 8, 64)` where feasible.
- [ ] Include skewed marginal/rate settings.
- [ ] Report complexity-relevant parameters:
      `d` for `LAMP`, `k` for `FSS`, FFT length for `SpectralFGN`.

---

## Priority 4: Next Generators

### MB2: Heavy-Tailed On/Off Doubly-Stochastic Markov Chain

Implement next if the goal is controllable local structure.

- [ ] Define constructor accepting `alphabet`, `transition_matrices`,
      `switching_matrix`, `alpha`, and `L_min`.
- [ ] Validate Markov matrices: square, row-stochastic, non-negative, matching alphabet.
- [ ] Generate regime sojourns from a heavy-tailed distribution.
- [ ] Emit symbols from the active regime's transition matrix.
- [ ] Test marginal and bigram control within regimes and in aggregate.
- [ ] Document that aggregate marginals depend on regime occupancy and transition
      stationary distributions.

### PB3: Wavelet-Cascade Driving a Markov State Machine

Implement after the local-structure specification is settled.

- [ ] Reuse `MarkovSpec` or equivalent transition-matrix interface.
- [ ] Generate a latent LRD driver.
- [ ] Map driver values to regimes.
- [ ] Emit via regime-specific transition matrices.
- [ ] Test controllability of bigrams conditional on regime and in aggregate.

### PB2: Latent Gaussian Categorical Model

Implement when we want another property-based baseline.

- [ ] Accept `alphabet` and `marginal`.
- [ ] Implement marginal calibration through latent means or thresholds.
- [ ] Start with a practical FFT/circulant approximation rather than exact Cholesky for
      large `n`.
- [ ] Test marginal control by simulation.

---

## Documentation Tasks

- [x] Update README to describe controllability:
      alphabet, marginal, bigram/trigram, and LRD mechanism separately.
- [x] Add a capability matrix to docs:
      method versus alphabet/marginal/bigram/trigram/LRD parameter controls.
- [x] Explain that strong LRD-parameter validation is deferred to external estimator
      packages.
- [x] Add examples for custom alphabets and non-uniform marginals.
- [ ] Add examples of what cannot be controlled by each method.

---

## Release Hygiene

- [ ] Keep CHANGELOG.md updated as new controls, tests, and methods are added.
- [ ] Ensure CI runs fast unit tests.
- [ ] Decide whether longer simulation studies run manually, nightly, or behind an
      environment variable.
- [ ] Revisit version number after the TODO and docs match the implemented package.
