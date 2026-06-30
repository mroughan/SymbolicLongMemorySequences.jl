# API

SymbolicLongMemorySequences.jl has two public construction paths:

- use `make_generator` for standard, named cases across all implemented methods;
- compose a `LatentSource` with a `Symbolizer` for property-based studies that
  need to vary the numerical driver and the symbolization rule separately;
- use the method-specific constructors when a study needs full control of the
  scientific parameters.

Both paths return an `LRDGenerator` and both use the same generation call:

```julia
seq = generate(g, n; rng)
```

Always pass an explicit `rng` when a sequence needs to be reproducible.

Property-based generators also support an additive validation helper:

```julia
seq, latent = generate_with_latent(g, n; rng)
```

`latent` is a `width × n` numerical matrix containing the LRD driver before
symbolization. Use this when a study needs to diagnose the numerical source and
the symbolic transform separately.

## Standard Method Factory

The factory is the cleanest way to list methods, inspect their standard
parameters, and build a generator for examples, smoke tests, and benchmark grids.

```julia
using SymbolicLongMemorySequences, StableRNGs

alphabet = [:a, :b, :c]
g = make_generator(:PB1, alphabet; H = 0.8, marginal = [0.2, 0.3, 0.5])
seq = generate(g, 8; rng = StableRNG(42))
```

Example output:

```julia
8-element Vector{Symbol}:
 :a
 :c
 :b
 :c
 :c
 :a
 :c
 :b
```

`id` may be a stable method identifier, a string, or an exported type-name alias:

```julia
make_generator(:PB1, alphabet)
make_generator("MB1c", alphabet)
make_generator(:SpectralFGN, alphabet)
```

Use `method_ids` and `method_info` to discover supported methods:

```julia
method_ids()
method_ids(family = :model_based)
method_info(:MB5).defaults
method_parameters(:MB5)
method_info(:SpectralFGN).id
```

The entries returned by `method_parameters(id)` describe the factory keywords
accepted by that method:

```julia
[(p.name, p.default, p.domain) for p in method_parameters(:PB1)]
```

Example output:

```julia
2-element Vector{Tuple{Symbol, Any, String}}:
 (:H, 0.8, "0.5 < H < 1")
 (:marginal, :uniform, "`:uniform` or a probability vector with one entry per symbol")
```

All factory methods also require the positional `alphabet` input. The generated
sequence length `n` is supplied later through `generate(g, n; rng)`.

Factory defaults are standard examples, not finite-sample scientific
calibrations. For validation studies that support method-specific claims, prefer
explicit constructors so the parameter assumptions remain visible in the script.

## Method IDs

| ID | Type | Standard cases |
|----|------|----------------|
| `:PB1` | `SpectralFGN` | `:standard` |
| `:PB2` | `LGCM` | `:standard` |
| `:PB3` | `WaveletMarkov` | `:persistent_regimes`, `:iid_regimes` |
| `:PB4` | `IntermittentMapSymbols` | `:standard` |
| `:MB1a` | `LAMP` | `:repeat`, `:iid` |
| `:MB1b` | `DyadicLAMP` | `:repeat`, `:iid` |
| `:MB1c` | `CalibratedAdditiveMarkov` | `:standard`, `:iid` |
| `:MB2` | `OnOffMarkov` | `:persistent_regimes`, `:iid_regimes` |
| `:MB3` | `FSS` | `:standard` |
| `:MB4` | `HawkesSymbol` | `:identity_excitation` |
| `:MB5` | `DuplicationMutation` | `:standard` |

Common keywords are shared where they mean the same thing:

- `marginal = :uniform` sets or approximates a target marginal when a method has
  a meaningful marginal-control knob;
- `case` selects one of the standard cases listed above;
- `H`, `beta`, `alpha`, or `z` select the method's nominal LRD-related
  parameter;
- `d` sets an explicit finite history cutoff for history-based methods.

## Composable Property-Based API

Property-based generators can be described as:

```julia
PropertyBasedGenerator(source, symbolizer)
```

The `source` creates one or more numerical latent LRD series. The `symbolizer`
turns those latent values into categorical symbols. This split makes the PB
taxonomy explicit: PB1, PB2, PB3, and PB4 are standard combinations rather than
unrelated mechanisms.

```julia
using SymbolicLongMemorySequences, StableRNGs

source = SpectralFGNSource(0.8)
symbolizer = QuantileSymbolizer([:a, :b, :c], [0.2, 0.3, 0.5])
g = PropertyBasedGenerator(source, symbolizer)

generate(g, 8; rng = StableRNG(42))
```

Example output:

```julia
8-element Vector{Symbol}:
 :a
 :c
 :b
 :c
 :c
 :a
 :c
 :b
```

Available latent sources are:

| Source | Latent behavior | Width |
|--------|-----------------|-------|
| `SpectralFGNSource(H)` | approximate spectral fGn | any positive width |
| `HaarLRDSource(H)` | Haar-like multiscale driver | any positive width |
| `IntermittentMapSource(z)` | intermittent-map dynamics | one stream |

Available symbolizers are:

| Symbolizer | Transformation | Required width |
|------------|----------------|----------------|
| `QuantileSymbolizer` | rank-bin one latent stream to target counts | 1 |
| `ArgmaxSymbolizer` | choose the largest offset latent stream | alphabet size |
| `MarkovRegimeSymbolizer` | rank-bin one latent stream to Markov regimes | 1 |

Some combinations are intentionally invalid. For example,
`IntermittentMapSource` cannot feed `ArgmaxSymbolizer`, because the intermittent
map source currently produces only one latent stream and argmax symbolization
needs one stream per symbol.

PB3's regime construction can also be built directly:

```julia
P1 = [0.95 0.05; 0.10 0.90]
P2 = [0.30 0.70; 0.70 0.30]
source = HaarLRDSource(0.8)
symbolizer = MarkovRegimeSymbolizer([:a, :b], [P1, P2])
g = PropertyBasedGenerator(source, symbolizer)
```

Use the composable API when the study is about transformations. Use
`make_generator` when the study only needs a standard named method.

To inspect the numerical driver used by a property-based generator:

```julia
seq, latent = generate_with_latent(g, 8; rng = StableRNG(42))
size(latent)
```

Example output:

```julia
(1, 8)
```

## Explicit Constructors

Use explicit constructors when a model needs exact transition matrices,
excitation matrices, rates, or calibration parameters.

```julia
using SymbolicLongMemorySequences, StableRNGs

P1 = [0.9 0.1; 0.2 0.8]
P2 = [0.3 0.7; 0.6 0.4]
Q = [0.2 0.8; 0.8 0.2]

g = OnOffMarkov(1.5, [:a, :b], [P1, P2], Q; L_min = 50.0)
generate(g, 6; rng = StableRNG(7))
```

Example output:

```julia
6-element Vector{Symbol}:
 :a
 :b
 :b
 :b
 :b
 :b
```

The constructors validate alphabets, probability vectors, transition matrices,
rates, and parameter ranges at construction time. Invalid construction should
raise `ArgumentError`.

## Controls

Every generator accepts an ordered alphabet and emits only values from that
alphabet. Duplicate alphabet entries are rejected because they make empirical
frequency tables ambiguous.

Use `control_capabilities` to ask what a generator claims to control:

```julia
g = make_generator(:MB1c, [:a, :b]; beta = 0.5, d = 500)
control_capabilities(g).marginal
```

Example output:

```julia
:centered
```

Use the empirical helpers for short, transparent diagnostics:

```julia
seq = [:a, :b, :a, :a]
empirical_marginal(seq, [:a, :b])
```

Example output:

```julia
2-element Vector{Float64}:
 0.75
 0.25
```

First-order local structure is represented by `MarkovSpec`:

```julia
spec = MarkovSpec([:a, :b], [0.9 0.1; 0.2 0.8])
local_structure_order(spec)
```

Example output:

```julia
1
```

`MarkovSpec` is currently used by `WaveletMarkov` and `OnOffMarkov`. The
abstract `LocalStructureSpec` is the extension point for a future trigram or
higher-order specification, but SymbolicLongMemorySequences.jl does not yet expose trigram control.

## Output And Provenance

Use `save_sequence` to write generated data with provenance metadata in INC
format:

```julia
g = make_generator(:PB1, [:a, :b]; H = 0.75)
seq = generate(g, 100; rng = StableRNG(11))
save_sequence("seq_pb1.inc", seq, g)
```

The file records the generator type, method identifier, sequence length,
creation date, and generator parameters. Retained validation outputs should use
the same provenance discipline.

## Tests

The main test path checks the public contracts:

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```

The fast tests cover construction errors, output length and element type,
alphabet membership, reproducible RNG use, provenance metadata, factory
discovery, and small bounded checks for declared controls. Broader empirical
studies belong on the [Validation and Benchmarks](@ref) page rather than in the
fast package test suite.
