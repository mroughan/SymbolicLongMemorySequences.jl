# API

S5.jl has two public construction paths:

- use `make_generator` for standard, named cases across all implemented methods;
- use the method-specific constructors when a study needs full control of the
  scientific parameters.

Both paths return an `LRDGenerator` and both use the same generation call:

```julia
seq = generate(g, n; rng)
```

Always pass an explicit `rng` when a sequence needs to be reproducible.

## Standard Method Factory

The factory is the cleanest way to list methods, inspect their standard
parameters, and build a generator for examples, smoke tests, and benchmark grids.

```julia
using S5, StableRNGs

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
method_info(:SpectralFGN).id
```

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

## Explicit Constructors

Use explicit constructors when a model needs exact transition matrices,
excitation matrices, rates, or calibration parameters.

```julia
using S5, StableRNGs

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
higher-order specification, but S5.jl does not yet expose trigram control.

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
