# Reference

## Abstract type

```@docs
SymbolicLongMemorySequences
LRDGenerator
```

## Common interface

```@docs
generate
generate_with_latent
MethodInfo
ParameterInfo
method_ids
method_info
method_parameters
make_generator
save_sequence
LocalStructureSpec
MarkovSpec
local_structure_order
ControlCapabilities
control_capabilities
```

## Property-based generators

```@docs
LatentSource
Symbolizer
PropertyBasedGenerator
SpectralFGNSource
HaarLRDSource
IntermittentMapSource
QuantileSymbolizer
ArgmaxSymbolizer
MarkovRegimeSymbolizer
latent_width
generate_latent
symbolize
SpectralFGN
LGCM
WaveletMarkov
IntermittentMapSymbols
```

## Model-based generators

```@docs
LAMP
DyadicLAMP
CalibratedAdditiveMarkov
OnOffMarkov
FSS
HawkesSymbol
DuplicationMutation
```

## Utilities

```@docs
target_marginal
empirical_marginal
empirical_bigram
empirical_trigram
bin_counts
total_variation
rowwise_total_variation
validate_transition_matrix
stationary_distribution
lamp_repeat_transition
SymbolicLongMemorySequences.quantize_to_symbols
SymbolicLongMemorySequences._fgn_spectral
SymbolicLongMemorySequences._pareto_sample
```
