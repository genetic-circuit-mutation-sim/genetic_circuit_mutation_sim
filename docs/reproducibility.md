# Reproducibility Appendix

This document provides the technical parameters, thresholds, and assumptions required to reproduce the simulation results presented in this project.

## 1. Simulation Parameters

The dynamical system is modeled using the standard toggle switch ODEs (deterministic) or Gillespie SSA (stochastic).

| Parameter | Symbol | Wild-type Value | Search Range | Description |
|-----------|--------|-----------------|--------------|-------------|
| Transcription Rate | $k_{tx}$ | 1.0 | [0.1, 1.5] | Maximum transcription rate |
| Translation Rate | $k_{tl}$ | 1.0 | [0.1, 1.5] | Maximum translation rate |
| Dissociation Const | $K_d$ | 1.0 | [0.1, 10.0] | Repressor-DNA affinity |
| Hill Coefficient | $n$ | 2.0 | [1.0, 4.0] | Cooperativity |
| Degradation (mRNA) | $d_m$ | 1.0 | Fixed | Fixed at 1.0 (relative time) |
| Degradation (Prot) | $d_p$ | 0.1 | Fixed | 10x more stable than mRNA |

## 2. Classification Thresholds

Phenotypes are classified based on steady-state concentrations (arbitrary units).

| Threshold | Value | Logic |
|-----------|-------|-------|
| **Bistability** | 0.5 | $|A^* - B^*| > 0.5$ for $highA$ and $highB$ initial conditions |
| **Leaky Offset** | 0.5 | $A_{off} < 0.5$ AND $B_{off} < 0.5$ |
| **Low Expression** | 0.1 | $A_{max} < 0.1$ OR $B_{max} < 0.1$ (classified as No Expression) |
| **Oscillation CV** | 0.3 | Coefficient of variation in late-stage trace > 0.3 |

## 3. Sequence-to-Parameter Mapping

Mappings use biophysically inspired heuristics with **Log-Normal Scaling**.

### Promoter Mapping
- **Consensus**: TTGACA (-35) and TATAAT (-10)
- **Neutral Zone**: Multipliers in range [0.7, 1.3] are buffered.
- **Noise**: Log-normal noise with $\sigma = 0.05$.

### RBS Mapping
- **Consensus**: AGGAGG (SD)
- **Spacing**: 5-9 bp optimal.
- **Neutral Zone**: [0.7, 1.3]

### CDS Mapping
- **Synonymous**: Explicitly neutral (Multiplier = 1.0).
- **Nonsense**: Truncation < 80% length results in multiplier = 0.0.
- **Missense**: 15% reduction per non-conservative change (buffered for 1-2 changes).

## 4. Random Seeds

For batch simulation reproducibility:
- **Default Seed**: 42
- **Ensemble Count**: 20 runs per variant for probabilistic classification.
- **Stochastic SSA**: Uses `numpy.random` with the global seed specified.

## 5. Model Assumptions

1. **Perfect Mixing**: Assumes cell contents are well-mixed (no spatial gradients).
2. **Infinite Resource**: No competition for RNA polymerase or ribosomes.
3. **Relative Rates**: All rates are normalized to wild-type = 1.0.
4. **Binary States**: "Robustness" is defined as the persistence of the bistable phenotype.
