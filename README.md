# genetic_circuit_mutation_sim
The Genetic Circuit Mutation Simulator predicts how genetic circuits degrade under mutations using probabilistic, stochastic, and biophysical modeling. It provides interactive analysis tools to visualize robustness, mutation paths, sensitivity, and stability landscapes with biologically validated accuracy.
# Genetic Circuit Mutation Simulator

A Python package for Monte Carlo simulation of genetic toggle switch mutations with failure mode classification and robustness analysis.

## Overview

This package simulates how random DNA sequence mutations affect the behavior of a genetic toggle switch circuit. It generates mutant variants, maps sequence changes to parameter perturbations, runs ODE/stochastic simulations, and classifies the resulting phenotypes.

### Key Features

- **Mutation Engine**: Generate controlled mutations (substitutions, insertions, deletions) in annotated DNA regions
- **Sequence-to-Parameter Mapping**: Convert promoter, RBS, CDS, and operator mutations to model parameter changes
- **Simulation**: Deterministic ODE (Tellurium/RoadRunner) and optional stochastic Gillespie SSA
- **Failure Classification**: Classify mutants as bistable, loss-of-bistability, leaky, no-expression, or oscillatory
- **Visualization**: Heatmaps, time courses, and summary figures

## Installation

```bash
# Clone or download the package
cd genetic_circuit_mutation_sim

# Install dependencies
pip install -r requirements.txt
```

### Dependencies

- **Core**: tellurium, libroadrunner, numpy, scipy, pandas
- **Visualization**: matplotlib, seaborn
- **Analysis**: scikit-learn
- **Optional**: SALib (for PRCC sensitivity analysis)

## Quick Start

```bash
# Run Monte Carlo with 1000 variants
python main.py --n 1000

# Run with custom settings
python main.py --n 500 --mutations-per-variant 2 --seed 42 --output-dir ./results
```

### Output Files

- `results.csv`: One row per mutant with sequence, mutations, parameters, and failure label
- `failure_heatmap.png`: Heatmap of failure modes by genomic region
- `robustness_summary.json`: Summary statistics (% bistable, % failed by mode)

## Package Structure

```
genetic_circuit_mutation_sim/
.
├── analysis/           # Probabilistic classifier, evolutionary walks, sensitivity
├── mapping/            # Biophysical mapping (Buffered & Log-Normal)
├── simulation/         # Deterministic ODE & Stochastic SSA engines
├── visualization/      # Static plot generation
├── frontend/           # Modern React/Vite dashboard
├── docs/               # Reproducibility appendix & technical specs
└── main.py             # Feature-rich CLI entry point
```

## Methods

### Toggle Switch Model

The toggle switch is a classic bistable genetic circuit with two mutually repressing genes (A and B). The Antimony model includes:

- **Transcription** with Hill-function repression: `tx_A · Kd^n / (Kd^n + B^n)`
- **Translation**: `tl_A · mRNA_A`
- **Degradation**: First-order decay for mRNA and proteins

Parameters exposed for mutation study:
- `tx_A`, `tx_B`: Transcription rates (AU/time)
- `tl_A`, `tl_B`: Translation rates (AU/time)
- `Kd`: Repressor dissociation constant (AU)
- `n`: Hill coefficient (cooperativity)

### Mutation-to-Parameter Mapping

| Region | Scoring Method | Effect |
|--------|---------------|--------|
| **Promoter** | -35/-10 box PWM match + spacer length | Multiplier for `tx_*` (0.1–1.5) |
| **RBS** | Shine-Dalgarno match + spacing to AUG | Multiplier for `tl_*` (0.1–1.5) |
| **CDS** | Stop codon / frameshift detection | `effect_multiplier` (0 for nonsense) |
| **Operator** | Mismatch count from consensus | Kd multiplier (2^mismatches) |

### Failure Classification

The classifier uses two-start simulations (high-A and high-B initial conditions) to determine bistability:

- **Bistable**: Different stable states from different starts
- **Loss of bistability**: Same final state regardless of start
- **Leaky**: High basal expression when repressed (off/on ratio > 0.2)
- **No expression**: All concentrations near zero
- **Oscillatory**: Sustained oscillations detected (CV > 0.3)

### Parameter Bounds

| Parameter | Range | Unit |
|-----------|-------|------|
| tx_A, tx_B | 1e-4 – 1e2 | AU/time |
| tl_A, tl_B | 1e-4 – 1e2 | AU/time |
| Kd | 1e-3 – 1e3 | AU |
| n | 1 – 6 | dimensionless |

## CLI Options

```
usage: main.py [-h] [--n N] [--mutations-per-variant M] [--stochastic]
               [--stochastic-runs R] [--output-dir DIR] [--seed SEED]
               [--substitution-rate RATE] [--insertion-rate RATE]
               [--deletion-rate RATE] [--quiet]

Options:
  --n N                   Number of mutant variants (default: 1000)
  --mutations-per-variant Number of mutations per variant (default: 1)
  --stochastic           Use Gillespie SSA instead of ODE
  --stochastic-runs R    Stochastic runs per variant (default: 100)
  --output-dir DIR       Output directory (default: .)
  --seed SEED            Random seed for reproducibility
  --quiet                Suppress progress output
```

## Running Tests

```bash
cd genetic_circuit_mutation_sim
python -m pytest tests/test_smoke.py -v
```

## Example Results

A typical run with 1000 variants might show:

```
ROBUSTNESS SUMMARY
==================
Total variants analyzed: 1000
Bistable (functional): 62.3%
Failed: 37.7%
Robustness score: 0.623

Failure mode distribution:
  bistable: 62.3%
  loss_of_bistability: 24.1%
  leaky: 8.5%
  no_expression: 4.2%
  oscillatory: 0.9%
```

## Biological Assumptions

1. **Promoter scoring**: Based on E. coli σ70 consensus sequences. Real promoters have ~100-fold strength variation.

2. **RBS scoring**: Shine-Dalgarno complementarity correlates with translation initiation rate in prokaryotes.

3. **CDS effects**: Nonsense mutations (premature stops) eliminate protein function. Frameshifts similarly disrupt the reading frame.

4. **Operator binding**: Each mismatch roughly doubles Kd (halves binding affinity).

5. **Parameter noise**: Small Gaussian noise (σ=0.05) models natural variation and measurement uncertainty.

## License

MIT License

## Citation

If you use this package in your research, please cite:

```
Genetic Circuit Mutation Simulator
A Monte Carlo framework for analyzing synthetic biology circuit robustness
```
