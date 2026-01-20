"""
Genetic Circuit Mutation Simulator Package

A Python package for Monte Carlo simulation of genetic toggle switch mutations
with failure mode classification and robustness analysis.

Main modules:
    - mutation_engine: Generate controlled DNA mutations
    - mapping: Map sequence changes to parameter perturbations
    - simulation: ODE/stochastic simulation engines
    - analysis: Failure classification and robustness metrics
    - visualization: Plotting functions
"""

__version__ = '1.0.0'
__author__ = 'Genetic Circuit Mutation Simulator'

# Convenience imports
from mutation_engine import mutate_sequence, AnnotatedRegion, MutationType
from mapping import promoter_to_tx_rate, rbs_to_tl_rate, cds_to_protein_effect
from analysis import classify_failure, FailureMode, compute_robustness
