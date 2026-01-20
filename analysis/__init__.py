"""
Analysis Module

Provides tools for analyzing toggle switch simulation results,
including failure classification, robustness metrics, sensitivity analysis,
probabilistic classification, evolutionary walks, and validation.

Public API:
    - classify_failure: Classify simulation result as failure mode
    - compute_robustness: Compute robustness metrics from batch results
    - prcc_sensitivity: PRCC-based sensitivity analysis (SALib wrapper)
    - classify_probabilistic: Ensemble-based probabilistic classification
    - PhenotypeDistribution: Probability distribution over phenotypes
    - simulate_mutation_walk: Multi-mutation trajectory simulation
    - run_validation_suite: Sanity checks against known biology
"""

from .failure_classifier import classify_failure, FailureMode
from .robustness_metrics import compute_robustness
from .sensitivity import prcc_sensitivity
from .probabilistic_classifier import (
    classify_probabilistic,
    PhenotypeDistribution,
    compute_ensemble_robustness,
    classify_with_threshold_sweep
)
from .evolutionary_walks import (
    simulate_mutation_walk,
    simulate_multiple_walks,
    compute_decay_curve,
    MutationTrajectory,
    DecayCurve
)
from .validation import (
    run_validation_suite,
    ValidationReport,
    ValidationResult,
    ValidationStatus
)

__all__ = [
    'classify_failure',
    'FailureMode',
    'compute_robustness',
    'prcc_sensitivity',
    # Probabilistic classification
    'classify_probabilistic',
    'PhenotypeDistribution',
    'compute_ensemble_robustness',
    'classify_with_threshold_sweep',
    # Evolutionary walks
    'simulate_mutation_walk',
    'simulate_multiple_walks',
    'compute_decay_curve',
    'MutationTrajectory',
    'DecayCurve',
    # Validation
    'run_validation_suite',
    'ValidationReport',
    'ValidationResult',
    'ValidationStatus'
]
