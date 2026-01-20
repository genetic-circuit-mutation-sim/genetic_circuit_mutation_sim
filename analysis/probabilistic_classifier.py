"""
Probabilistic Phenotype Classifier Module

Replaces hard-threshold binary classification with probability distributions
based on ensemble simulations. This provides:
- Continuous robustness quantification
- Reduced threshold sensitivity
- More biologically realistic behavior modeling

Public API:
    classify_probabilistic(simulation_fn, params, n_runs) -> PhenotypeDistribution
    compute_ensemble_robustness(results) -> float
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, Callable, List, Tuple
from enum import Enum
import warnings

from .failure_classifier import (
    FailureMode, 
    classify_failure, 
    ClassificationResult,
    BISTABILITY_THRESHOLD,
    OSCILLATION_CV_THRESHOLD,
    LEAKY_THRESHOLD
)


@dataclass
class PhenotypeDistribution:
    """
    Probability distribution over phenotype outcomes.
    
    Attributes:
        bistable_prob: Probability of bistable behavior
        loss_of_bistability_prob: Probability of monostable behavior
        leaky_prob: Probability of leaky expression
        no_expression_prob: Probability of no expression
        oscillatory_prob: Probability of oscillatory behavior
        simulation_failed_prob: Probability of simulation failure
        
        n_runs: Number of ensemble runs performed
        confidence_interval: 95% CI for bistable probability
        robustness_score: Continuous robustness (0.0 to 1.0)
        details: Additional statistics from ensemble
    """
    bistable_prob: float = 0.0
    loss_of_bistability_prob: float = 0.0
    leaky_prob: float = 0.0
    no_expression_prob: float = 0.0
    oscillatory_prob: float = 0.0
    simulation_failed_prob: float = 0.0
    
    n_runs: int = 0
    confidence_interval: Tuple[float, float] = (0.0, 0.0)
    robustness_score: float = 0.0
    details: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            'bistable_prob': self.bistable_prob,
            'loss_of_bistability_prob': self.loss_of_bistability_prob,
            'leaky_prob': self.leaky_prob,
            'no_expression_prob': self.no_expression_prob,
            'oscillatory_prob': self.oscillatory_prob,
            'simulation_failed_prob': self.simulation_failed_prob,
            'n_runs': self.n_runs,
            'confidence_interval': list(self.confidence_interval),
            'robustness_score': self.robustness_score,
            'details': self.details
        }
    
    @property
    def dominant_phenotype(self) -> FailureMode:
        """Return the most likely phenotype."""
        probs = {
            FailureMode.BISTABLE: self.bistable_prob,
            FailureMode.LOSS_OF_BISTABILITY: self.loss_of_bistability_prob,
            FailureMode.LEAKY: self.leaky_prob,
            FailureMode.NO_EXPRESSION: self.no_expression_prob,
            FailureMode.OSCILLATORY: self.oscillatory_prob,
            FailureMode.SIMULATION_FAILED: self.simulation_failed_prob
        }
        return max(probs, key=probs.get)
    
    @property 
    def is_functional(self) -> bool:
        """Return True if bistable probability > 50%."""
        return self.bistable_prob > 0.5


def classify_probabilistic(
    simulation_fn: Callable[[Dict[str, float]], Dict[str, Any]],
    base_params: Dict[str, float],
    n_runs: int = 20,
    noise_sigma: float = 0.1,
    noise_type: str = 'parameter',
    seed: Optional[int] = None
) -> PhenotypeDistribution:
    """
    Perform probabilistic phenotype classification using ensemble simulations.
    
    Runs multiple simulations with noise variations and aggregates results
    into probability distributions for each phenotype.
    
    Args:
        simulation_fn: Function that takes parameters dict and returns 
                      simulation result dict (compatible with classify_failure)
        base_params: Base parameter values to perturb
        n_runs: Number of ensemble runs (default: 20)
        noise_sigma: Relative noise standard deviation (default: 0.1 = 10%)
        noise_type: Type of noise to apply:
            - 'parameter': Gaussian noise on parameter values
            - 'initial': Variation in initial conditions
            - 'both': Both parameter and initial condition noise
        seed: Random seed for reproducibility
        
    Returns:
        PhenotypeDistribution with probability for each phenotype
        
    Notes:
        - More runs = more accurate probability estimates
        - 20 runs gives ~10% precision on probabilities
        - 100 runs gives ~5% precision
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Track phenotype counts
    phenotype_counts = {mode: 0 for mode in FailureMode}
    classification_results: List[ClassificationResult] = []
    
    # Store steady-state values for variance analysis
    steady_states = {
        'A_highA': [], 'B_highA': [],
        'A_highB': [], 'B_highB': []
    }
    
    for run_idx in range(n_runs):
        # Apply noise to parameters
        noisy_params = _apply_parameter_noise(
            base_params, noise_sigma, noise_type
        )
        
        try:
            # Run simulation
            sim_result = simulation_fn(noisy_params)
            
            # Classify phenotype
            classification = classify_failure(sim_result)
            classification_results.append(classification)
            phenotype_counts[classification.mode] += 1
            
            # Store steady-states if available
            if sim_result.get('simulation_successful', False):
                steady_states['A_highA'].append(
                    sim_result.get('steady_state_A_highA', 0)
                )
                steady_states['B_highA'].append(
                    sim_result.get('steady_state_B_highA', 0)
                )
                steady_states['A_highB'].append(
                    sim_result.get('steady_state_A_highB', 0)
                )
                steady_states['B_highB'].append(
                    sim_result.get('steady_state_B_highB', 0)
                )
                
        except Exception as e:
            phenotype_counts[FailureMode.SIMULATION_FAILED] += 1
            warnings.warn(f"Simulation failed in run {run_idx}: {e}")
    
    # Calculate probabilities
    probs = {mode: count / n_runs for mode, count in phenotype_counts.items()}
    
    # Calculate confidence interval for bistable probability using Wilson score
    bistable_ci = _wilson_confidence_interval(
        phenotype_counts[FailureMode.BISTABLE], n_runs
    )
    
    # Calculate robustness score (weighted by confidence)
    robustness = _compute_weighted_robustness(probs, classification_results)
    
    # Compute additional statistics
    details = _compute_ensemble_statistics(steady_states, classification_results)
    
    return PhenotypeDistribution(
        bistable_prob=probs[FailureMode.BISTABLE],
        loss_of_bistability_prob=probs[FailureMode.LOSS_OF_BISTABILITY],
        leaky_prob=probs[FailureMode.LEAKY],
        no_expression_prob=probs[FailureMode.NO_EXPRESSION],
        oscillatory_prob=probs[FailureMode.OSCILLATORY],
        simulation_failed_prob=probs[FailureMode.SIMULATION_FAILED],
        n_runs=n_runs,
        confidence_interval=bistable_ci,
        robustness_score=robustness,
        details=details
    )


def _apply_parameter_noise(
    params: Dict[str, float],
    sigma: float,
    noise_type: str
) -> Dict[str, float]:
    """
    Apply log-normal noise to parameters.
    
    Uses log-normal distribution to ensure parameters stay positive
    and changes are multiplicative (biologically realistic).
    """
    noisy_params = params.copy()
    
    if noise_type in ('parameter', 'both'):
        for key, value in params.items():
            if key in ('tx_A', 'tx_B', 'tl_A', 'tl_B', 'Kd'):
                # Log-normal multiplicative noise
                multiplier = np.exp(np.random.normal(0, sigma))
                noisy_params[key] = value * multiplier
            elif key == 'n':
                # Hill coefficient: additive bounded noise
                noise = np.random.normal(0, sigma * value)
                noisy_params[key] = np.clip(value + noise, 1.0, 6.0)
    
    return noisy_params


def _wilson_confidence_interval(
    successes: int, 
    total: int, 
    confidence: float = 0.95
) -> Tuple[float, float]:
    """
    Calculate Wilson score confidence interval for a proportion.
    
    More accurate than normal approximation for small sample sizes.
    """
    if total == 0:
        return (0.0, 1.0)
    
    from scipy import stats
    
    p = successes / total
    z = stats.norm.ppf(1 - (1 - confidence) / 2)
    
    denominator = 1 + z**2 / total
    center = (p + z**2 / (2 * total)) / denominator
    margin = z * np.sqrt(p * (1 - p) / total + z**2 / (4 * total**2)) / denominator
    
    return (max(0, center - margin), min(1, center + margin))


def _compute_weighted_robustness(
    probs: Dict[FailureMode, float],
    results: List[ClassificationResult]
) -> float:
    """
    Compute weighted robustness score.
    
    Weights bistable probability by average confidence of classifications.
    """
    if not results:
        return 0.0
    
    # Base robustness = bistable probability
    base_robustness = probs[FailureMode.BISTABLE]
    
    # Weight by average classification confidence
    bistable_results = [r for r in results if r.mode == FailureMode.BISTABLE]
    if bistable_results:
        avg_confidence = np.mean([r.confidence for r in bistable_results])
        weighted_robustness = base_robustness * avg_confidence
    else:
        weighted_robustness = 0.0
    
    return float(weighted_robustness)


def _compute_ensemble_statistics(
    steady_states: Dict[str, List[float]],
    results: List[ClassificationResult]
) -> Dict[str, Any]:
    """
    Compute additional statistics from ensemble results.
    """
    stats = {}
    
    # Steady-state variability
    for key, values in steady_states.items():
        if values:
            stats[f'{key}_mean'] = float(np.mean(values))
            stats[f'{key}_std'] = float(np.std(values))
            stats[f'{key}_cv'] = float(np.std(values) / (np.mean(values) + 1e-10))
    
    # Mode transition analysis (how often does phenotype change?)
    if len(results) > 1:
        modes = [r.mode for r in results]
        transitions = sum(1 for i in range(1, len(modes)) if modes[i] != modes[i-1])
        stats['phenotype_variability'] = transitions / (len(modes) - 1)
    
    # Average confidence across all classifications
    if results:
        stats['avg_confidence'] = float(np.mean([r.confidence for r in results]))
    
    return stats


def compute_ensemble_robustness(
    distributions: List[PhenotypeDistribution]
) -> Dict[str, Any]:
    """
    Compute aggregate robustness metrics from multiple ensemble results.
    
    Args:
        distributions: List of PhenotypeDistribution from multiple mutants
        
    Returns:
        Dict with aggregate robustness metrics
    """
    if not distributions:
        return {
            'total_mutants': 0,
            'mean_robustness': 0.0,
            'median_robustness': 0.0,
            'std_robustness': 0.0,
            'pct_functional': 0.0,
            'pct_high_confidence_functional': 0.0
        }
    
    robustness_scores = [d.robustness_score for d in distributions]
    bistable_probs = [d.bistable_prob for d in distributions]
    
    # Count high-confidence functional (bistable_prob > 0.8 and CI doesn't include 0.5)
    high_conf_functional = sum(
        1 for d in distributions 
        if d.bistable_prob > 0.8 and d.confidence_interval[0] > 0.5
    )
    
    return {
        'total_mutants': len(distributions),
        'mean_robustness': float(np.mean(robustness_scores)),
        'median_robustness': float(np.median(robustness_scores)),
        'std_robustness': float(np.std(robustness_scores)),
        'pct_functional': float(np.mean([d.is_functional for d in distributions]) * 100),
        'pct_high_confidence_functional': float(high_conf_functional / len(distributions) * 100),
        'robustness_distribution': {
            'q25': float(np.percentile(robustness_scores, 25)),
            'q50': float(np.percentile(robustness_scores, 50)),
            'q75': float(np.percentile(robustness_scores, 75)),
            'q90': float(np.percentile(robustness_scores, 90))
        }
    }


def classify_with_threshold_sweep(
    simulation_result: Dict[str, Any],
    threshold_variations: Optional[Dict[str, List[float]]] = None
) -> Dict[str, Any]:
    """
    Classify phenotype across a range of threshold values.
    
    Useful for sensitivity analysis - shows how classification
    changes with threshold adjustments.
    
    Args:
        simulation_result: Single simulation result
        threshold_variations: Dict mapping threshold names to list of values to test
            Default: ±20% variation on each threshold
            
    Returns:
        Dict with classification at each threshold combination
    """
    if threshold_variations is None:
        # Default ±20% variation
        threshold_variations = {
            'bistability': [
                BISTABILITY_THRESHOLD * 0.8,
                BISTABILITY_THRESHOLD,
                BISTABILITY_THRESHOLD * 1.2
            ],
            'oscillation_cv': [
                OSCILLATION_CV_THRESHOLD * 0.8,
                OSCILLATION_CV_THRESHOLD,
                OSCILLATION_CV_THRESHOLD * 1.2
            ],
            'leaky': [
                LEAKY_THRESHOLD * 0.8,
                LEAKY_THRESHOLD,
                LEAKY_THRESHOLD * 1.2
            ]
        }
    
    results = []
    
    for bistab in threshold_variations.get('bistability', [BISTABILITY_THRESHOLD]):
        for osc_cv in threshold_variations.get('oscillation_cv', [OSCILLATION_CV_THRESHOLD]):
            for leaky in threshold_variations.get('leaky', [LEAKY_THRESHOLD]):
                thresholds = {
                    'bistability': bistab,
                    'oscillation_cv': osc_cv,
                    'leaky': leaky
                }
                
                classification = classify_failure(simulation_result, thresholds)
                
                results.append({
                    'thresholds': thresholds,
                    'phenotype': classification.mode.value,
                    'confidence': classification.confidence
                })
    
    # Analyze threshold sensitivity
    phenotypes = [r['phenotype'] for r in results]
    unique_phenotypes = set(phenotypes)
    
    return {
        'classifications': results,
        'is_threshold_sensitive': len(unique_phenotypes) > 1,
        'phenotype_count': {p: phenotypes.count(p) for p in unique_phenotypes},
        'dominant_phenotype': max(unique_phenotypes, key=phenotypes.count)
    }
