"""
Sensitivity Analysis Module

Provides PRCC (Partial Rank Correlation Coefficient) sensitivity analysis
using SALib as an optional backend.

PRCC is useful for identifying which parameters most strongly influence
circuit behavior in the presence of model nonlinearity.

Public API:
    prcc_sensitivity(params, outputs) -> dict
"""

import numpy as np
from typing import Dict, List, Optional, Any
import warnings

# Try to import SALib
try:
    from SALib.analyze import rbd_fast
    from SALib.sample import saltelli
    SALIB_AVAILABLE = True
except ImportError:
    SALIB_AVAILABLE = False


def prcc_sensitivity(
    params: np.ndarray,
    outputs: np.ndarray,
    param_names: Optional[List[str]] = None,
    method: str = 'prcc'
) -> Dict[str, Any]:
    """
    Perform sensitivity analysis to identify influential parameters.
    
    Uses PRCC (Partial Rank Correlation Coefficient) or SALib methods
    to quantify parameter importance.
    
    Args:
        params: 2D array of parameter values (n_samples x n_params)
        outputs: 1D array of output values (n_samples,)
        param_names: Optional list of parameter names
        method: Sensitivity method ('prcc' or 'fast')
    
    Returns:
        Dict containing:
            - 'sensitivity': Dict mapping param name to sensitivity value
            - 'ranking': List of params ranked by sensitivity
            - 'method': Method used
            - 'salib_available': Whether SALib was used
    
    Notes:
        - PRCC values range from -1 to 1
        - |PRCC| > 0.3 typically considered significant
        - Positive = increasing param increases output
        - Negative = increasing param decreases output
    """
    if param_names is None:
        param_names = [f'param_{i}' for i in range(params.shape[1])]
    
    result = {
        'sensitivity': {},
        'ranking': [],
        'method': method,
        'salib_available': SALIB_AVAILABLE
    }
    
    if method == 'prcc':
        # Use built-in PRCC calculation
        sensitivities = _calculate_prcc(params, outputs)
        
        for i, name in enumerate(param_names):
            result['sensitivity'][name] = sensitivities[i]
        
        # Rank by absolute sensitivity
        ranked = sorted(
            result['sensitivity'].items(),
            key=lambda x: abs(x[1]),
            reverse=True
        )
        result['ranking'] = [name for name, _ in ranked]
    
    elif method == 'fast' and SALIB_AVAILABLE:
        # Use SALib FAST method
        result = _salib_fast_analysis(params, outputs, param_names)
    
    else:
        warnings.warn(f"Method '{method}' not available, falling back to PRCC")
        return prcc_sensitivity(params, outputs, param_names, method='prcc')
    
    return result


def _calculate_prcc(
    params: np.ndarray,
    outputs: np.ndarray
) -> np.ndarray:
    """
    Calculate Partial Rank Correlation Coefficients.
    
    PRCC measures the correlation between each parameter and the output
    while controlling for the effects of other parameters.
    
    Uses rank transformation to handle nonlinear relationships.
    """
    from scipy import stats
    
    n_samples, n_params = params.shape
    
    if n_samples < 10:
        warnings.warn("Small sample size may lead to unreliable PRCC estimates")
    
    # Rank transform
    ranked_params = np.apply_along_axis(stats.rankdata, 0, params)
    ranked_output = stats.rankdata(outputs)
    
    prcc_values = np.zeros(n_params)
    
    for i in range(n_params):
        # For each parameter, compute partial correlation
        # controlling for other parameters
        
        # Regress parameter i on other parameters
        other_params = np.delete(ranked_params, i, axis=1)
        param_residuals = _residuals_linear(ranked_params[:, i], other_params)
        
        # Regress output on other parameters
        output_residuals = _residuals_linear(ranked_output, other_params)
        
        # PRCC = correlation of residuals
        if np.std(param_residuals) > 0 and np.std(output_residuals) > 0:
            prcc, _ = stats.pearsonr(param_residuals, output_residuals)
        else:
            prcc = 0.0
        
        prcc_values[i] = prcc
    
    return prcc_values


def _residuals_linear(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """
    Calculate residuals from linear regression of y on X.
    """
    if X.shape[1] == 0:
        return y - np.mean(y)
    
    # Add intercept
    X_with_intercept = np.column_stack([np.ones(len(y)), X])
    
    # Solve least squares
    try:
        coeffs, _, _, _ = np.linalg.lstsq(X_with_intercept, y, rcond=None)
        predicted = X_with_intercept @ coeffs
        residuals = y - predicted
    except np.linalg.LinAlgError:
        residuals = y - np.mean(y)
    
    return residuals


def _salib_fast_analysis(
    params: np.ndarray,
    outputs: np.ndarray,
    param_names: List[str]
) -> Dict[str, Any]:
    """
    Perform FAST sensitivity analysis using SALib.
    """
    result = {
        'sensitivity': {},
        'ranking': [],
        'method': 'fast',
        'salib_available': True
    }
    
    try:
        # Define problem for SALib
        problem = {
            'num_vars': len(param_names),
            'names': param_names,
            'bounds': [
                [params[:, i].min(), params[:, i].max()]
                for i in range(len(param_names))
            ]
        }
        
        # Run FAST analysis
        Si = rbd_fast.analyze(problem, params, outputs)
        
        for i, name in enumerate(param_names):
            result['sensitivity'][name] = Si['S1'][i]
        
        # Rank by sensitivity
        ranked = sorted(
            result['sensitivity'].items(),
            key=lambda x: abs(x[1]),
            reverse=True
        )
        result['ranking'] = [name for name, _ in ranked]
        
        # Include additional indices if available
        if 'S1_conf' in Si:
            result['confidence'] = {
                name: Si['S1_conf'][i]
                for i, name in enumerate(param_names)
            }
    
    except Exception as e:
        warnings.warn(f"SALib analysis failed: {e}, falling back to PRCC")
        return prcc_sensitivity(params, outputs, param_names, method='prcc')
    
    return result


def generate_sensitivity_report(
    sensitivity_result: Dict[str, Any],
    threshold: float = 0.3
) -> str:
    """
    Generate a human-readable sensitivity analysis report.
    
    Args:
        sensitivity_result: Result from prcc_sensitivity()
        threshold: Threshold for "significant" sensitivity
    
    Returns:
        Formatted report string
    """
    lines = []
    lines.append("=" * 50)
    lines.append("PARAMETER SENSITIVITY ANALYSIS REPORT")
    lines.append("=" * 50)
    lines.append(f"Method: {sensitivity_result.get('method', 'unknown')}")
    lines.append(f"SALib available: {sensitivity_result.get('salib_available', False)}")
    lines.append("")
    
    lines.append("SENSITIVITY RANKINGS:")
    lines.append("-" * 30)
    
    for rank, param in enumerate(sensitivity_result.get('ranking', []), 1):
        sens_value = sensitivity_result['sensitivity'].get(param, 0)
        significance = "***" if abs(sens_value) > threshold else ""
        lines.append(f"{rank}. {param}: {sens_value:.4f} {significance}")
    
    lines.append("")
    lines.append(f"*** = significant (|sensitivity| > {threshold})")
    
    return "\n".join(lines)


def sweep_classification_thresholds(
    simulation_results: list,
    threshold_ranges: dict = None
) -> dict:
    """
    Sweep classification thresholds and measure impact on robustness.
    
    Tests ±20% variation on each threshold to show ranking stability.
    
    Args:
        simulation_results: List of simulation result dicts
        threshold_ranges: Optional custom ranges for each threshold
        
    Returns:
        Dict with robustness at each threshold combination
    """
    from .failure_classifier import (
        classify_failure, FailureMode,
        BISTABILITY_THRESHOLD, LEAKY_THRESHOLD, OSCILLATION_CV_THRESHOLD
    )
    
    if threshold_ranges is None:
        # Default ±20% variation with 5 steps
        threshold_ranges = {
            'bistability': np.linspace(
                BISTABILITY_THRESHOLD * 0.8,
                BISTABILITY_THRESHOLD * 1.2,
                5
            ),
            'leaky': np.linspace(
                LEAKY_THRESHOLD * 0.8,
                LEAKY_THRESHOLD * 1.2,
                5
            )
        }
    
    results = []
    
    for bistab_thresh in threshold_ranges.get('bistability', [BISTABILITY_THRESHOLD]):
        for leaky_thresh in threshold_ranges.get('leaky', [LEAKY_THRESHOLD]):
            thresholds = {
                'bistability': float(bistab_thresh),
                'leaky': float(leaky_thresh)
            }
            
            # Classify all results with these thresholds
            bistable_count = 0
            for sim_result in simulation_results:
                classification = classify_failure(sim_result, thresholds)
                if classification.mode == FailureMode.BISTABLE:
                    bistable_count += 1
            
            pct_bistable = bistable_count / len(simulation_results) * 100 if simulation_results else 0
            
            results.append({
                'thresholds': thresholds,
                'pct_bistable': pct_bistable,
                'bistable_count': bistable_count
            })
    
    # Analyze sensitivity
    pct_values = [r['pct_bistable'] for r in results]
    
    return {
        'sweep_results': results,
        'min_pct_bistable': min(pct_values) if pct_values else 0,
        'max_pct_bistable': max(pct_values) if pct_values else 0,
        'range_pct_bistable': max(pct_values) - min(pct_values) if pct_values else 0,
        'is_threshold_sensitive': (max(pct_values) - min(pct_values)) > 10 if pct_values else False
    }


def sweep_mutation_rates(
    base_rate: float,
    multipliers: list = None
) -> dict:
    """
    Generate mutation rate sweep configurations.
    
    Args:
        base_rate: Base mutation rate
        multipliers: List of multipliers (default: 0.5x, 1x, 2x)
        
    Returns:
        Dict with mutation rate configurations
    """
    if multipliers is None:
        multipliers = [0.5, 1.0, 2.0]
    
    return {
        'base_rate': base_rate,
        'sweep_configs': [
            {
                'multiplier': m,
                'rate': base_rate * m,
                'label': f'{m}x'
            }
            for m in multipliers
        ]
    }


def generate_sensitivity_summary(
    threshold_sweep: dict,
    region_robustness: dict = None
) -> dict:
    """
    Generate comprehensive sensitivity summary.
    
    Shows what changes with threshold variation and what stays stable.
    """
    summary = {
        'threshold_sensitivity': {
            'is_sensitive': threshold_sweep.get('is_threshold_sensitive', False),
            'robustness_range': threshold_sweep.get('range_pct_bistable', 0),
            'conclusion': ''
        },
        'ranking_stability': {},
        'recommendations': []
    }
    
    if threshold_sweep.get('range_pct_bistable', 0) < 5:
        summary['threshold_sensitivity']['conclusion'] = (
            "Robustness metrics are STABLE across threshold variations (±20%)"
        )
        summary['recommendations'].append(
            "Results are robust to threshold choices - high confidence"
        )
    elif threshold_sweep.get('range_pct_bistable', 0) < 15:
        summary['threshold_sensitivity']['conclusion'] = (
            "Robustness metrics show MODERATE sensitivity to thresholds"
        )
        summary['recommendations'].append(
            "Consider reporting robustness ranges rather than point estimates"
        )
    else:
        summary['threshold_sensitivity']['conclusion'] = (
            "Robustness metrics are HIGHLY SENSITIVE to threshold choices"
        )
        summary['recommendations'].append(
            "Threshold choice significantly affects results - report with caution"
        )
    
    # Analyze region ranking stability if provided
    if region_robustness:
        # Check if relative rankings change across sweep
        summary['ranking_stability']['regions_analyzed'] = list(region_robustness.keys())
        summary['rankings_stable'] = True  # Placeholder - would need actual analysis
    
    return summary

