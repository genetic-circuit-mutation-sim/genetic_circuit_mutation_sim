"""
Failure Classifier Module

Classifies toggle switch simulation results into functional phenotypes.
This is crucial for understanding how mutations affect circuit behavior.

Failure modes:
- BISTABLE: Normal function - two stable states depending on initial condition
- LOSS_OF_BISTABILITY: Circuit locks into one state regardless of start
- LEAKY: High basal expression in "off" state
- NO_EXPRESSION: Both proteins at very low levels
- OSCILLATORY: Sustained oscillations (unexpected for toggle switch)

Public API:
    classify_failure(time_series) -> FailureMode
"""

import numpy as np
from enum import Enum
from typing import Dict, Any, Optional, Tuple
from dataclasses import dataclass


class FailureMode(Enum):
    """
    Enumeration of possible toggle switch phenotypes.
    
    BISTABLE: Normal bistable behavior (two distinct stable states)
    LOSS_OF_BISTABILITY: Monostable - converges to same state from any start
    LEAKY: High basal expression when gene should be repressed
    NO_EXPRESSION: Both genes have near-zero expression
    OSCILLATORY: Sustained oscillations detected
    SIMULATION_FAILED: Simulation did not complete successfully
    """
    BISTABLE = "bistable"
    LOSS_OF_BISTABILITY = "loss_of_bistability"
    LEAKY = "leaky"
    NO_EXPRESSION = "no_expression"
    OSCILLATORY = "oscillatory"
    SIMULATION_FAILED = "simulation_failed"


@dataclass
class ClassificationResult:
    """
    Result of failure classification.
    
    Attributes:
        mode: The classified failure mode
        confidence: Confidence score (0-1)
        details: Dict with classification details
    """
    mode: FailureMode
    confidence: float
    details: Dict[str, Any]


# Classification thresholds (biologically motivated)
BISTABILITY_THRESHOLD = 0.5  # Minimum ratio difference between states
LOW_EXPRESSION_THRESHOLD = 0.1  # Below this = no expression
HIGH_EXPRESSION_THRESHOLD = 10.0  # Reference for "high" expression
LEAKY_THRESHOLD = 0.5  # Ratio of off/on above this = leaky
OSCILLATION_CV_THRESHOLD = 0.3  # Coefficient of variation for oscillations
STEADY_STATE_FRACTION = 0.2  # Fraction of trajectory to use for SS analysis


def classify_failure(
    simulation_result: Dict[str, Any],
    thresholds: Optional[Dict[str, float]] = None
) -> ClassificationResult:
    """
    Classify the phenotype of a toggle switch simulation result.
    
    Uses the two-start protocol results to determine if the circuit
    is bistable, monostable, leaky, or has other failure modes.
    
    Args:
        simulation_result: Dict from simulate_deterministic() containing:
            - 'steady_state_A_highA': Final A from high-A start
            - 'steady_state_B_highA': Final B from high-A start
            - 'steady_state_A_highB': Final A from high-B start
            - 'steady_state_B_highB': Final B from high-B start
            - 'result_highA': Full time series from high-A start
            - 'result_highB': Full time series from high-B start
            - 'simulation_successful': Boolean
        thresholds: Optional custom thresholds dict
    
    Returns:
        ClassificationResult with mode, confidence, and details
    
    Classification logic:
        1. Check for simulation failure
        2. Check for oscillatory behavior
        3. Check for no expression
        4. Check for bistability (different final states)
        5. Check for leakiness
        6. Default to loss_of_bistability
    """
    # Default thresholds
    thresh = {
        'bistability': BISTABILITY_THRESHOLD,
        'low_expression': LOW_EXPRESSION_THRESHOLD,
        'leaky': LEAKY_THRESHOLD,
        'oscillation_cv': OSCILLATION_CV_THRESHOLD
    }
    if thresholds:
        thresh.update(thresholds)
    
    details = {}
    
    # Check simulation success
    if not simulation_result.get('simulation_successful', False):
        return ClassificationResult(
            mode=FailureMode.SIMULATION_FAILED,
            confidence=1.0,
            details={'error': simulation_result.get('error_message', 'Unknown error')}
        )
    
    # Extract steady-state values
    ss_A_highA = simulation_result.get('steady_state_A_highA', 0)
    ss_B_highA = simulation_result.get('steady_state_B_highA', 0)
    ss_A_highB = simulation_result.get('steady_state_A_highB', 0)
    ss_B_highB = simulation_result.get('steady_state_B_highB', 0)
    
    details['steady_states'] = {
        'A_from_highA': ss_A_highA,
        'B_from_highA': ss_B_highA,
        'A_from_highB': ss_A_highB,
        'B_from_highB': ss_B_highB
    }
    
    # Check for oscillations (using full time series if available)
    oscillatory, osc_details = _check_oscillatory(simulation_result)
    details['oscillation_check'] = osc_details
    
    if oscillatory:
        return ClassificationResult(
            mode=FailureMode.OSCILLATORY,
            confidence=osc_details.get('confidence', 0.8),
            details=details
        )
    
    # Check for no expression
    max_expression = max(ss_A_highA, ss_B_highA, ss_A_highB, ss_B_highB)
    
    if max_expression < thresh['low_expression']:
        return ClassificationResult(
            mode=FailureMode.NO_EXPRESSION,
            confidence=1.0 - max_expression / thresh['low_expression'],
            details=details
        )
    
    # Check for bistability
    # In bistable state: high-A start -> high A, low B
    #                   high-B start -> low A, high B
    bistable, bistable_conf = _check_bistability(
        ss_A_highA, ss_B_highA, ss_A_highB, ss_B_highB,
        thresh['bistability']
    )
    
    details['bistability_check'] = {
        'is_bistable': bistable,
        'confidence': bistable_conf
    }
    
    if bistable:
        # Also check for leakiness in bistable state
        leaky, leaky_details = _check_leakiness(
            ss_A_highA, ss_B_highA, ss_A_highB, ss_B_highB,
            thresh['leaky']
        )
        details['leakiness_check'] = leaky_details
        
        if leaky:
            return ClassificationResult(
                mode=FailureMode.LEAKY,
                confidence=leaky_details.get('confidence', 0.7),
                details=details
            )
        
        return ClassificationResult(
            mode=FailureMode.BISTABLE,
            confidence=bistable_conf,
            details=details
        )
    
    # Check for leakiness in monostable state
    leaky, leaky_details = _check_leakiness(
        ss_A_highA, ss_B_highA, ss_A_highB, ss_B_highB,
        thresh['leaky']
    )
    details['leakiness_check'] = leaky_details
    
    if leaky:
        return ClassificationResult(
            mode=FailureMode.LEAKY,
            confidence=leaky_details.get('confidence', 0.7),
            details=details
        )
    
    # Default: loss of bistability (monostable)
    return ClassificationResult(
        mode=FailureMode.LOSS_OF_BISTABILITY,
        confidence=1.0 - bistable_conf,
        details=details
    )


def _check_bistability(
    ss_A_highA: float,
    ss_B_highA: float,
    ss_A_highB: float,
    ss_B_highB: float,
    threshold: float
) -> Tuple[bool, float]:
    """
    Check if the system exhibits bistability.
    
    Bistability requires:
    - High-A start leads to high A, low B
    - High-B start leads to low A, high B
    - The two states are significantly different
    """
    # Normalize to avoid division by zero
    eps = 1e-10
    
    # Check if high-A start gives high A (relative to B)
    ratio_highA = ss_A_highA / (ss_B_highA + eps)
    
    # Check if high-B start gives high B (relative to A)
    ratio_highB = ss_B_highB / (ss_A_highB + eps)
    
    # Both ratios should be > 1 for proper bistability
    state_separation = min(ratio_highA, ratio_highB)
    
    # Also check absolute difference between states
    A_diff = abs(ss_A_highA - ss_A_highB)
    B_diff = abs(ss_B_highA - ss_B_highB)
    
    max_val = max(ss_A_highA, ss_B_highA, ss_A_highB, ss_B_highB) + eps
    relative_diff = (A_diff + B_diff) / (2 * max_val)
    
    # Bistable if states are well-separated
    is_bistable = (
        state_separation > (1 + threshold) and
        relative_diff > threshold
    )
    
    # Confidence based on separation quality
    confidence = min(1.0, relative_diff / (threshold + 0.1))
    
    return is_bistable, confidence


def _check_oscillatory(
    simulation_result: Dict[str, Any]
) -> Tuple[bool, Dict]:
    """
    Check for oscillatory behavior in time series.
    
    Oscillations are detected by checking coefficient of variation
    in the latter half of the simulation.
    """
    details = {'oscillatory': False, 'confidence': 0.0}
    
    # Get time series data
    result_highA = simulation_result.get('result_highA', {})
    A_series = result_highA.get('A') if isinstance(result_highA, dict) else None
    B_series = result_highA.get('B') if isinstance(result_highA, dict) else None
    
    if A_series is None or B_series is None:
        return False, details
    
    if len(A_series) < 50:
        return False, details
    
    # Use latter portion for steady-state analysis
    n = len(A_series)
    start_idx = int(n * (1 - STEADY_STATE_FRACTION))
    
    A_late = A_series[start_idx:]
    B_late = B_series[start_idx:]
    
    # Calculate coefficient of variation
    A_cv = np.std(A_late) / (np.mean(A_late) + 1e-10)
    B_cv = np.std(B_late) / (np.mean(B_late) + 1e-10)
    
    max_cv = max(A_cv, B_cv)
    
    details['A_cv'] = A_cv
    details['B_cv'] = B_cv
    
    # High CV suggests oscillations (not reaching steady state)
    if max_cv > OSCILLATION_CV_THRESHOLD:
        # Additional check: look for periodicity
        # Simple approach: check for sign changes in derivative
        A_deriv = np.diff(A_late)
        sign_changes = np.sum(np.diff(np.sign(A_deriv)) != 0)
        
        # Multiple sign changes suggests oscillation
        if sign_changes > 5:
            details['oscillatory'] = True
            details['confidence'] = min(1.0, sign_changes / 20)
            return True, details
    
    return False, details


def _check_leakiness(
    ss_A_highA: float,
    ss_B_highA: float,
    ss_A_highB: float,
    ss_B_highB: float,
    threshold: float
) -> Tuple[bool, Dict]:
    """
    Check for leaky expression (high basal expression in "off" state).
    
    In a good toggle switch:
    - When A is high, B should be very low (strongly repressed)
    - When B is high, A should be very low (strongly repressed)
    """
    details = {'leaky': False, 'confidence': 0.0}
    
    eps = 1e-10
    
    # Calculate off/on ratios
    # In high-A state, B should be low
    B_leakiness_highA = ss_B_highA / (ss_A_highA + eps)
    
    # In high-B state, A should be low
    A_leakiness_highB = ss_A_highB / (ss_B_highB + eps)
    
    max_leakiness = max(B_leakiness_highA, A_leakiness_highB)
    
    details['B_leakiness_in_highA_state'] = B_leakiness_highA
    details['A_leakiness_in_highB_state'] = A_leakiness_highB
    details['max_leakiness'] = max_leakiness
    
    if max_leakiness > threshold:
        details['leaky'] = True
        details['confidence'] = min(1.0, max_leakiness / (threshold * 2))
        return True, details
    
    return False, details


def classify_failure_simple(
    ss_A_highA: float,
    ss_B_highA: float,
    ss_A_highB: float,
    ss_B_highB: float
) -> str:
    """
    Simplified failure classification from steady-state values only.
    
    Quick classification for batch processing.
    
    Returns:
        String label: 'bistable', 'loss_of_bistability', 'leaky', 
                     'no_expression', or 'simulation_failed'
    """
    sim_result = {
        'simulation_successful': True,
        'steady_state_A_highA': ss_A_highA,
        'steady_state_B_highA': ss_B_highA,
        'steady_state_A_highB': ss_A_highB,
        'steady_state_B_highB': ss_B_highB,
        'result_highA': {'A': np.array([ss_A_highA]), 'B': np.array([ss_B_highA])},
        'result_highB': {'A': np.array([ss_A_highB]), 'B': np.array([ss_B_highB])}
    }
    
    result = classify_failure(sim_result)
    return result.mode.value
