"""
Promoter Mapping Module

Maps promoter DNA sequences to transcription rate multipliers using
heuristic scoring based on sigma70 consensus sequences.

Biological Background:
    Bacterial promoters recognized by sigma70 contain two conserved elements:
    - The -35 box (consensus: TTGACA, ~17 bp upstream of -10)
    - The -10 box (consensus: TATAAT, ~10 bp upstream of transcription start)
    
    Mutations in these boxes reduce RNA polymerase binding affinity,
    thereby reducing transcription rate.

Public API:
    promoter_to_tx_rate(seq, score_params) -> float
"""

import numpy as np
from typing import Dict, Optional, Tuple
import re


# Sigma70 consensus sequences
MINUS_35_CONSENSUS = "TTGACA"
MINUS_10_CONSENSUS = "TATAAT"

# Position weight matrix derived from E. coli promoter studies
# Values represent log-odds of each nucleotide at each position
# Reference: consensus sequence data from various promoter studies
MINUS_35_PWM = {
    0: {'T': 1.0, 'A': -0.5, 'G': -1.0, 'C': -1.0},
    1: {'T': 1.0, 'A': -0.5, 'G': -1.0, 'C': -1.0},
    2: {'G': 1.0, 'A': 0.3, 'T': -0.5, 'C': -1.0},
    3: {'A': 1.0, 'G': 0.2, 'T': -0.5, 'C': -0.5},
    4: {'C': 1.0, 'T': 0.3, 'A': -0.3, 'G': -0.5},
    5: {'A': 1.0, 'G': 0.1, 'T': -0.3, 'C': -0.5},
}

MINUS_10_PWM = {
    0: {'T': 1.0, 'A': 0.2, 'G': -1.0, 'C': -1.0},
    1: {'A': 1.0, 'T': 0.1, 'G': -0.5, 'C': -0.5},
    2: {'T': 1.0, 'A': 0.3, 'G': -0.5, 'C': -0.5},
    3: {'A': 1.0, 'T': 0.3, 'G': -0.3, 'C': -0.5},
    4: {'A': 1.0, 'T': 0.5, 'G': -0.5, 'C': -0.5},
    5: {'T': 1.0, 'A': 0.3, 'G': -0.5, 'C': -0.5},
}

# Optimal spacer length between -35 and -10 is 17 bp
OPTIMAL_SPACER = 17
SPACER_TOLERANCE = 2  # Allow ±2 bp variation

# Neutral zone for buffered parameter mapping
# Mutations within this zone have minimal effect on transcription
NEUTRAL_ZONE = (0.7, 1.3)  # Multipliers within this range are biologically tolerable

# Mutation severity profiles
MUTATION_SEVERITY = {
    'conservative': {'scale': 0.5, 'neutral_width': 0.4},   # Smaller effects
    'standard': {'scale': 1.0, 'neutral_width': 0.3},        # Default
    'aggressive': {'scale': 1.5, 'neutral_width': 0.2}       # Larger effects
}


def promoter_to_tx_rate(
    seq: str,
    score_params: Optional[Dict] = None,
    add_noise: bool = True,
    noise_sigma: float = 0.05,
    use_log_normal: bool = True,
    severity_profile: str = 'standard',
    apply_neutral_zone: bool = True
) -> float:
    """
    Map a promoter sequence to a transcription rate multiplier.
    
    Uses heuristic scoring based on:
    1. Match to -35 box consensus (TTGACA)
    2. Match to -10 box consensus (TATAAT)
    3. Spacer length between boxes (optimal: 17 bp)
    
    Args:
        seq: Promoter DNA sequence (should be ~50-60 bp covering both boxes)
        score_params: Optional dict to override scoring parameters
            - 'minus35_weight': Weight for -35 box score (default: 0.35)
            - 'minus10_weight': Weight for -10 box score (default: 0.45)
            - 'spacer_weight': Weight for spacer score (default: 0.20)
        add_noise: Whether to add biological noise (default: True)
        noise_sigma: Standard deviation of Gaussian noise (default: 0.05)
    
    Returns:
        Transcription rate multiplier (0.1 to 1.5 range, clamped)
        - 1.0 = wild-type consensus activity
        - < 1.0 = reduced activity due to mutations
        - > 1.0 = enhanced activity (rare, via improved box sequences)
    
    Biological Notes:
        - The -10 box is more critical than -35 for most promoters
        - Spacer length affects helical phasing and polymerase binding
        - Natural promoters vary widely in strength (100-fold range)
    """
    seq = seq.upper()
    
    # Default scoring parameters
    params = {
        'minus35_weight': 0.35,
        'minus10_weight': 0.45,
        'spacer_weight': 0.20
    }
    if score_params:
        params.update(score_params)
    
    # Find best matches to -35 and -10 boxes
    minus35_match = _find_best_box_match(seq, MINUS_35_CONSENSUS, MINUS_35_PWM)
    minus10_match = _find_best_box_match(seq, MINUS_10_CONSENSUS, MINUS_10_PWM)
    
    # Calculate spacer score if both boxes found
    spacer_score = 0.0
    if minus35_match[1] is not None and minus10_match[1] is not None:
        spacer_length = minus10_match[1] - (minus35_match[1] + 6)  # -35 is 6 bp
        spacer_score = _score_spacer(spacer_length)
    
    # Combine scores
    combined_score = (
        params['minus35_weight'] * minus35_match[0] +
        params['minus10_weight'] * minus10_match[0] +
        params['spacer_weight'] * spacer_score
    )
    
    # Get severity profile
    severity = MUTATION_SEVERITY.get(severity_profile, MUTATION_SEVERITY['standard'])
    scale = severity['scale']
    
    # Convert score to multiplier
    # Score range: -1 to 1, map to multiplier range
    if use_log_normal:
        # Log-normal scaling: changes are multiplicative
        # This ensures parameters stay positive and changes scale properly
        log_multiplier = combined_score * scale * 0.5  # Map to ~[-0.5, 0.5] log space
        multiplier = np.exp(log_multiplier)
    else:
        # Original sigmoid-like transformation
        multiplier = 0.1 + 1.4 * (1 / (1 + np.exp(-3 * combined_score * scale)))
    
    # Apply neutral zone - mild mutations have minimal effect
    if apply_neutral_zone:
        neutral_low, neutral_high = NEUTRAL_ZONE
        neutral_width = severity['neutral_width']
        
        # If multiplier is close to 1.0, pull it back toward 1.0
        if neutral_low <= multiplier <= neutral_high:
            # Dampen the effect when in neutral zone
            deviation = multiplier - 1.0
            damped_deviation = deviation * (1 - neutral_width)
            multiplier = 1.0 + damped_deviation
    
    # Add biological noise (log-normal for multiplicative effects)
    if add_noise:
        if use_log_normal:
            noise_multiplier = np.exp(np.random.normal(0, noise_sigma))
            multiplier *= noise_multiplier
        else:
            noise = np.random.normal(0, noise_sigma)
            multiplier += noise
    
    # Clamp to biologically plausible range
    multiplier = np.clip(multiplier, 0.1, 1.5)
    
    return float(multiplier)


def _find_best_box_match(
    seq: str, 
    consensus: str, 
    pwm: Dict
) -> Tuple[float, Optional[int]]:
    """
    Find the best matching position for a consensus box in the sequence.
    
    Returns:
        Tuple of (normalized score [-1 to 1], position of best match)
    """
    if len(seq) < len(consensus):
        return -1.0, None
    
    best_score = float('-inf')
    best_pos = None
    
    for i in range(len(seq) - len(consensus) + 1):
        window = seq[i:i + len(consensus)]
        score = _score_box(window, pwm)
        if score > best_score:
            best_score = score
            best_pos = i
    
    # Normalize score to [-1, 1] range
    # Max possible score: sum of max values in PWM
    max_score = sum(max(pwm[i].values()) for i in range(len(consensus)))
    min_score = sum(min(pwm[i].values()) for i in range(len(consensus)))
    
    if max_score == min_score:
        normalized = 0.0
    else:
        normalized = 2 * (best_score - min_score) / (max_score - min_score) - 1
    
    return normalized, best_pos


def _score_box(window: str, pwm: Dict) -> float:
    """Score a sequence window against a PWM."""
    score = 0.0
    for i, base in enumerate(window):
        if base in pwm.get(i, {}):
            score += pwm[i][base]
        else:
            score += -1.0  # Unknown base penalty
    return score


def _score_spacer(length: int) -> float:
    """
    Score spacer length between -35 and -10 boxes.
    
    Optimal length is 17 bp with tolerance of ±2 bp.
    Score decreases linearly outside tolerance range.
    """
    deviation = abs(length - OPTIMAL_SPACER)
    
    if deviation <= SPACER_TOLERANCE:
        return 1.0
    else:
        # Linear decrease, bottoms out at 0 for deviation > 7
        return max(0.0, 1.0 - (deviation - SPACER_TOLERANCE) * 0.15)


def calculate_promoter_consensus_match(seq: str) -> Dict[str, float]:
    """
    Calculate detailed consensus matching scores for a promoter sequence.
    
    Useful for debugging and detailed analysis.
    
    Returns:
        Dict with keys: 'minus35_score', 'minus10_score', 'spacer_score', 
                       'minus35_pos', 'minus10_pos', 'total_score'
    """
    seq = seq.upper()
    
    minus35_match = _find_best_box_match(seq, MINUS_35_CONSENSUS, MINUS_35_PWM)
    minus10_match = _find_best_box_match(seq, MINUS_10_CONSENSUS, MINUS_10_PWM)
    
    spacer_score = 0.0
    spacer_length = None
    if minus35_match[1] is not None and minus10_match[1] is not None:
        spacer_length = minus10_match[1] - (minus35_match[1] + 6)
        spacer_score = _score_spacer(spacer_length)
    
    return {
        'minus35_score': minus35_match[0],
        'minus10_score': minus10_match[0],
        'spacer_score': spacer_score,
        'spacer_length': spacer_length,
        'minus35_pos': minus35_match[1],
        'minus10_pos': minus10_match[1],
        'total_score': 0.35 * minus35_match[0] + 0.45 * minus10_match[0] + 0.20 * spacer_score
    }
