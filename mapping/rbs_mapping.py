"""
RBS Mapping Module

Maps ribosome binding site (RBS) sequences to translation rate multipliers
using heuristic scoring based on Shine-Dalgarno sequence complementarity.

Biological Background:
    The ribosome binding site in bacteria contains:
    - Shine-Dalgarno (SD) sequence: complementary to 16S rRNA 3' end
    - Consensus SD: AGGAGG (anti-SD in 16S rRNA: CCUCCU)
    - Optimal spacing: 5-9 nucleotides between SD and start codon (AUG)
    
    Translation initiation rate depends on:
    1. Complementarity to anti-SD sequence
    2. Spacing to start codon
    3. Local mRNA secondary structure (simplified in this model)

Public API:
    rbs_to_tl_rate(seq) -> float
"""

import numpy as np
from typing import Dict, Tuple, Optional


# Shine-Dalgarno consensus sequence
SD_CONSENSUS = "AGGAGG"

# Anti-Shine-Dalgarno (16S rRNA complement) for scoring
ANTI_SD = "TCCCTC"  # DNA equivalent of UCCUCC

# Optimal spacing between SD and start codon
OPTIMAL_SPACING_MIN = 5
OPTIMAL_SPACING_MAX = 9
SPACING_TOLERANCE = 2

# Neutral zone for buffered parameter mapping
NEUTRAL_ZONE = (0.7, 1.3)

# Mutation severity profiles
MUTATION_SEVERITY = {
    'conservative': {'scale': 0.5, 'neutral_width': 0.4},
    'standard': {'scale': 1.0, 'neutral_width': 0.3},
    'aggressive': {'scale': 1.5, 'neutral_width': 0.2}
}


# Position weight matrix for SD sequence
# Based on natural RBS variation studies
SD_PWM = {
    0: {'A': 1.0, 'G': 0.3, 'T': -0.5, 'C': -0.5},
    1: {'G': 1.0, 'A': 0.3, 'T': -0.5, 'C': -0.5},
    2: {'G': 1.0, 'A': 0.3, 'T': -0.5, 'C': -0.5},
    3: {'A': 1.0, 'G': 0.5, 'T': -0.3, 'C': -0.5},
    4: {'G': 1.0, 'A': 0.3, 'T': -0.5, 'C': -0.3},
    5: {'G': 1.0, 'A': 0.3, 'T': -0.3, 'C': -0.5},
}


def rbs_to_tl_rate(
    seq: str,
    add_noise: bool = True,
    noise_sigma: float = 0.05,
    use_log_normal: bool = True,
    severity_profile: str = 'standard',
    apply_neutral_zone: bool = True
) -> float:
    """
    Map an RBS sequence to a translation rate multiplier.
    
    Uses heuristic scoring based on:
    1. Shine-Dalgarno sequence match (AGGAGG consensus)
    2. Spacing between SD and start codon (AUG)
    3. AU content around start codon (accessibility)
    
    Args:
        seq: RBS DNA sequence (should include SD region and start codon context)
             Typically ~20-30 bp covering the region upstream of start codon
        add_noise: Whether to add biological noise (default: True)
        noise_sigma: Standard deviation of Gaussian noise (default: 0.05)
        use_log_normal: Use log-normal scaling (default: True)
        severity_profile: Mutation severity level (default: 'standard')
        apply_neutral_zone: Buffer mild mutations (default: True)
    
    Returns:
        Translation rate multiplier (0.1 to 1.5 range, clamped)
        - 1.0 = optimal RBS activity
        - < 1.0 = reduced translation initiation
        - > 1.0 = enhanced translation (rare)
    """
    seq = seq.upper()
    
    # Find Shine-Dalgarno match and position
    sd_score, sd_pos = _find_sd_match(seq)
    
    # Find start codon (AUG/ATG)
    start_pos = _find_start_codon(seq)
    
    # Calculate spacing score
    spacing_score = 0.0
    if sd_pos is not None and start_pos is not None:
        spacing = start_pos - (sd_pos + len(SD_CONSENSUS))
        spacing_score = _score_spacing(spacing)
    
    # Calculate AU content score (accessibility proxy)
    au_score = _calculate_au_content(seq)
    
    # Combine scores
    combined_score = (
        0.55 * sd_score +
        0.30 * spacing_score +
        0.15 * au_score
    )
    
    # Get severity profile
    severity = MUTATION_SEVERITY.get(severity_profile, MUTATION_SEVERITY['standard'])
    scale = severity['scale']
    
    # Convert to multiplier
    if use_log_normal:
        log_multiplier = combined_score * scale * 0.5
        multiplier = np.exp(log_multiplier)
    else:
        multiplier = 0.1 + 1.4 * (1 / (1 + np.exp(-3 * combined_score * scale)))
    
    # Apply neutral zone
    if apply_neutral_zone:
        neutral_low, neutral_high = NEUTRAL_ZONE
        neutral_width = severity['neutral_width']
        
        if neutral_low <= multiplier <= neutral_high:
            deviation = multiplier - 1.0
            damped_deviation = deviation * (1 - neutral_width)
            multiplier = 1.0 + damped_deviation
    
    # Add biological noise
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


def _find_sd_match(seq: str) -> Tuple[float, Optional[int]]:
    """
    Find best Shine-Dalgarno match in the sequence.
    
    Returns:
        Tuple of (normalized score [-1 to 1], position of best match)
    """
    if len(seq) < len(SD_CONSENSUS):
        return -1.0, None
    
    best_score = float('-inf')
    best_pos = None
    
    for i in range(len(seq) - len(SD_CONSENSUS) + 1):
        window = seq[i:i + len(SD_CONSENSUS)]
        score = _score_sd(window)
        if score > best_score:
            best_score = score
            best_pos = i
    
    # Normalize score
    max_score = sum(max(SD_PWM[i].values()) for i in range(len(SD_CONSENSUS)))
    min_score = sum(min(SD_PWM[i].values()) for i in range(len(SD_CONSENSUS)))
    
    if max_score == min_score:
        normalized = 0.0
    else:
        normalized = 2 * (best_score - min_score) / (max_score - min_score) - 1
    
    return normalized, best_pos


def _score_sd(window: str) -> float:
    """Score a sequence window against the SD PWM."""
    score = 0.0
    for i, base in enumerate(window):
        if base in SD_PWM.get(i, {}):
            score += SD_PWM[i][base]
        else:
            score += -1.0
    return score


def _find_start_codon(seq: str) -> Optional[int]:
    """Find the position of the start codon (ATG) in the sequence."""
    # Look for ATG
    pos = seq.find('ATG')
    if pos != -1:
        return pos
    
    # Also check for alternative start codons (GTG, TTG)
    for alt_start in ['GTG', 'TTG']:
        pos = seq.find(alt_start)
        if pos != -1:
            return pos
    
    return None


def _score_spacing(spacing: int) -> float:
    """
    Score the spacing between SD and start codon.
    
    Optimal spacing is 5-9 nucleotides.
    """
    if OPTIMAL_SPACING_MIN <= spacing <= OPTIMAL_SPACING_MAX:
        return 1.0
    elif spacing < OPTIMAL_SPACING_MIN:
        # Too close - ribosome can't bind properly
        deviation = OPTIMAL_SPACING_MIN - spacing
        return max(-1.0, 1.0 - deviation * 0.3)
    else:
        # Too far - reduced pairing efficiency
        deviation = spacing - OPTIMAL_SPACING_MAX
        return max(-1.0, 1.0 - deviation * 0.15)


def _calculate_au_content(seq: str) -> float:
    """
    Calculate AU content as a proxy for mRNA accessibility.
    
    Higher AU content generally means less secondary structure,
    making the RBS more accessible to ribosomes.
    
    Returns score from -1 to 1.
    """
    if not seq:
        return 0.0
    
    au_count = seq.count('A') + seq.count('T')  # DNA equivalent of A and U
    au_fraction = au_count / len(seq)
    
    # Optimal AU content is around 50-60%
    # Score linearly from this optimum
    optimal = 0.55
    deviation = abs(au_fraction - optimal)
    
    return 1.0 - deviation * 2  # Score from -1 to 1


def calculate_rbs_detailed_score(seq: str) -> Dict[str, float]:
    """
    Calculate detailed RBS scoring information.
    
    Useful for debugging and analysis.
    
    Returns:
        Dict with detailed scoring components
    """
    seq = seq.upper()
    
    sd_score, sd_pos = _find_sd_match(seq)
    start_pos = _find_start_codon(seq)
    
    spacing = None
    spacing_score = 0.0
    if sd_pos is not None and start_pos is not None:
        spacing = start_pos - (sd_pos + len(SD_CONSENSUS))
        spacing_score = _score_spacing(spacing)
    
    au_score = _calculate_au_content(seq)
    
    return {
        'sd_score': sd_score,
        'sd_position': sd_pos,
        'start_codon_position': start_pos,
        'spacing': spacing,
        'spacing_score': spacing_score,
        'au_content_score': au_score,
        'total_score': 0.55 * sd_score + 0.30 * spacing_score + 0.15 * au_score
    }
