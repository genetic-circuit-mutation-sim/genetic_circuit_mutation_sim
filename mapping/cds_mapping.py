"""
CDS Mapping Module

Maps coding sequence (CDS) mutations to protein function effects.
Handles detection of nonsense mutations, frameshifts, and operator binding effects.

Biological Background:
    CDS mutations can affect protein function through:
    - Nonsense mutations: Create premature stop codons, truncating the protein
    - Frameshift mutations: Indels not divisible by 3, scrambling downstream sequence
    - Missense mutations: Amino acid substitutions (effect varies)
    - Silent mutations: No amino acid change (usually neutral)
    
    Operator mutations affect repressor binding:
    - Perfect match = strong binding (low Kd)
    - Mismatches = weaker binding (higher Kd)

Public API:
    cds_to_protein_effect(seq) -> dict
    operator_to_binding_params(seq, original_seq) -> dict
"""

import numpy as np
from typing import Dict, List, Tuple, Optional


# Standard genetic code (DNA codons to amino acids)
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Stop codons
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

# Start codon
START_CODON = 'ATG'

# Operator consensus sequence (lac operator-like)
OPERATOR_CONSENSUS = "AATTGTGAGCGGATAACAATT"

# Neutral zone for buffered parameter mapping
NEUTRAL_ZONE = (0.7, 1.3)

# Mutation severity profiles
MUTATION_SEVERITY = {
    'conservative': {'scale': 0.5, 'neutral_width': 0.4},
    'standard': {'scale': 1.0, 'neutral_width': 0.3},
    'aggressive': {'scale': 1.5, 'neutral_width': 0.2}
}

# Codon table for synonymous mutation detection
# Group codons by amino acid
CODON_TO_AA = GENETIC_CODE.copy()
AA_TO_CODONS = {}
for codon, aa in CODON_TO_AA.items():
    if aa not in AA_TO_CODONS:
        AA_TO_CODONS[aa] = []
    AA_TO_CODONS[aa].append(codon)


def cds_to_protein_effect(
    seq: str,
    original_seq: Optional[str] = None,
    add_noise: bool = True,
    noise_sigma: float = 0.02,
    use_log_normal: bool = True,
    apply_neutral_zone: bool = True
) -> Dict:
    """
    Map CDS sequence mutations to protein function effects.
    
    Analyzes the coding sequence to detect:
    - Synonymous mutations (neutral - no amino acid change)
    - Premature stop codons (nonsense mutations)
    - Frameshifts from length changes
    - Missense mutations (amino acid substitutions)
    
    Args:
        seq: Mutated CDS DNA sequence (should start with ATG)
        original_seq: Original wild-type sequence (for frameshift detection)
        add_noise: Whether to add biological noise (default: True)
        noise_sigma: Standard deviation of noise (default: 0.02)
        use_log_normal: Use log-normal noise (default: True)
        apply_neutral_zone: Buffer mild mutations (default: True)
    
    Returns:
        Dict containing:
            - 'effect_multiplier': Protein activity multiplier (0.0 to 1.0)
            - 'effect_type': 'functional', 'synonymous', 'nonsense', 'frameshift', 'missense'
            - 'stop_position': Position of premature stop (if any)
            - 'protein_length': Predicted protein length
            - 'original_length': Original protein length (if original_seq provided)
            - 'synonymous_count': Number of synonymous changes
            - 'missense_count': Number of missense changes
    """
    seq = seq.upper()
    
    result = {
        'effect_multiplier': 1.0,
        'effect_type': 'functional',
        'stop_position': None,
        'protein_length': 0,
        'original_length': None,
        'amino_acid_changes': [],
        'synonymous_count': 0,
        'missense_count': 0
    }
    
    # Check for frameshift (length not multiple of 3)
    if len(seq) % 3 != 0:
        result['effect_type'] = 'frameshift'
        result['effect_multiplier'] = 0.0
        return result
    
    # Compare to original if provided
    if original_seq:
        original_seq = original_seq.upper()
        len_diff = len(seq) - len(original_seq)
        
        if len_diff % 3 != 0:
            result['effect_type'] = 'frameshift'
            result['effect_multiplier'] = 0.0
            return result
        
        result['original_length'] = _count_codons_to_stop(original_seq)
    
    # Translate and find premature stops
    protein_length, stop_pos = _find_first_stop(seq)
    result['protein_length'] = protein_length
    
    # Check for premature stop (nonsense mutation)
    expected_length = len(seq) // 3 - 1
    
    if stop_pos is not None and stop_pos < expected_length - 1:
        result['effect_type'] = 'nonsense'
        result['stop_position'] = stop_pos
        truncation_fraction = stop_pos / max(expected_length, 1)
        result['effect_multiplier'] = 0.0 if truncation_fraction < 0.8 else 0.1
        return result
    
    # Analyze codon changes: separate synonymous from missense
    if original_seq and len(original_seq) == len(seq):
        synonymous, missense = _classify_codon_changes(seq, original_seq)
        result['synonymous_count'] = len(synonymous)
        result['missense_count'] = len(missense)
        result['amino_acid_changes'] = missense
        
        # Synonymous mutations are neutral
        if len(missense) == 0 and len(synonymous) > 0:
            result['effect_type'] = 'synonymous'
            result['effect_multiplier'] = 1.0  # Explicitly neutral
            return result
        
        if missense:
            # Estimate effect based on missense count only
            # Apply neutral zone buffering for mild effects
            base_penalty = len(missense) * 0.15
            multiplier = max(0.2, 1.0 - base_penalty)
            
            if apply_neutral_zone:
                # Mild mutations (1-2 missense) get buffered
                if len(missense) <= 2 and multiplier >= NEUTRAL_ZONE[0]:
                    # Dampen the effect
                    deviation = multiplier - 1.0
                    damped = deviation * 0.7  # 30% dampening
                    multiplier = 1.0 + damped
            
            result['effect_multiplier'] = multiplier
            result['effect_type'] = 'missense'
    
    # Add biological noise
    if add_noise and result['effect_multiplier'] > 0:
        if use_log_normal:
            noise_mult = np.exp(np.random.normal(0, noise_sigma))
            result['effect_multiplier'] = np.clip(
                result['effect_multiplier'] * noise_mult,
                0.0, 1.0
            )
        else:
            noise = np.random.normal(0, noise_sigma)
            result['effect_multiplier'] = np.clip(
                result['effect_multiplier'] + noise,
                0.0, 1.0
            )
    
    return result


def _classify_codon_changes(
    mutant_seq: str, 
    original_seq: str
) -> Tuple[List[Dict], List[Dict]]:
    """
    Classify codon changes as synonymous or missense.
    
    Returns:
        Tuple of (synonymous_changes, missense_changes)
    """
    synonymous = []
    missense = []
    
    min_len = min(len(mutant_seq), len(original_seq))
    
    for i in range(0, min_len - 2, 3):
        orig_codon = original_seq[i:i+3]
        mut_codon = mutant_seq[i:i+3]
        
        if orig_codon != mut_codon:
            orig_aa = GENETIC_CODE.get(orig_codon, 'X')
            mut_aa = GENETIC_CODE.get(mut_codon, 'X')
            
            change = {
                'position': i // 3,
                'original_codon': orig_codon,
                'mutant_codon': mut_codon,
                'original_aa': orig_aa,
                'mutant_aa': mut_aa
            }
            
            if orig_aa == mut_aa:
                # Synonymous (silent) mutation
                synonymous.append(change)
            else:
                # Missense mutation
                missense.append(change)
    
    return synonymous, missense


def _find_first_stop(seq: str) -> Tuple[int, Optional[int]]:
    """
    Find the first stop codon in a sequence.
    
    Returns:
        Tuple of (number of codons before stop, position of stop codon)
    """
    codon_count = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in STOP_CODONS:
            return codon_count, codon_count
        codon_count += 1
    
    return codon_count, None


def _count_codons_to_stop(seq: str) -> int:
    """Count codons until first stop codon."""
    count, _ = _find_first_stop(seq)
    return count


def _find_amino_acid_changes(mutant_seq: str, original_seq: str) -> List[Dict]:
    """
    Compare two sequences and find amino acid changes.
    
    Returns:
        List of dicts with 'position', 'original_aa', 'mutant_aa'
    """
    changes = []
    
    min_len = min(len(mutant_seq), len(original_seq))
    
    for i in range(0, min_len - 2, 3):
        orig_codon = original_seq[i:i+3]
        mut_codon = mutant_seq[i:i+3]
        
        if orig_codon != mut_codon:
            orig_aa = GENETIC_CODE.get(orig_codon, 'X')
            mut_aa = GENETIC_CODE.get(mut_codon, 'X')
            
            if orig_aa != mut_aa:
                changes.append({
                    'position': i // 3,
                    'original_aa': orig_aa,
                    'mutant_aa': mut_aa,
                    'original_codon': orig_codon,
                    'mutant_codon': mut_codon
                })
    
    return changes


def operator_to_binding_params(
    seq: str,
    original_seq: Optional[str] = None,
    add_noise: bool = True,
    noise_sigma: float = 0.05
) -> Dict:
    """
    Map operator sequence to repressor binding parameters.
    
    Operator mutations affect DNA-binding protein (repressor) affinity,
    which translates to changes in Kd (dissociation constant) and
    potentially Hill coefficient n.
    
    Args:
        seq: Mutated operator DNA sequence
        original_seq: Original wild-type operator sequence
        add_noise: Whether to add biological noise
        noise_sigma: Standard deviation of noise
    
    Returns:
        Dict containing:
            - 'kd_multiplier': Multiplier for dissociation constant Kd
                              > 1.0 = weaker binding (mutations disrupt)
                              < 1.0 = stronger binding (rare, improved)
            - 'n_multiplier': Multiplier for Hill coefficient
            - 'mismatch_count': Number of mismatches from consensus
            - 'binding_score': Overall binding strength score
    
    Biological Notes:
        - Operator sequences are palindromic for dimeric repressors
        - Mismatches in critical positions severely reduce binding
        - Complete loss of binding = constitutive expression
    """
    seq = seq.upper()
    
    if original_seq is None:
        original_seq = OPERATOR_CONSENSUS
    else:
        original_seq = original_seq.upper()
    
    # Calculate mismatch count
    mismatches = _count_mismatches(seq, original_seq)
    
    # Score sequence similarity to consensus
    consensus_score = _score_operator_match(seq, OPERATOR_CONSENSUS)
    
    # Calculate Kd multiplier
    # Each mismatch increases Kd (weaker binding)
    # Strong effect: 2-fold increase per mismatch
    kd_multiplier = 2.0 ** mismatches
    
    # Clamp to biologically plausible range
    kd_multiplier = np.clip(kd_multiplier, 0.5, 1000.0)
    
    # Hill coefficient generally decreases with weaker binding
    # due to reduced cooperativity
    if mismatches > 2:
        n_multiplier = max(0.5, 1.0 - mismatches * 0.1)
    else:
        n_multiplier = 1.0
    
    # Add noise
    if add_noise:
        kd_noise = np.random.normal(0, noise_sigma * kd_multiplier)
        kd_multiplier = np.clip(kd_multiplier + kd_noise, 0.5, 1000.0)
        
        n_noise = np.random.normal(0, noise_sigma * 0.5)
        n_multiplier = np.clip(n_multiplier + n_noise, 0.5, 1.2)
    
    return {
        'kd_multiplier': float(kd_multiplier),
        'n_multiplier': float(n_multiplier),
        'mismatch_count': mismatches,
        'binding_score': consensus_score
    }


def _count_mismatches(seq1: str, seq2: str) -> int:
    """Count number of mismatched positions between two sequences."""
    min_len = min(len(seq1), len(seq2))
    mismatches = sum(1 for i in range(min_len) if seq1[i] != seq2[i])
    
    # Length differences also count as mismatches
    mismatches += abs(len(seq1) - len(seq2))
    
    return mismatches


def _score_operator_match(seq: str, consensus: str) -> float:
    """
    Score operator sequence match to consensus.
    
    Returns score from 0.0 (no match) to 1.0 (perfect match).
    """
    if not seq or not consensus:
        return 0.0
    
    min_len = min(len(seq), len(consensus))
    if min_len == 0:
        return 0.0
    
    matches = sum(1 for i in range(min_len) if seq[i] == consensus[i])
    
    # Penalize length differences
    len_penalty = abs(len(seq) - len(consensus)) / max(len(seq), len(consensus))
    
    score = (matches / min_len) * (1 - len_penalty)
    
    return max(0.0, score)


def translate_sequence(dna_seq: str) -> str:
    """
    Translate a DNA sequence to amino acid sequence.
    
    Args:
        dna_seq: DNA sequence (will be upper-cased)
    
    Returns:
        Amino acid sequence (single letter code)
    """
    dna_seq = dna_seq.upper()
    protein = []
    
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        aa = GENETIC_CODE.get(codon, 'X')
        
        if aa == '*':  # Stop codon
            break
        protein.append(aa)
    
    return ''.join(protein)
