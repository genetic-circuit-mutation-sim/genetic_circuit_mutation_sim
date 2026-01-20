"""
Sequence Mutator Module

Core mutation engine for introducing controlled mutations into DNA sequences.
Supports substitutions, insertions, and deletions with user-defined rates
and region-specific constraints.

Public API:
    mutate_sequence(seq, regions, n_mut, rates) -> (str, List[Mutation])
"""

import random
from typing import List, Tuple, Dict, Optional
import numpy as np

from .mutation_types import (
    Mutation, 
    MutationType, 
    AnnotatedRegion, 
    MutationRates,
    RegionType
)


# Standard DNA nucleotides
NUCLEOTIDES = ['A', 'T', 'G', 'C']

# Transition mutations (purine<->purine, pyrimidine<->pyrimidine) are more common
TRANSITIONS = {'A': 'G', 'G': 'A', 'T': 'C', 'C': 'T'}

# Transversion mutations (purine<->pyrimidine)
TRANSVERSIONS = {
    'A': ['T', 'C'],
    'G': ['T', 'C'],
    'T': ['A', 'G'],
    'C': ['A', 'G']
}


def mutate_sequence(
    seq: str,
    regions: List[AnnotatedRegion],
    n_mut: int = 1,
    rates: Optional[Dict[str, float]] = None,
    per_base_rate: Optional[float] = None,
    transition_bias: float = 2.0,
    seed: Optional[int] = None
) -> Tuple[str, List[Mutation]]:
    """
    Introduce controlled mutations into a DNA sequence.
    
    Mutations are restricted to annotated regions only. The function supports
    three mutation types: substitutions, insertions, and deletions, with
    user-controllable probabilities.
    
    Args:
        seq: Input DNA sequence (uppercase A, T, G, C)
        regions: List of AnnotatedRegion objects defining mutable regions
        n_mut: Number of mutations to introduce (default: 1)
        rates: Dict with keys 'substitution_rate', 'insertion_rate', 'deletion_rate'
               Values represent relative probabilities (will be normalized)
               Default: {'substitution_rate': 0.7, 'insertion_rate': 0.15, 'deletion_rate': 0.15}
        per_base_rate: Optional per-base mutation probability (overrides n_mut if set)
        transition_bias: Ratio of transitions to transversions (default: 2.0)
                        Biological rationale: transitions are ~2x more common in nature
        seed: Random seed for reproducibility
    
    Returns:
        Tuple of:
            - Mutated sequence string
            - List of Mutation objects describing each mutation
    
    Raises:
        ValueError: If sequence is empty, no regions provided, or invalid rates
    
    Example:
        >>> seq = "ATGCATGCATGC"
        >>> regions = [AnnotatedRegion("test", 0, 12, RegionType.CDS)]
        >>> mutated, mutations = mutate_sequence(seq, regions, n_mut=2)
    
    Biological Notes:
        - Substitutions are weighted toward transitions (e.g., A<->G, T<->C)
        - Insertions add 1-3 random nucleotides (mimicking replication errors)
        - Deletions remove 1-3 nucleotides (mimicking replication slippage)
        - Mutations in CDS can cause frameshifts if indel length % 3 != 0
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    # Input validation
    if not seq:
        raise ValueError("Sequence cannot be empty")
    if not regions:
        raise ValueError("At least one annotated region must be provided")
    
    seq = seq.upper()
    
    # Set default rates
    if rates is None:
        rates = {
            'substitution_rate': 0.7,
            'insertion_rate': 0.15,
            'deletion_rate': 0.15
        }
    
    mutation_rates = MutationRates(
        substitution_rate=rates.get('substitution_rate', 0.7),
        insertion_rate=rates.get('insertion_rate', 0.15),
        deletion_rate=rates.get('deletion_rate', 0.15)
    )
    
    # Build list of mutable positions (positions within annotated regions)
    mutable_positions = []
    position_to_region = {}
    
    for region in regions:
        for pos in range(region.start, min(region.end, len(seq))):
            if pos not in position_to_region:  # Avoid duplicates from overlapping regions
                mutable_positions.append(pos)
                position_to_region[pos] = region
    
    if not mutable_positions:
        raise ValueError("No mutable positions found within provided regions")
    
    # Determine number of mutations
    if per_base_rate is not None:
        # Calculate expected mutations based on per-base rate
        n_mut = np.random.binomial(len(mutable_positions), per_base_rate)
        n_mut = max(1, n_mut)  # Ensure at least 1 mutation
    
    # Cap mutations at number of available positions
    n_mut = min(n_mut, len(mutable_positions))
    
    # Select mutation positions
    mutation_positions = random.sample(mutable_positions, n_mut)
    mutation_positions.sort()  # Process in order for correct offset tracking
    
    # Generate mutations
    mutations = []
    mutated_seq = list(seq)
    offset = 0  # Track position offset due to indels
    
    for original_pos in mutation_positions:
        current_pos = original_pos + offset
        region = position_to_region[original_pos]
        
        # Select mutation type based on rates
        mutation_type = _select_mutation_type(mutation_rates)
        
        if mutation_type == MutationType.SUBSTITUTION:
            mutation, mutated_seq = _apply_substitution(
                mutated_seq, current_pos, original_pos, region, transition_bias
            )
        elif mutation_type == MutationType.INSERTION:
            mutation, mutated_seq, delta = _apply_insertion(
                mutated_seq, current_pos, original_pos, region
            )
            offset += delta
        else:  # DELETION
            mutation, mutated_seq, delta = _apply_deletion(
                mutated_seq, current_pos, original_pos, region
            )
            offset += delta
        
        mutations.append(mutation)
    
    return ''.join(mutated_seq), mutations


def _select_mutation_type(rates: MutationRates) -> MutationType:
    """Select mutation type based on normalized rates."""
    norm_rates = rates.normalized
    r = random.random()
    
    if r < norm_rates[0]:
        return MutationType.SUBSTITUTION
    elif r < norm_rates[0] + norm_rates[1]:
        return MutationType.INSERTION
    else:
        return MutationType.DELETION


def _apply_substitution(
    seq: List[str],
    current_pos: int,
    original_pos: int,
    region: AnnotatedRegion,
    transition_bias: float
) -> Tuple[Mutation, List[str]]:
    """
    Apply a substitution mutation at the specified position.
    
    Biological rationale: Transitions (purine<->purine, pyrimidine<->pyrimidine)
    are more common than transversions due to chemical similarity.
    """
    original_base = seq[current_pos]
    
    # Choose between transition and transversion
    if random.random() < transition_bias / (transition_bias + 1):
        # Transition
        new_base = TRANSITIONS.get(original_base, random.choice(NUCLEOTIDES))
    else:
        # Transversion
        new_base = random.choice(TRANSVERSIONS.get(original_base, NUCLEOTIDES))
    
    # Ensure we actually change the base
    while new_base == original_base:
        new_base = random.choice([n for n in NUCLEOTIDES if n != original_base])
    
    seq[current_pos] = new_base
    
    mutation = Mutation(
        position=original_pos,
        mutation_type=MutationType.SUBSTITUTION,
        original=original_base,
        replacement=new_base,
        region=region
    )
    
    return mutation, seq


def _apply_insertion(
    seq: List[str],
    current_pos: int,
    original_pos: int,
    region: AnnotatedRegion
) -> Tuple[Mutation, List[str], int]:
    """
    Apply an insertion mutation at the specified position.
    
    Biological rationale: Small insertions (1-3 bp) are common due to
    DNA polymerase slippage during replication, especially in repetitive regions.
    """
    # Insert 1-3 nucleotides (weighted toward smaller insertions)
    insert_length = random.choices([1, 2, 3], weights=[0.6, 0.3, 0.1])[0]
    inserted_bases = ''.join(random.choices(NUCLEOTIDES, k=insert_length))
    
    # Insert after current position
    seq = seq[:current_pos + 1] + list(inserted_bases) + seq[current_pos + 1:]
    
    mutation = Mutation(
        position=original_pos,
        mutation_type=MutationType.INSERTION,
        original='',
        replacement=inserted_bases,
        region=region
    )
    
    return mutation, seq, insert_length


def _apply_deletion(
    seq: List[str],
    current_pos: int,
    original_pos: int,
    region: AnnotatedRegion
) -> Tuple[Mutation, List[str], int]:
    """
    Apply a deletion mutation at the specified position.
    
    Biological rationale: Small deletions (1-3 bp) are common due to
    DNA polymerase slippage, similar to insertions.
    """
    # Delete 1-3 nucleotides (weighted toward smaller deletions)
    max_delete = min(3, len(seq) - current_pos)  # Don't delete past end
    if max_delete <= 0:
        max_delete = 1
    
    weights = [0.6, 0.3, 0.1][:max_delete]
    weights = [w / sum(weights) for w in weights]  # Renormalize
    
    delete_length = random.choices(range(1, max_delete + 1), weights=weights)[0]
    deleted_bases = ''.join(seq[current_pos:current_pos + delete_length])
    
    # Remove nucleotides
    seq = seq[:current_pos] + seq[current_pos + delete_length:]
    
    mutation = Mutation(
        position=original_pos,
        mutation_type=MutationType.DELETION,
        original=deleted_bases,
        replacement='',
        region=region
    )
    
    return mutation, seq, -delete_length


def generate_toggle_switch_sequence(gene_a_params: dict = None, gene_b_params: dict = None) -> Tuple[str, List[AnnotatedRegion]]:
    """
    Generate a synthetic toggle switch DNA sequence with annotated regions.
    
    This creates a simplified representation of the toggle switch genetic circuit
    with standard regulatory elements for both genes A and B.
    
    Structure for each gene:
        [Operator][Promoter][-35 box][spacer][-10 box][spacer][RBS][CDS]
    
    Returns:
        Tuple of:
            - Complete DNA sequence
            - List of annotated regions
    
    Biological Notes:
        - Promoter contains -35 and -10 boxes (sigma70 consensus)
        - RBS contains Shine-Dalgarno sequence
        - CDS is a simplified coding sequence
        - Operator contains repressor binding site
    """
    # Default consensus sequences
    operator_consensus = "AATTGTGAGCGGATAACAATT"  # lac operator-like
    promoter_minus35 = "TTGACA"  # -35 box consensus
    promoter_minus10 = "TATAAT"  # -10 box consensus
    rbs_sd = "AGGAGG"  # Shine-Dalgarno
    
    # Build gene A cassette
    gene_a_operator = operator_consensus
    gene_a_promoter = f"NNNNNN{promoter_minus35}NNNNNNNNNNNNNNNNN{promoter_minus10}NNNNNN"
    gene_a_rbs = f"NNNNNN{rbs_sd}NNNNNNNATG"  # Includes start codon context
    gene_a_cds = "ATGGCTAGCAAAGAATTCGGATCCCTGCAGGAATTCGATATCAAGCTTATCGATACCGTCGACTAA"  # Simplified CDS
    
    # Build gene B cassette (same structure)
    gene_b_operator = operator_consensus
    gene_b_promoter = f"NNNNNN{promoter_minus35}NNNNNNNNNNNNNNNNN{promoter_minus10}NNNNNN"
    gene_b_rbs = f"NNNNNN{rbs_sd}NNNNNNNATG"
    gene_b_cds = "ATGGATCCGAATTCAAGCTTGATATCGAATTCCTGCAGGGATCCGAATTCTTTGCTAGCCATTAA"
    
    # Replace N's with random nucleotides
    def fill_random(seq):
        return ''.join(random.choice(NUCLEOTIDES) if c == 'N' else c for c in seq)
    
    gene_a_operator = fill_random(gene_a_operator)
    gene_a_promoter = fill_random(gene_a_promoter)
    gene_a_rbs = fill_random(gene_a_rbs)
    gene_b_operator = fill_random(gene_b_operator)
    gene_b_promoter = fill_random(gene_b_promoter)
    gene_b_rbs = fill_random(gene_b_rbs)
    
    # Spacer between genes
    spacer = "NNNNNNNNNNNNNNNNNNNN"
    spacer = fill_random(spacer)
    
    # Assemble full sequence
    full_seq = (
        gene_a_operator + gene_a_promoter + gene_a_rbs + gene_a_cds +
        spacer +
        gene_b_operator + gene_b_promoter + gene_b_rbs + gene_b_cds
    )
    
    # Calculate region positions
    pos = 0
    regions = []
    
    # Gene A regions
    regions.append(AnnotatedRegion("operator_A", pos, pos + len(gene_a_operator), RegionType.OPERATOR, "A"))
    pos += len(gene_a_operator)
    
    regions.append(AnnotatedRegion("promoter_A", pos, pos + len(gene_a_promoter), RegionType.PROMOTER, "A"))
    pos += len(gene_a_promoter)
    
    regions.append(AnnotatedRegion("rbs_A", pos, pos + len(gene_a_rbs), RegionType.RBS, "A"))
    pos += len(gene_a_rbs)
    
    regions.append(AnnotatedRegion("cds_A", pos, pos + len(gene_a_cds), RegionType.CDS, "A"))
    pos += len(gene_a_cds)
    
    pos += len(spacer)  # Skip spacer (not annotated)
    
    # Gene B regions
    regions.append(AnnotatedRegion("operator_B", pos, pos + len(gene_b_operator), RegionType.OPERATOR, "B"))
    pos += len(gene_b_operator)
    
    regions.append(AnnotatedRegion("promoter_B", pos, pos + len(gene_b_promoter), RegionType.PROMOTER, "B"))
    pos += len(gene_b_promoter)
    
    regions.append(AnnotatedRegion("rbs_B", pos, pos + len(gene_b_rbs), RegionType.RBS, "B"))
    pos += len(gene_b_rbs)
    
    regions.append(AnnotatedRegion("cds_B", pos, pos + len(gene_b_cds), RegionType.CDS, "B"))
    
    return full_seq, regions
