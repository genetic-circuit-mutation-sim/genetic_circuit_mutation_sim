"""
Mapping Module

Provides functions to convert DNA sequence mutations into parameter perturbations
for the toggle switch model. Each mapping function uses heuristic scoring based
on consensus sequence motifs.

Public API:
    - promoter_to_tx_rate: Map promoter sequence to transcription rate multiplier
    - rbs_to_tl_rate: Map RBS sequence to translation rate multiplier
    - cds_to_protein_effect: Map CDS mutations to protein function effects
"""

from .promoter_mapping import promoter_to_tx_rate
from .rbs_mapping import rbs_to_tl_rate
from .cds_mapping import cds_to_protein_effect, operator_to_binding_params

__all__ = [
    'promoter_to_tx_rate',
    'rbs_to_tl_rate', 
    'cds_to_protein_effect',
    'operator_to_binding_params'
]
