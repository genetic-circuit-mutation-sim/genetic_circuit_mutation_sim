"""
Mutation Engine Module

Provides tools for introducing controlled mutations into DNA sequences,
with support for different mutation types and region-specific constraints.

Public API:
    - mutate_sequence: Main function to generate mutated sequences
    - MutationType: Enum of mutation types
    - Mutation: Dataclass representing a single mutation
    - AnnotatedRegion: Dataclass for sequence region annotations
"""

from .mutation_types import MutationType, Mutation, AnnotatedRegion
from .sequence_mutator import mutate_sequence

__all__ = ['mutate_sequence', 'MutationType', 'Mutation', 'AnnotatedRegion']
