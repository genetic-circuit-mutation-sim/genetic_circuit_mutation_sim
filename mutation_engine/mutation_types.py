"""
Mutation Types Module

Defines data structures for representing mutations and annotated genomic regions.
These types are used throughout the mutation engine to track and classify
sequence modifications.

Biological context:
- Substitutions: Single nucleotide changes (transitions/transversions)
- Insertions: Addition of nucleotides (can cause frameshifts in CDS)
- Deletions: Removal of nucleotides (can cause frameshifts in CDS)
"""

from enum import Enum
from dataclasses import dataclass, field
from typing import Optional


class MutationType(Enum):
    """
    Enumeration of possible mutation types.
    
    SUBSTITUTION: Single nucleotide replacement (most common mutation type)
    INSERTION: Addition of one or more nucleotides
    DELETION: Removal of one or more nucleotides
    """
    SUBSTITUTION = "substitution"
    INSERTION = "insertion"
    DELETION = "deletion"


class RegionType(Enum):
    """
    Enumeration of functional genomic regions.
    
    Each region type has different effects when mutated:
    - PROMOTER: Affects transcription rate (tx)
    - RBS: Affects translation rate (tl)
    - CDS: Affects protein function (nonsense, missense, frameshift)
    - OPERATOR: Affects repressor binding (Kd, Hill coefficient n)
    """
    PROMOTER = "promoter"
    RBS = "rbs"
    CDS = "cds"
    OPERATOR = "operator"


@dataclass
class AnnotatedRegion:
    """
    Represents a functional region within a DNA sequence.
    
    Attributes:
        name: Descriptive name (e.g., "promoter_A", "RBS_B")
        start: Start position in sequence (0-indexed, inclusive)
        end: End position in sequence (0-indexed, exclusive)
        region_type: Type of functional region
        gene: Associated gene name (e.g., "A" or "B" for toggle switch)
    
    Example:
        >>> region = AnnotatedRegion("promoter_A", 0, 50, RegionType.PROMOTER, "A")
    """
    name: str
    start: int
    end: int
    region_type: RegionType
    gene: Optional[str] = None
    
    def __post_init__(self):
        """Validate region boundaries."""
        if self.start < 0:
            raise ValueError(f"Region start must be >= 0, got {self.start}")
        if self.end <= self.start:
            raise ValueError(f"Region end must be > start, got end={self.end}, start={self.start}")
    
    def contains(self, position: int) -> bool:
        """Check if a position falls within this region."""
        return self.start <= position < self.end
    
    @property
    def length(self) -> int:
        """Return the length of the region."""
        return self.end - self.start


@dataclass
class Mutation:
    """
    Represents a single mutation event.
    
    Attributes:
        position: Position in the original sequence where mutation occurs
        mutation_type: Type of mutation (substitution, insertion, deletion)
        original: Original nucleotide(s) at this position
        replacement: New nucleotide(s) (empty string for deletions)
        region: The annotated region where this mutation occurred
        
    Biological notes:
    - For substitutions: len(original) == 1 and len(replacement) == 1
    - For insertions: original can be empty, replacement contains inserted bases
    - For deletions: replacement is empty, original contains deleted bases
    """
    position: int
    mutation_type: MutationType
    original: str
    replacement: str
    region: Optional[AnnotatedRegion] = None
    
    def __str__(self) -> str:
        """Human-readable mutation notation."""
        region_str = f" in {self.region.name}" if self.region else ""
        
        if self.mutation_type == MutationType.SUBSTITUTION:
            return f"{self.original}{self.position+1}{self.replacement}{region_str}"
        elif self.mutation_type == MutationType.INSERTION:
            return f"ins{self.position+1}{self.replacement}{region_str}"
        else:  # DELETION
            return f"del{self.position+1}{self.original}{region_str}"
    
    @property
    def is_frameshift(self) -> bool:
        """
        Check if mutation causes a frameshift (for CDS mutations).
        
        Frameshifts occur when indels are not multiples of 3,
        disrupting the reading frame.
        """
        if self.mutation_type == MutationType.SUBSTITUTION:
            return False
        
        len_diff = len(self.replacement) - len(self.original)
        return len_diff % 3 != 0


@dataclass 
class MutationRates:
    """
    Container for mutation rate parameters.
    
    Attributes:
        substitution_rate: Probability of substitution (0-1)
        insertion_rate: Probability of insertion (0-1)
        deletion_rate: Probability of deletion (0-1)
        
    Note: Rates are normalized internally so they sum to 1.
    """
    substitution_rate: float = 0.7
    insertion_rate: float = 0.15
    deletion_rate: float = 0.15
    
    def __post_init__(self):
        """Validate and normalize rates."""
        if any(r < 0 for r in [self.substitution_rate, self.insertion_rate, self.deletion_rate]):
            raise ValueError("All rates must be non-negative")
        
        total = self.substitution_rate + self.insertion_rate + self.deletion_rate
        if total == 0:
            raise ValueError("At least one rate must be > 0")
    
    @property
    def normalized(self) -> tuple:
        """Return normalized rates that sum to 1."""
        total = self.substitution_rate + self.insertion_rate + self.deletion_rate
        return (
            self.substitution_rate / total,
            self.insertion_rate / total,
            self.deletion_rate / total
        )
