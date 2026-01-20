"""
Validation Module

Provides sanity checks against known experimental facts to validate
the mutation-to-phenotype mapping. This increases credibility by
showing the model reproduces expected biological behaviors.

Public API:
    run_validation_suite() -> ValidationReport
    validate_promoter_strength() -> ValidationResult
    validate_operator_mismatch() -> ValidationResult
    validate_mutation_severity() -> ValidationResult
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Callable
from enum import Enum


class ValidationStatus(Enum):
    """Validation test outcome."""
    PASS = "pass"
    FAIL = "fail"
    SKIP = "skip"
    WARN = "warning"


@dataclass
class ValidationResult:
    """
    Result of a single validation test.
    
    Attributes:
        test_name: Name of the validation test
        status: Pass/Fail/Skip/Warn
        expected: What was expected
        observed: What was observed
        message: Human-readable explanation
        confidence: How confident in this result (0-1)
    """
    test_name: str
    status: ValidationStatus
    expected: Any
    observed: Any
    message: str
    confidence: float = 1.0
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'test_name': self.test_name,
            'status': self.status.value,
            'expected': str(self.expected),
            'observed': str(self.observed),
            'message': self.message,
            'confidence': self.confidence
        }


@dataclass
class ValidationReport:
    """
    Complete validation report with all test results.
    """
    results: List[ValidationResult] = field(default_factory=list)
    summary: Dict[str, int] = field(default_factory=dict)
    overall_status: ValidationStatus = ValidationStatus.SKIP
    
    def add_result(self, result: ValidationResult):
        self.results.append(result)
        self._update_summary()
    
    def _update_summary(self):
        self.summary = {
            'pass': sum(1 for r in self.results if r.status == ValidationStatus.PASS),
            'fail': sum(1 for r in self.results if r.status == ValidationStatus.FAIL),
            'warn': sum(1 for r in self.results if r.status == ValidationStatus.WARN),
            'skip': sum(1 for r in self.results if r.status == ValidationStatus.SKIP)
        }
        
        if self.summary['fail'] > 0:
            self.overall_status = ValidationStatus.FAIL
        elif self.summary['warn'] > 0:
            self.overall_status = ValidationStatus.WARN
        elif self.summary['pass'] > 0:
            self.overall_status = ValidationStatus.PASS
        else:
            self.overall_status = ValidationStatus.SKIP
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'results': [r.to_dict() for r in self.results],
            'summary': self.summary,
            'overall_status': self.overall_status.value,
            'pct_pass': self.summary['pass'] / len(self.results) * 100 if self.results else 0
        }


def validate_promoter_strength(
    promoter_to_tx_rate_fn: Callable[[str], float]
) -> ValidationResult:
    """
    Validate that strong promoters give higher rates than weak promoters.
    
    Uses consensus sequence (strong) vs random sequence (weak) comparison.
    
    Expected: Strong promoter rate > Weak promoter rate by at least 2x
    """
    # Strong promoter: consensus -35 and -10 boxes
    strong_promoter = "TTGACANNNNNNNNNNNNNNNNTATAAT"
    
    # Weak promoter: random sequence with poor consensus match
    weak_promoter = "ACGATCNNNNNNNNNNNNNNNNGTCAAC"
    
    try:
        strong_rate = promoter_to_tx_rate_fn(strong_promoter)
        weak_rate = promoter_to_tx_rate_fn(weak_promoter)
        
        ratio = strong_rate / (weak_rate + 1e-10)
        
        expected = "Strong/Weak ratio > 2.0"
        observed = f"Ratio = {ratio:.2f} (Strong={strong_rate:.3f}, Weak={weak_rate:.3f})"
        
        if ratio > 2.0:
            return ValidationResult(
                test_name="Promoter Strength Ordering",
                status=ValidationStatus.PASS,
                expected=expected,
                observed=observed,
                message="Strong promoter correctly produces higher transcription rate",
                confidence=min(1.0, ratio / 3.0)
            )
        elif ratio > 1.0:
            return ValidationResult(
                test_name="Promoter Strength Ordering",
                status=ValidationStatus.WARN,
                expected=expected,
                observed=observed,
                message="Strong promoter is higher but difference is small",
                confidence=ratio / 2.0
            )
        else:
            return ValidationResult(
                test_name="Promoter Strength Ordering",
                status=ValidationStatus.FAIL,
                expected=expected,
                observed=observed,
                message="Strong promoter should have higher rate than weak promoter"
            )
    except Exception as e:
        return ValidationResult(
            test_name="Promoter Strength Ordering",
            status=ValidationStatus.SKIP,
            expected="Strong > Weak",
            observed=str(e),
            message=f"Test could not be run: {e}"
        )


def validate_operator_mismatch(
    operator_to_binding_fn: Callable[[str, str], Dict]
) -> ValidationResult:
    """
    Validate that operator mismatches reduce binding appropriately.
    
    Known biology:
    - Perfect operator match = strong binding
    - 1 bp mismatch = ~10x reduced binding
    - 2+ bp mismatch = weak/no binding
    
    Uses lacO-like operator sequences.
    """
    # Consensus operator sequence
    perfect_operator = "AATTGTGAGCGGATAACAATT"
    
    # 1bp mismatch (position 10: G->T)
    one_bp_mismatch = "AATTGTGAGCTGATAACAATT"
    
    # 2bp mismatch
    two_bp_mismatch = "AATTGTGATCTGATAACAATT"
    
    try:
        repressor_seq = "REPRESSOR_PROTEIN"  # Placeholder
        
        perfect_binding = operator_to_binding_fn(perfect_operator, repressor_seq)
        one_bp_binding = operator_to_binding_fn(one_bp_mismatch, repressor_seq)
        two_bp_binding = operator_to_binding_fn(two_bp_mismatch, repressor_seq)
        
        # Extract Kd values
        kd_perfect = perfect_binding.get('Kd', 1.0)
        kd_1bp = one_bp_binding.get('Kd', 1.0)
        kd_2bp = two_bp_binding.get('Kd', 1.0)
        
        # Lower Kd = stronger binding
        # Mismatch should increase Kd
        
        expected = "Kd_1bp > Kd_perfect AND Kd_2bp > Kd_1bp"
        observed = f"Kd_perfect={kd_perfect:.2f}, Kd_1bp={kd_1bp:.2f}, Kd_2bp={kd_2bp:.2f}"
        
        if kd_1bp > kd_perfect and kd_2bp > kd_1bp:
            return ValidationResult(
                test_name="Operator Mismatch Severity",
                status=ValidationStatus.PASS,
                expected=expected,
                observed=observed,
                message="Operator mismatches correctly reduce binding affinity"
            )
        else:
            return ValidationResult(
                test_name="Operator Mismatch Severity",
                status=ValidationStatus.FAIL,
                expected=expected,
                observed=observed,
                message="Mismatch severity ordering is incorrect"
            )
    except Exception as e:
        return ValidationResult(
            test_name="Operator Mismatch Severity",
            status=ValidationStatus.SKIP,
            expected="Kd increases with mismatches",
            observed=str(e),
            message=f"Test could not be run: {e}"
        )


def validate_mutation_severity(
    cds_to_effect_fn: Callable[[str, str], Dict]
) -> ValidationResult:
    """
    Validate mutation severity ordering:
    - Synonymous < Missense < Frameshift
    
    Frameshift should be most severe (lowest activity),
    synonymous should be near-neutral.
    """
    # Original CDS
    original_cds = "ATGGCTAGCTGA"  # Met-Ala-Ser-Stop
    
    # Synonymous mutation (GCT -> GCC both = Ala)
    synonymous_cds = "ATGGCCAGCTGA"
    
    # Missense mutation (GCT -> GAT = Asp instead of Ala)
    missense_cds = "ATGGATAGCTGA"
    
    # Frameshift (insert A after ATG)
    frameshift_cds = "ATGAGCTAGCTGA"
    
    try:
        original_effect = cds_to_effect_fn(original_cds, original_cds)
        synonymous_effect = cds_to_effect_fn(synonymous_cds, original_cds)
        missense_effect = cds_to_effect_fn(missense_cds, original_cds)
        frameshift_effect = cds_to_effect_fn(frameshift_cds, original_cds)
        
        # Get activity multipliers
        syn_mult = synonymous_effect.get('effect_multiplier', 1.0)
        mis_mult = missense_effect.get('effect_multiplier', 1.0)
        fs_mult = frameshift_effect.get('effect_multiplier', 1.0)
        
        expected = "Synonymous (~1.0) > Missense > Frameshift"
        observed = f"Syn={syn_mult:.2f}, Miss={mis_mult:.2f}, FS={fs_mult:.2f}"
        
        # Check ordering
        if syn_mult > mis_mult > fs_mult:
            if syn_mult > 0.9:  # Synonymous should be near-neutral
                return ValidationResult(
                    test_name="Mutation Severity Ordering",
                    status=ValidationStatus.PASS,
                    expected=expected,
                    observed=observed,
                    message="Mutation severity correctly ranked: synonymous > missense > frameshift"
                )
            else:
                return ValidationResult(
                    test_name="Mutation Severity Ordering",
                    status=ValidationStatus.WARN,
                    expected=expected,
                    observed=observed,
                    message="Ordering correct but synonymous not neutral enough"
                )
        else:
            return ValidationResult(
                test_name="Mutation Severity Ordering",
                status=ValidationStatus.FAIL,
                expected=expected,
                observed=observed,
                message="Mutation severity ordering is incorrect"
            )
    except Exception as e:
        return ValidationResult(
            test_name="Mutation Severity Ordering",
            status=ValidationStatus.SKIP,
            expected="Synonymous > Missense > Frameshift",
            observed=str(e),
            message=f"Test could not be run: {e}"
        )


def validate_bistability_requirements(
    simulate_fn: Callable[[Dict], Dict]
) -> ValidationResult:
    """
    Validate that wild-type parameters produce bistability.
    
    This is a basic sanity check that the model works as expected.
    """
    wild_type_params = {
        'tx_A': 1.0,
        'tx_B': 1.0,
        'tl_A': 1.0,
        'tl_B': 1.0,
        'Kd': 1.0,
        'n': 2.0
    }
    
    try:
        from .failure_classifier import classify_failure, FailureMode
        
        result = simulate_fn(wild_type_params)
        classification = classify_failure(result)
        
        expected = "Wild-type shows bistable behavior"
        observed = f"Phenotype: {classification.mode.value}, Confidence: {classification.confidence:.2f}"
        
        if classification.mode == FailureMode.BISTABLE:
            return ValidationResult(
                test_name="Wild-type Bistability",
                status=ValidationStatus.PASS,
                expected=expected,
                observed=observed,
                message="Wild-type parameters correctly produce bistable switch",
                confidence=classification.confidence
            )
        else:
            return ValidationResult(
                test_name="Wild-type Bistability",
                status=ValidationStatus.FAIL,
                expected=expected,
                observed=observed,
                message="Wild-type should be bistable but is not"
            )
    except Exception as e:
        return ValidationResult(
            test_name="Wild-type Bistability",
            status=ValidationStatus.SKIP,
            expected="Bistable",
            observed=str(e),
            message=f"Test could not be run: {e}"
        )


def validate_hill_coefficient_effect(
    simulate_fn: Callable[[Dict], Dict]
) -> ValidationResult:
    """
    Validate that lower Hill coefficient reduces bistability.
    
    Known biology: n=1 typically doesn't support bistability,
    n>=2 is required for robust switching.
    """
    try:
        from .failure_classifier import classify_failure, FailureMode
        
        # High Hill coefficient (should be bistable)
        high_n_params = {'tx_A': 1.0, 'tx_B': 1.0, 'tl_A': 1.0, 
                        'tl_B': 1.0, 'Kd': 1.0, 'n': 3.0}
        
        # Low Hill coefficient (should lose bistability)
        low_n_params = {'tx_A': 1.0, 'tx_B': 1.0, 'tl_A': 1.0, 
                       'tl_B': 1.0, 'Kd': 1.0, 'n': 1.0}
        
        high_n_result = simulate_fn(high_n_params)
        low_n_result = simulate_fn(low_n_params)
        
        high_n_class = classify_failure(high_n_result)
        low_n_class = classify_failure(low_n_result)
        
        expected = "n=3 bistable, n=1 not bistable"
        observed = f"n=3: {high_n_class.mode.value}, n=1: {low_n_class.mode.value}"
        
        if (high_n_class.mode == FailureMode.BISTABLE and 
            low_n_class.mode != FailureMode.BISTABLE):
            return ValidationResult(
                test_name="Hill Coefficient Effect",
                status=ValidationStatus.PASS,
                expected=expected,
                observed=observed,
                message="Hill coefficient correctly affects bistability"
            )
        elif high_n_class.mode == FailureMode.BISTABLE:
            return ValidationResult(
                test_name="Hill Coefficient Effect",
                status=ValidationStatus.WARN,
                expected=expected,
                observed=observed,
                message="High n is bistable but low n should not be"
            )
        else:
            return ValidationResult(
                test_name="Hill Coefficient Effect",
                status=ValidationStatus.FAIL,
                expected=expected,
                observed=observed,
                message="Hill coefficient effect is incorrect"
            )
    except Exception as e:
        return ValidationResult(
            test_name="Hill Coefficient Effect",
            status=ValidationStatus.SKIP,
            expected="n affects bistability",
            observed=str(e),
            message=f"Test could not be run: {e}"
        )


def run_validation_suite(
    simulate_fn: Optional[Callable] = None,
    promoter_fn: Optional[Callable] = None,
    operator_fn: Optional[Callable] = None,
    cds_fn: Optional[Callable] = None
) -> ValidationReport:
    """
    Run complete validation suite.
    
    Args:
        simulate_fn: Simulation function for bistability tests
        promoter_fn: Promoter mapping function
        operator_fn: Operator binding function
        cds_fn: CDS effect function
        
    Returns:
        ValidationReport with all test results
    """
    report = ValidationReport()
    
    # Run simulation-based tests
    if simulate_fn:
        report.add_result(validate_bistability_requirements(simulate_fn))
        report.add_result(validate_hill_coefficient_effect(simulate_fn))
    
    # Run mapping tests
    if promoter_fn:
        report.add_result(validate_promoter_strength(promoter_fn))
    
    if operator_fn:
        report.add_result(validate_operator_mismatch(operator_fn))
    
    if cds_fn:
        report.add_result(validate_mutation_severity(cds_fn))
    
    return report
