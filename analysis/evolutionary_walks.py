"""
Evolutionary Walks Module

Simulates multi-mutation trajectories to understand how robustness
decays as mutations accumulate. This provides:
- Evolutionary resilience assessment
- Robustness decay curves
- Region-specific mutation tolerance

Public API:
    simulate_mutation_walk(sequence, regions, max_mutations) -> MutationTrajectory
    compute_decay_curve(trajectories) -> DecayCurve
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Callable, Tuple
import warnings
import copy

from .failure_classifier import FailureMode, classify_failure


@dataclass
class MutationStep:
    """
    Single step in a mutation walk.
    
    Attributes:
        step_number: Number of accumulated mutations
        mutation_info: Info about the mutation added at this step
        region_affected: Which circuit region was mutated
        phenotype: Resulting phenotype after this mutation
        is_functional: Whether circuit still exhibits bistability
        robustness_score: Continuous robustness measure
        parameters: Current parameter values after mutation
    """
    step_number: int
    mutation_info: Dict[str, Any]
    region_affected: str
    phenotype: str
    is_functional: bool
    robustness_score: float
    parameters: Dict[str, float]
    confidence: float = 0.0


@dataclass
class MutationTrajectory:
    """
    Complete mutation walk from wild-type to failure.
    
    Attributes:
        steps: List of MutationStep objects
        mutations_to_failure: Number of mutations until bistability lost
        final_phenotype: Final phenotype at trajectory end
        total_mutations: Total mutations in trajectory
        region_breakdown: Dict of mutations per region
    """
    steps: List[MutationStep] = field(default_factory=list)
    mutations_to_failure: int = 0
    final_phenotype: str = ""
    total_mutations: int = 0
    region_breakdown: Dict[str, int] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            'mutations_to_failure': self.mutations_to_failure,
            'final_phenotype': self.final_phenotype,
            'total_mutations': self.total_mutations,
            'region_breakdown': self.region_breakdown,
            'steps': [
                {
                    'step': s.step_number,
                    'region': s.region_affected,
                    'phenotype': s.phenotype,
                    'is_functional': s.is_functional,
                    'robustness': s.robustness_score
                }
                for s in self.steps
            ],
            'robustness_curve': [s.robustness_score for s in self.steps]
        }


@dataclass 
class DecayCurve:
    """
    Aggregated robustness decay statistics across multiple trajectories.
    
    Attributes:
        mutation_counts: List of mutation counts (x-axis)
        mean_robustness: Mean robustness at each mutation count
        std_robustness: Standard deviation at each mutation count
        pct_functional: Percentage still functional at each count
        n_trajectories: Number of trajectories averaged
    """
    mutation_counts: List[int] = field(default_factory=list)
    mean_robustness: List[float] = field(default_factory=list)
    std_robustness: List[float] = field(default_factory=list)
    pct_functional: List[float] = field(default_factory=list)
    n_trajectories: int = 0
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'mutation_counts': self.mutation_counts,
            'mean_robustness': self.mean_robustness,
            'std_robustness': self.std_robustness,
            'pct_functional': self.pct_functional,
            'n_trajectories': self.n_trajectories
        }


def simulate_mutation_walk(
    mutate_fn: Callable,
    simulate_fn: Callable[[Dict[str, float]], Dict[str, Any]],
    map_mutations_fn: Callable,
    base_sequence: str,
    regions: List,
    base_params: Dict[str, float],
    max_mutations: int = 10,
    stop_on_failure: bool = True,
    seed: Optional[int] = None
) -> MutationTrajectory:
    """
    Simulate sequential mutation accumulation.
    
    Starts with wild-type and adds mutations one at a time,
    tracking robustness decay until bistability is lost.
    
    Args:
        mutate_fn: Function to generate a single mutation
        simulate_fn: Function to run simulation with parameters
        map_mutations_fn: Function to map mutations to parameters
        base_sequence: Wild-type DNA sequence
        regions: List of annotated regions
        base_params: Wild-type parameters
        max_mutations: Maximum mutations to accumulate
        stop_on_failure: If True, stop when bistability is lost
        seed: Random seed
        
    Returns:
        MutationTrajectory with all steps recorded
    """
    if seed is not None:
        np.random.seed(seed)
    
    trajectory = MutationTrajectory()
    current_sequence = base_sequence
    current_params = base_params.copy()
    accumulated_mutations = []
    region_counts = {}
    
    # Record initial (wild-type) state
    try:
        wt_result = simulate_fn(current_params)
        wt_classification = classify_failure(wt_result)
        
        initial_step = MutationStep(
            step_number=0,
            mutation_info={'type': 'wild_type'},
            region_affected='none',
            phenotype=wt_classification.mode.value,
            is_functional=(wt_classification.mode == FailureMode.BISTABLE),
            robustness_score=1.0 if wt_classification.mode == FailureMode.BISTABLE else 0.0,
            parameters=current_params.copy(),
            confidence=wt_classification.confidence
        )
        trajectory.steps.append(initial_step)
    except Exception as e:
        warnings.warn(f"Wild-type simulation failed: {e}")
        return trajectory
    
    # Accumulate mutations
    for step in range(1, max_mutations + 1):
        try:
            # Generate one new mutation
            mutated_seq, new_mutations = mutate_fn(
                current_sequence, 
                regions, 
                n_mut=1
            )
            
            if not new_mutations:
                continue
            
            mutation = new_mutations[0]
            accumulated_mutations.append(mutation)
            
            # Track region
            region_name = getattr(mutation, 'region_name', 'unknown')
            region_counts[region_name] = region_counts.get(region_name, 0) + 1
            
            # Map all accumulated mutations to parameters
            new_params = map_mutations_fn(
                accumulated_mutations,
                mutated_seq,
                base_sequence,
                regions,
                base_params
            )
            
            # Run simulation
            sim_result = simulate_fn(new_params)
            classification = classify_failure(sim_result)
            
            is_functional = (classification.mode == FailureMode.BISTABLE)
            
            # Calculate robustness score based on classification confidence
            if is_functional:
                robustness = classification.confidence
            else:
                robustness = 0.0
            
            step_record = MutationStep(
                step_number=step,
                mutation_info={
                    'position': getattr(mutation, 'position', 0),
                    'type': getattr(mutation, 'mutation_type', 'unknown'),
                    'original': getattr(mutation, 'original_base', ''),
                    'mutated': getattr(mutation, 'mutated_base', '')
                },
                region_affected=region_name,
                phenotype=classification.mode.value,
                is_functional=is_functional,
                robustness_score=robustness,
                parameters=new_params.copy(),
                confidence=classification.confidence
            )
            trajectory.steps.append(step_record)
            
            # Update current state
            current_sequence = mutated_seq
            current_params = new_params
            
            # Track first failure
            if not is_functional and trajectory.mutations_to_failure == 0:
                trajectory.mutations_to_failure = step
                
                if stop_on_failure:
                    break
                    
        except Exception as e:
            warnings.warn(f"Mutation step {step} failed: {e}")
            continue
    
    # Finalize trajectory
    if trajectory.steps:
        trajectory.final_phenotype = trajectory.steps[-1].phenotype
        trajectory.total_mutations = len(trajectory.steps) - 1  # Exclude wild-type
        
    trajectory.region_breakdown = region_counts
    
    if trajectory.mutations_to_failure == 0 and trajectory.steps:
        # Never lost bistability
        trajectory.mutations_to_failure = trajectory.total_mutations
    
    return trajectory


def simulate_multiple_walks(
    mutate_fn: Callable,
    simulate_fn: Callable,
    map_mutations_fn: Callable,
    base_sequence: str,
    regions: List,
    base_params: Dict[str, float],
    n_walks: int = 10,
    max_mutations: int = 10,
    seed: Optional[int] = None
) -> List[MutationTrajectory]:
    """
    Run multiple mutation walk simulations.
    
    Returns list of trajectories for statistical analysis.
    """
    if seed is not None:
        np.random.seed(seed)
    
    trajectories = []
    
    for i in range(n_walks):
        trajectory = simulate_mutation_walk(
            mutate_fn=mutate_fn,
            simulate_fn=simulate_fn,
            map_mutations_fn=map_mutations_fn,
            base_sequence=base_sequence,
            regions=regions,
            base_params=base_params,
            max_mutations=max_mutations,
            stop_on_failure=False,  # Continue for full curve
            seed=None  # Use overall seed control
        )
        trajectories.append(trajectory)
    
    return trajectories


def compute_decay_curve(
    trajectories: List[MutationTrajectory]
) -> DecayCurve:
    """
    Compute mean robustness decay curve from multiple trajectories.
    
    Args:
        trajectories: List of MutationTrajectory objects
        
    Returns:
        DecayCurve with aggregated statistics
    """
    if not trajectories:
        return DecayCurve()
    
    # Find maximum mutation count
    max_steps = max(len(t.steps) for t in trajectories)
    
    # Aggregate robustness at each step
    robustness_by_step = [[] for _ in range(max_steps)]
    functional_by_step = [[] for _ in range(max_steps)]
    
    for trajectory in trajectories:
        for step in trajectory.steps:
            if step.step_number < max_steps:
                robustness_by_step[step.step_number].append(step.robustness_score)
                functional_by_step[step.step_number].append(1.0 if step.is_functional else 0.0)
    
    # Compute statistics
    mutation_counts = []
    mean_robustness = []
    std_robustness = []
    pct_functional = []
    
    for i in range(max_steps):
        if robustness_by_step[i]:
            mutation_counts.append(i)
            mean_robustness.append(float(np.mean(robustness_by_step[i])))
            std_robustness.append(float(np.std(robustness_by_step[i])))
            pct_functional.append(float(np.mean(functional_by_step[i]) * 100))
    
    return DecayCurve(
        mutation_counts=mutation_counts,
        mean_robustness=mean_robustness,
        std_robustness=std_robustness,
        pct_functional=pct_functional,
        n_trajectories=len(trajectories)
    )


def analyze_region_vulnerability(
    trajectories: List[MutationTrajectory]
) -> Dict[str, Any]:
    """
    Analyze which regions are most vulnerable to mutations.
    
    Returns breakdown of failures by region and average mutations
    tolerated per region.
    """
    region_failures = {}
    region_mutations = {}
    
    for trajectory in trajectories:
        for step in trajectory.steps:
            region = step.region_affected
            if region == 'none':
                continue
                
            if region not in region_mutations:
                region_mutations[region] = []
                region_failures[region] = 0
            
            region_mutations[region].append(step.step_number)
            
            # Count failures
            if not step.is_functional and step.step_number > 0:
                # Check if this was the failure point
                prev_step = trajectory.steps[step.step_number - 1] if step.step_number > 0 else None
                if prev_step and prev_step.is_functional:
                    region_failures[region] += 1
    
    # Compute vulnerability scores
    vulnerability = {}
    for region in region_mutations:
        n_mutations = len(region_mutations[region])
        n_failures = region_failures.get(region, 0)
        
        vulnerability[region] = {
            'total_mutations': n_mutations,
            'failures_caused': n_failures,
            'failure_rate': n_failures / n_mutations if n_mutations > 0 else 0.0
        }
    
    # Rank by vulnerability
    ranked = sorted(
        vulnerability.items(),
        key=lambda x: x[1]['failure_rate'],
        reverse=True
    )
    
    return {
        'region_vulnerability': vulnerability,
        'ranking': [r[0] for r in ranked],
        'most_vulnerable': ranked[0][0] if ranked else None,
        'least_vulnerable': ranked[-1][0] if ranked else None
    }
