"""
Robustness Metrics Module

Computes robustness metrics from batch mutation simulation results.
These metrics quantify how tolerant the toggle switch circuit is to mutations.

Public API:
    compute_robustness(results) -> dict
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional
from collections import Counter


def compute_robustness(
    results: pd.DataFrame,
    label_column: str = 'failure_label',
    mutation_columns: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    Compute robustness metrics from batch simulation results.
    
    Analyzes the distribution of failure modes and identifies
    vulnerable regions or mutation types.
    
    Args:
        results: DataFrame with simulation results
            Required columns: failure_label
            Optional: mutation_list, region, mutation_type, parameters
        label_column: Name of column containing failure labels
        mutation_columns: Optional list of mutation info columns
    
    Returns:
        Dict containing:
            - 'total_variants': Total number of variants analyzed
            - 'pct_bistable': Percentage of variants that remain bistable
            - 'pct_failed': Percentage of variants that failed
            - 'failure_distribution': Dict of failure mode counts
            - 'failure_percentages': Dict of failure mode percentages
            - 'region_vulnerability': Dict of failure rates by region
            - 'mutation_type_effects': Dict of effects by mutation type
            - 'robustness_score': Overall robustness score (0-1)
            - 'summary_stats': Additional summary statistics
    
    Notes:
        Robustness score = fraction of variants that remain bistable
        Higher score = more robust circuit
    """
    if results.empty:
        return _empty_robustness_result()
    
    total = len(results)
    
    # Count failure modes
    failure_counts = results[label_column].value_counts().to_dict()
    
    # Calculate percentages
    failure_pcts = {
        mode: (count / total) * 100 
        for mode, count in failure_counts.items()
    }
    
    # Calculate key metrics
    bistable_count = failure_counts.get('bistable', 0)
    pct_bistable = (bistable_count / total) * 100 if total > 0 else 0
    pct_failed = 100 - pct_bistable
    
    # Robustness score (0-1, higher = more robust)
    robustness_score = bistable_count / total if total > 0 else 0
    
    # Analyze region vulnerability if data available
    region_vulnerability = {}
    if 'region' in results.columns:
        region_vulnerability = _compute_region_vulnerability(
            results, label_column
        )
    
    # Analyze mutation type effects if data available
    mutation_type_effects = {}
    if 'mutation_type' in results.columns:
        mutation_type_effects = _compute_mutation_type_effects(
            results, label_column
        )
    
    # Summary statistics
    summary_stats = _compute_summary_stats(results, label_column)
    
    return {
        'total_variants': total,
        'pct_bistable': pct_bistable,
        'pct_failed': pct_failed,
        'failure_distribution': failure_counts,
        'failure_percentages': failure_pcts,
        'region_vulnerability': region_vulnerability,
        'mutation_type_effects': mutation_type_effects,
        'robustness_score': robustness_score,
        'summary_stats': summary_stats
    }


def _empty_robustness_result() -> Dict[str, Any]:
    """Return empty robustness result for edge cases."""
    return {
        'total_variants': 0,
        'pct_bistable': 0.0,
        'pct_failed': 0.0,
        'failure_distribution': {},
        'failure_percentages': {},
        'region_vulnerability': {},
        'mutation_type_effects': {},
        'robustness_score': 0.0,
        'summary_stats': {}
    }


def _compute_region_vulnerability(
    results: pd.DataFrame,
    label_column: str
) -> Dict[str, Dict[str, float]]:
    """
    Compute failure rates for each genomic region.
    
    Returns dict mapping region -> {failure_rate, dominant_failure_mode}
    """
    region_stats = {}
    
    for region in results['region'].dropna().unique():
        region_data = results[results['region'] == region]
        total_region = len(region_data)
        
        if total_region == 0:
            continue
        
        # Count failures in this region
        failure_counts = region_data[label_column].value_counts()
        bistable_count = failure_counts.get('bistable', 0)
        
        failure_rate = 1 - (bistable_count / total_region)
        
        # Find dominant failure mode
        non_bistable = failure_counts.drop('bistable', errors='ignore')
        dominant_mode = non_bistable.idxmax() if len(non_bistable) > 0 else 'none'
        
        region_stats[region] = {
            'failure_rate': failure_rate,
            'total_mutations': total_region,
            'dominant_failure_mode': dominant_mode,
            'failure_breakdown': failure_counts.to_dict()
        }
    
    return region_stats


def _compute_mutation_type_effects(
    results: pd.DataFrame,
    label_column: str
) -> Dict[str, Dict[str, float]]:
    """
    Compute failure rates for each mutation type.
    """
    type_stats = {}
    
    for mut_type in results['mutation_type'].dropna().unique():
        type_data = results[results['mutation_type'] == mut_type]
        total_type = len(type_data)
        
        if total_type == 0:
            continue
        
        failure_counts = type_data[label_column].value_counts()
        bistable_count = failure_counts.get('bistable', 0)
        
        failure_rate = 1 - (bistable_count / total_type)
        
        type_stats[str(mut_type)] = {
            'failure_rate': failure_rate,
            'total_mutations': total_type,
            'failure_breakdown': failure_counts.to_dict()
        }
    
    return type_stats


def _compute_summary_stats(
    results: pd.DataFrame,
    label_column: str
) -> Dict[str, Any]:
    """
    Compute summary statistics from results.
    """
    stats = {}
    
    # Mean number of mutations per variant
    if 'n_mutations' in results.columns:
        stats['mean_mutations_per_variant'] = results['n_mutations'].mean()
        stats['std_mutations_per_variant'] = results['n_mutations'].std()
    
    # Parameter statistics for failed vs successful variants
    param_cols = [col for col in results.columns if col.startswith('param_')]
    
    if param_cols and label_column in results.columns:
        bistable_mask = results[label_column] == 'bistable'
        
        for param in param_cols:
            if param in results.columns:
                bistable_vals = results.loc[bistable_mask, param]
                failed_vals = results.loc[~bistable_mask, param]
                
                stats[f'{param}_bistable_mean'] = bistable_vals.mean() if len(bistable_vals) > 0 else np.nan
                stats[f'{param}_failed_mean'] = failed_vals.mean() if len(failed_vals) > 0 else np.nan
    
    return stats


def compute_position_failure_rates(
    results: pd.DataFrame,
    position_column: str = 'mutation_positions',
    label_column: str = 'failure_label',
    seq_length: int = None
) -> np.ndarray:
    """
    Compute failure rate at each sequence position.
    
    Used for generating heatmaps of vulnerable regions.
    
    Args:
        results: DataFrame with mutation data
        position_column: Column containing mutation positions
        label_column: Column containing failure labels
        seq_length: Length of sequence (inferred if not provided)
    
    Returns:
        Array of failure rates per position
    """
    if position_column not in results.columns:
        return np.array([])
    
    # Infer sequence length if not provided
    if seq_length is None:
        all_positions = []
        for pos_list in results[position_column].dropna():
            if isinstance(pos_list, (list, tuple)):
                all_positions.extend(pos_list)
            elif isinstance(pos_list, str):
                # Parse string representation
                try:
                    positions = eval(pos_list)
                    if isinstance(positions, (list, tuple)):
                        all_positions.extend(positions)
                except:
                    pass
        
        seq_length = max(all_positions) + 1 if all_positions else 0
    
    if seq_length == 0:
        return np.array([])
    
    # Count mutations and failures at each position
    mutation_counts = np.zeros(seq_length)
    failure_counts = np.zeros(seq_length)
    
    for idx, row in results.iterrows():
        positions = row.get(position_column, [])
        label = row.get(label_column, '')
        
        if isinstance(positions, str):
            try:
                positions = eval(positions)
            except:
                continue
        
        if not isinstance(positions, (list, tuple)):
            continue
        
        for pos in positions:
            if 0 <= pos < seq_length:
                mutation_counts[pos] += 1
                if label != 'bistable':
                    failure_counts[pos] += 1
    
    # Calculate failure rates
    with np.errstate(divide='ignore', invalid='ignore'):
        failure_rates = np.where(
            mutation_counts > 0,
            failure_counts / mutation_counts,
            0
        )
    
    return failure_rates


def compute_region_failure_matrix(
    results: pd.DataFrame,
    regions: List[Dict],
    label_column: str = 'failure_label'
) -> Dict[str, np.ndarray]:
    """
    Compute failure rate matrix by region and failure mode.
    
    Args:
        results: DataFrame with simulation results
        regions: List of region dicts with 'name', 'start', 'end'
        label_column: Column containing failure labels
    
    Returns:
        Dict with 'matrix' (regions x failure_modes) and 'labels'
    """
    failure_modes = [
        'bistable', 'loss_of_bistability', 'leaky', 
        'no_expression', 'oscillatory', 'simulation_failed'
    ]
    
    region_names = [r['name'] if isinstance(r, dict) else r.name for r in regions]
    n_regions = len(region_names)
    n_modes = len(failure_modes)
    
    matrix = np.zeros((n_regions, n_modes))
    
    if 'region' in results.columns:
        for i, region_name in enumerate(region_names):
            region_data = results[results['region'] == region_name]
            total = len(region_data)
            
            if total > 0:
                for j, mode in enumerate(failure_modes):
                    count = len(region_data[region_data[label_column] == mode])
                    matrix[i, j] = count / total
    
    return {
        'matrix': matrix,
        'region_labels': region_names,
        'mode_labels': failure_modes
    }
