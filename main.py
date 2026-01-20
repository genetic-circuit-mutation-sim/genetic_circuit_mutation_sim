#!/usr/bin/env python3
"""
Genetic Circuit Mutation Simulator - Main CLI

This is the main entry point for running Monte Carlo mutation simulations
on the toggle switch genetic circuit.

The pipeline:
1. Generate N mutant sequences from the wild-type toggle switch
2. Map sequence mutations to parameter perturbations
3. Simulate each mutant using deterministic ODE (or optional stochastic)
4. Classify failure modes
5. Output results.csv, failure_heatmap.png, and robustness_summary.json

Usage:
    python main.py --n 1000
    python main.py --n 100 --stochastic --output-dir ./results

Author: Genetic Circuit Mutation Simulator
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Any, Optional
import warnings

import numpy as np
import pandas as pd

# Suppress matplotlib backend warnings
warnings.filterwarnings('ignore', category=UserWarning)

# Import package modules
from mutation_engine import mutate_sequence, AnnotatedRegion, MutationType
from mutation_engine.mutation_types import RegionType
from mutation_engine.sequence_mutator import generate_toggle_switch_sequence

from mapping import promoter_to_tx_rate, rbs_to_tl_rate, cds_to_protein_effect
from mapping.cds_mapping import operator_to_binding_params

from simulation.deterministic import (
    simulate_deterministic, 
    get_default_model,
    create_mutant_model
)
from simulation.stochastic import simulate_gillespie

from analysis import (
    classify_failure, 
    FailureMode, 
    compute_robustness,
    classify_probabilistic,
    simulate_mutation_walk,
    run_validation_suite
)
from analysis.robustness_metrics import compute_position_failure_rates
from analysis.sensitivity import sweep_classification_thresholds, generate_sensitivity_summary

from visualization.heatmaps import (
    plot_failure_heatmap,
    plot_region_failure_heatmap,
    plot_failure_mode_pie,
    create_summary_figure
)


# Default parameter values for wild-type toggle switch
DEFAULT_PARAMS = {
    'tx_A': 1.0,
    'tx_B': 1.0,
    'tl_A': 1.0,
    'tl_B': 1.0,
    'Kd': 1.0,
    'n': 2.0
}

# Parameter bounds (biologically plausible ranges)
PARAM_BOUNDS = {
    'tx_A': (1e-4, 1e2),
    'tx_B': (1e-4, 1e2),
    'tl_A': (1e-4, 1e2),
    'tl_B': (1e-4, 1e2),
    'Kd': (1e-3, 1e3),
    'n': (1.0, 6.0)
}


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Run Monte Carlo mutation simulation on toggle switch circuit',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python main.py --n 1000
    python main.py --n 500 --mutations-per-variant 3
    python main.py --n 100 --stochastic --stochastic-runs 50
    python main.py --n 1000 --output-dir ./my_results --seed 42
        """
    )
    
    parser.add_argument(
        '--n', type=int, default=1000,
        help='Number of mutant variants to generate (default: 1000)'
    )
    
    parser.add_argument(
        '--mutations-per-variant', type=int, default=1,
        help='Number of mutations per variant (default: 1)'
    )
    
    parser.add_argument(
        '--stochastic', action='store_true',
        help='Use stochastic Gillespie simulation instead of deterministic ODE'
    )
    
    parser.add_argument(
        '--stochastic-runs', type=int, default=100,
        help='Number of stochastic runs per variant (default: 100)'
    )
    
    parser.add_argument(
        '--output-dir', type=str, default='.',
        help='Output directory for results (default: current directory)'
    )
    
    parser.add_argument(
        '--seed', type=int, default=None,
        help='Random seed for reproducibility'
    )
    
    parser.add_argument(
        '--substitution-rate', type=float, default=0.7,
        help='Relative rate of substitution mutations (default: 0.7)'
    )
    
    parser.add_argument(
        '--insertion-rate', type=float, default=0.15,
        help='Relative rate of insertion mutations (default: 0.15)'
    )
    
    parser.add_argument(
        '--deletion-rate', type=float, default=0.15,
        help='Relative rate of deletion mutations (default: 0.15)'
    )
    
    parser.add_argument(
        '--quiet', action='store_true',
        help='Suppress progress output'
    )

    parser.add_argument(
        '--probabilistic', action='store_true',
        help='Use ensemble-based probabilistic classification'
    )

    parser.add_argument(
        '--ensemble-runs', type=int, default=20,
        help='Number of runs for ensemble classification (default: 20)'
    )

    parser.add_argument(
        '--multi-mutation', action='store_true',
        help='Run evolutionary walk (multi-mutation trajectory) analysis'
    )

    parser.add_argument(
        '--sensitivity', action='store_true',
        help='Run threshold sensitivity analysis'
    )

    # Base Parameter Overrides
    parser.add_argument('--base-tx', type=float, default=1.0, help='Base transcription rate')
    parser.add_argument('--base-tl', type=float, default=1.0, help='Base translation rate')
    parser.add_argument('--base-Kd', type=float, default=1.0, help='Base dissociation constant')
    parser.add_argument('--base-n', type=float, default=2.0, help='Base Hill coefficient')
    
    return parser.parse_args()


def map_mutations_to_parameters(
    mutations: List,
    mutated_seq: str,
    original_seq: str,
    regions: List[AnnotatedRegion],
    base_params: Dict[str, float]
) -> Dict[str, float]:
    """
    Map sequence mutations to parameter perturbations.
    
    For each mutation:
    - Promoter mutations: affect tx_* (transcription)
    - RBS mutations: affect tl_* (translation)
    - CDS mutations: affect protein function (can reduce tx/tl to 0)
    - Operator mutations: affect Kd and n (repressor binding)
    
    Args:
        mutations: List of Mutation objects
        mutated_seq: Full mutated sequence
        original_seq: Original wild-type sequence
        regions: List of AnnotatedRegion objects
        base_params: Base (wild-type) parameter dict
    
    Returns:
        Dict of perturbed parameters
    """
    params = base_params.copy()
    
    # Build region lookup
    region_lookup = {r.name: r for r in regions}
    
    for mutation in mutations:
        if mutation.region is None:
            continue
        
        region = mutation.region
        region_type = region.region_type
        gene = region.gene  # 'A' or 'B'
        
        # Extract region sequence from mutated sequence
        region_seq = mutated_seq[region.start:region.end]
        original_region_seq = original_seq[region.start:region.end]
        
        if region_type == RegionType.PROMOTER:
            # Promoter mutation -> affects transcription rate
            multiplier = promoter_to_tx_rate(region_seq, add_noise=True)
            
            param_name = f'tx_{gene}' if gene else 'tx_A'
            params[param_name] = params[param_name] * multiplier
            
        elif region_type == RegionType.RBS:
            # RBS mutation -> affects translation rate
            multiplier = rbs_to_tl_rate(region_seq, add_noise=True)
            
            param_name = f'tl_{gene}' if gene else 'tl_A'
            params[param_name] = params[param_name] * multiplier
            
        elif region_type == RegionType.CDS:
            # CDS mutation -> can cause nonsense/frameshift
            effect = cds_to_protein_effect(
                region_seq, 
                original_seq=original_region_seq,
                add_noise=True
            )
            
            effect_mult = effect['effect_multiplier']
            
            # CDS damage reduces both transcription and translation
            # (modeling mRNA decay triggered by premature stop)
            param_tx = f'tx_{gene}' if gene else 'tx_A'
            param_tl = f'tl_{gene}' if gene else 'tl_A'
            
            params[param_tx] = params[param_tx] * effect_mult
            params[param_tl] = params[param_tl] * effect_mult
            
        elif region_type == RegionType.OPERATOR:
            # Operator mutation -> affects repressor binding
            binding = operator_to_binding_params(
                region_seq,
                original_seq=original_region_seq,
                add_noise=True
            )
            
            params['Kd'] = params['Kd'] * binding['kd_multiplier']
            params['n'] = params['n'] * binding['n_multiplier']
    
    # Clamp parameters to biologically plausible bounds
    for param_name, (min_val, max_val) in PARAM_BOUNDS.items():
        if param_name in params:
            params[param_name] = np.clip(params[param_name], min_val, max_val)
    
    return params


def run_monte_carlo(
    n_variants: int,
    n_mutations: int,
    mutation_rates: Dict[str, float],
    base_params: Dict[str, float],
    use_stochastic: bool,
    stochastic_runs: int,
    seed: Optional[int],
    quiet: bool,
    args: Any = None
) -> pd.DataFrame:
    """
    Run Monte Carlo mutation simulation.
    
    Args:
        n_variants: Number of mutant variants to generate
        n_mutations: Mutations per variant
        mutation_rates: Dict with substitution/insertion/deletion rates
        base_params: Base parameter values
        use_stochastic: Whether to use stochastic simulation
        stochastic_runs: Number of runs for stochastic sim
        seed: Random seed
        quiet: Suppress progress output
    
    Returns:
        DataFrame with simulation results
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Generate wild-type sequence and regions
    wt_sequence, regions = generate_toggle_switch_sequence()
    
    # Get base model
    base_model = get_default_model()
    
    # Storage for results
    results = []
    
    start_time = time.time()
    
    for i in range(n_variants):
        if not quiet and (i + 1) % 100 == 0:
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed
            eta = (n_variants - i - 1) / rate
            print(f"Progress: {i+1}/{n_variants} ({100*(i+1)/n_variants:.1f}%) "
                  f"- ETA: {eta:.1f}s")
        
        try:
            # Generate mutant
            mutated_seq, mutations = mutate_sequence(
                wt_sequence,
                regions,
                n_mut=n_mutations,
                rates=mutation_rates,
                seed=seed + i if seed else None
            )
            
            # Map to parameters
            perturbed_params = map_mutations_to_parameters(
                mutations,
                mutated_seq,
                wt_sequence,
                regions,
                base_params
            )
            
            # Run simulation
            if use_stochastic:
                sim_result = simulate_gillespie(
                    base_model,
                    n_runs=stochastic_runs,
                    parameters=perturbed_params
                )
                # Convert stochastic result to deterministic-like format
                sim_result['simulation_successful'] = sim_result.get('simulation_successful', False)
                sim_result['steady_state_A_highA'] = np.mean(sim_result.get('final_A_distribution', [0]))
                sim_result['steady_state_B_highA'] = np.mean(sim_result.get('final_B_distribution', [0]))
                sim_result['steady_state_A_highB'] = sim_result['steady_state_A_highA']
                sim_result['steady_state_B_highB'] = sim_result['steady_state_B_highA']
            else:
                sim_result = simulate_deterministic(
                    base_model,
                    parameters=perturbed_params
                )
            
            # Classify failure mode
            if hasattr(args, 'probabilistic') and args.probabilistic:
                prob_result = classify_probabilistic(
                    base_model,
                    perturbed_params,
                    n_runs=args.ensemble_runs
                )
                classification = prob_result['most_likely_phenotype']
                # Store probabilistic details in sim_result for later
                sim_result['probabilities'] = prob_result['distribution']
            else:
                classification = classify_failure(sim_result)
            
            # Build result row
            mutation_strs = [str(m) for m in mutations]
            mutation_positions = [m.position for m in mutations]
            mutation_types = [m.mutation_type.value for m in mutations]
            first_region = mutations[0].region.name if mutations and mutations[0].region else ''
            first_type = mutations[0].mutation_type.value if mutations else ''
            
            result_row = {
                'variant_id': i,
                'n_mutations': len(mutations),
                'mutations': ';'.join(mutation_strs),
                'mutation_positions': mutation_positions,
                'mutation_type': first_type,
                'region': first_region,
                'mutated_sequence': mutated_seq[:100] + '...',  # Truncate for storage
                
                # Parameters
                'param_tx_A': perturbed_params['tx_A'],
                'param_tx_B': perturbed_params['tx_B'],
                'param_tl_A': perturbed_params['tl_A'],
                'param_tl_B': perturbed_params['tl_B'],
                'param_Kd': perturbed_params['Kd'],
                'param_n': perturbed_params['n'],
                
                # Steady states
                'ss_A_highA': sim_result.get('steady_state_A_highA', 0),
                'ss_B_highA': sim_result.get('steady_state_B_highA', 0),
                'ss_A_highB': sim_result.get('steady_state_A_highB', 0),
                'ss_B_highB': sim_result.get('steady_state_B_highB', 0),
                
                # Classification
                'failure_label': classification.mode.value,
                'classification_confidence': classification.confidence,
                
                # Simulation status
                'simulation_successful': sim_result.get('simulation_successful', False)
            }
            
            results.append(result_row)
            
        except Exception as e:
            # Record failed simulation
            results.append({
                'variant_id': i,
                'n_mutations': 0,
                'mutations': '',
                'mutation_positions': [],
                'mutation_type': '',
                'region': '',
                'mutated_sequence': '',
                'param_tx_A': np.nan,
                'param_tx_B': np.nan,
                'param_tl_A': np.nan,
                'param_tl_B': np.nan,
                'param_Kd': np.nan,
                'param_n': np.nan,
                'ss_A_highA': 0,
                'ss_B_highA': 0,
                'ss_A_highB': 0,
                'ss_B_highB': 0,
                'failure_label': 'simulation_failed',
                'classification_confidence': 1.0,
                'simulation_successful': False,
                'error': str(e)
            })
    
    return pd.DataFrame(results)


def main():
    """Main entry point."""
    args = parse_args()
    
    print("=" * 60)
    print("GENETIC CIRCUIT MUTATION SIMULATOR")
    print("=" * 60)
    print(f"Variants to generate: {args.n}")
    print(f"Mutations per variant: {args.mutations_per_variant}")
    print(f"Simulation mode: {'Stochastic' if args.stochastic else 'Deterministic'}")
    print(f"Output directory: {args.output_dir}")
    if args.seed:
        print(f"Random seed: {args.seed}")
    print("=" * 60)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Mutation rates
    mutation_rates = {
        'substitution_rate': args.substitution_rate,
        'insertion_rate': args.insertion_rate,
        'deletion_rate': args.deletion_rate
    }
    
    # Run Monte Carlo simulation
    print("\nRunning Monte Carlo simulation...")
    start_time = time.time()
    
    # Base parameters
    base_params = {
        'tx_A': args.base_tx,
        'tx_B': args.base_tx,
        'tl_A': args.base_tl,
        'tl_B': args.base_tl,
        'Kd': args.base_Kd,
        'n': args.base_n
    }
    
    # Global args for run_monte_carlo to access probabilistic settings
    global_args = args
    
    if args.multi_mutation:
        print("\nRunning multi-mutation trajectory analysis...")
        wt_sequence, regions = generate_toggle_switch_sequence()
        model = get_default_model()
        walk_result = simulate_mutation_walk(
            wt_sequence,
            regions,
            model,
            max_mutations=10,
            n_walks=args.n // 10 if args.n > 100 else 10
        )
        # Save walk results
        with open(output_dir / 'evolutionary_walks.json', 'w') as f:
            json.dump(walk_result, f, indent=2)
        print(f"Saved evolutionary walks to: {output_dir / 'evolutionary_walks.json'}")

    if args.sensitivity:
        print("\nRunning threshold sensitivity analysis...")
        # We need some simulation results first
        # For simplicity, if --sensitivity is on, we run a smaller sweep
        pass # Will be handled after simulation
    
    results_df = run_monte_carlo(
        n_variants=args.n,
        n_mutations=args.mutations_per_variant,
        mutation_rates=mutation_rates,
        base_params=base_params,
        use_stochastic=args.stochastic,
        stochastic_runs=args.stochastic_runs,
        seed=args.seed,
        quiet=args.quiet,
        args=args # Pass args to run_monte_carlo
    )
    
    elapsed = time.time() - start_time
    print(f"\nSimulation completed in {elapsed:.2f} seconds")
    print(f"Rate: {args.n / elapsed:.1f} variants/second")
    
    # Save results CSV
    results_path = output_dir / 'results.csv'
    results_df.to_csv(results_path, index=False)
    print(f"\nSaved results to: {results_path}")
    
    # Compute robustness metrics
    print("\nComputing robustness metrics...")
    robustness = compute_robustness(results_df, label_column='failure_label')
    
    # Save robustness summary
    robustness_path = output_dir / 'robustness_summary.json'
    
    # Convert numpy types for JSON serialization
    def convert_for_json(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_for_json(i) for i in obj]
        return obj
    
    robustness_json = convert_for_json(robustness)
    
    with open(robustness_path, 'w') as f:
        json.dump(robustness_json, f, indent=2)
    
    print(f"Saved robustness summary to: {robustness_path}")
    
    # Sensitivity analysis summary
    if args.sensitivity:
        print("\nAnalyzing threshold sensitivity...")
        sweep_data = sweep_classification_thresholds(results_df.to_dict('records'))
        summary = generate_sensitivity_summary(sweep_data)
        
        sensitivity_path = output_dir / 'sensitivity_analysis.json'
        with open(sensitivity_path, 'w') as f:
            json.dump(convert_for_json(summary), f, indent=2)
        print(f"Saved sensitivity analysis to: {sensitivity_path}")

    
    # Print summary
    print("\n" + "=" * 60)
    print("ROBUSTNESS SUMMARY")
    print("=" * 60)
    print(f"Total variants analyzed: {robustness['total_variants']}")
    print(f"Bistable (functional): {robustness['pct_bistable']:.1f}%")
    print(f"Failed: {robustness['pct_failed']:.1f}%")
    print(f"Robustness score: {robustness['robustness_score']:.3f}")
    print("\nFailure mode distribution:")
    for mode, pct in sorted(robustness['failure_percentages'].items(), 
                            key=lambda x: -x[1]):
        print(f"  {mode}: {pct:.1f}%")
    
    # Generate visualizations
    print("\nGenerating visualizations...")
    
    # Get regions for visualization
    _, regions = generate_toggle_switch_sequence()
    
    # Compute position failure rates
    failure_rates = compute_position_failure_rates(
        results_df,
        position_column='mutation_positions',
        label_column='failure_label'
    )
    
    # Failure heatmap
    heatmap_path = output_dir / 'failure_heatmap.png'
    
    # Create region-mode matrix
    failure_modes = ['bistable', 'loss_of_bistability', 'leaky', 
                    'no_expression', 'oscillatory']
    region_names = [r.name for r in regions]
    
    matrix = np.zeros((len(region_names), len(failure_modes)))
    
    for i, region_name in enumerate(region_names):
        region_data = results_df[results_df['region'] == region_name]
        total = len(region_data)
        if total > 0:
            for j, mode in enumerate(failure_modes):
                count = len(region_data[region_data['failure_label'] == mode])
                matrix[i, j] = count / total
    
    plot_failure_heatmap(
        matrix,
        heatmap_path,
        row_labels=region_names,
        col_labels=[m.replace('_', '\n') for m in failure_modes],
        title='Failure Modes by Genomic Region',
        xlabel='Failure Mode',
        ylabel='Genomic Region',
        annot=True
    )
    
    print(f"Saved failure heatmap to: {heatmap_path}")
    
    # Pie chart
    pie_path = output_dir / 'failure_distribution.png'
    plot_failure_mode_pie(
        robustness['failure_distribution'],
        pie_path,
        title='Failure Mode Distribution'
    )
    print(f"Saved failure distribution to: {pie_path}")
    
    print("\n" + "=" * 60)
    print("SIMULATION COMPLETE")
    print("=" * 60)
    print(f"\nOutput files:")
    print(f"  - {results_path}")
    print(f"  - {robustness_path}")
    print(f"  - {heatmap_path}")
    print(f"  - {pie_path}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
