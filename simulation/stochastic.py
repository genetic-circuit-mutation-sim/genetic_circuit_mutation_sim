"""
Stochastic Simulation Module

Provides Gillespie algorithm (SSA) simulation for the toggle switch circuit.
This captures inherent stochasticity in gene expression that can affect
circuit behavior, especially near bifurcation points.

Public API:
    simulate_gillespie(model_sbml, n_runs) -> dict
"""

import numpy as np
from typing import Dict, List, Optional, Any
import os

# Try to import simulation backends
try:
    import tellurium as te
    TELLURIUM_AVAILABLE = True
except ImportError:
    TELLURIUM_AVAILABLE = False

try:
    import roadrunner
    ROADRUNNER_AVAILABLE = True
except ImportError:
    ROADRUNNER_AVAILABLE = False


# Default stochastic simulation parameters
DEFAULT_SIM_TIME = 500.0  # Shorter for stochastic (computationally expensive)
DEFAULT_NUM_RUNS = 100
DEFAULT_NUM_POINTS = 501


def simulate_gillespie(
    model_sbml: str,
    n_runs: int = DEFAULT_NUM_RUNS,
    parameters: Optional[Dict[str, float]] = None,
    sim_time: float = DEFAULT_SIM_TIME,
    num_points: int = DEFAULT_NUM_POINTS,
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Run stochastic Gillespie (SSA) simulations of the toggle switch.
    
    Performs multiple runs to capture the distribution of outcomes
    due to stochastic gene expression.
    
    Args:
        model_sbml: SBML model string or Antimony model (will be converted)
        n_runs: Number of stochastic runs (default: 100)
        parameters: Optional parameter modifications
        sim_time: Simulation time per run (default: 500)
        num_points: Number of time points to record (default: 501)
        seed: Random seed for reproducibility
    
    Returns:
        Dict containing:
            - 'time': Common time points array
            - 'runs': List of individual run results
            - 'mean_A': Mean A trajectory across runs
            - 'mean_B': Mean B trajectory across runs
            - 'std_A': Std deviation of A across runs
            - 'std_B': Std deviation of B across runs
            - 'final_A_distribution': Final A values across runs
            - 'final_B_distribution': Final B values across runs
            - 'n_runs': Number of runs completed
            - 'simulation_successful': Boolean success flag
            - 'error_message': Error message if failed
    
    Biological Notes:
        - Stochastic effects are most important at low molecule counts
        - Bistability can be "noisy" - switching between states
        - Multiple runs needed to characterize behavior distribution
    """
    result = {
        'time': None,
        'runs': [],
        'mean_A': None,
        'mean_B': None,
        'std_A': None,
        'std_B': None,
        'final_A_distribution': [],
        'final_B_distribution': [],
        'n_runs': 0,
        'simulation_successful': False,
        'error_message': None
    }
    
    if seed is not None:
        np.random.seed(seed)
    
    # Check for simulation backend
    if not TELLURIUM_AVAILABLE and not ROADRUNNER_AVAILABLE:
        result['error_message'] = "Neither Tellurium nor RoadRunner available"
        return result
    
    try:
        # Load model
        model = _load_stochastic_model(model_sbml)
        
        if model is None:
            result['error_message'] = "Failed to load model for stochastic simulation"
            return result
        
        # Apply parameters
        if parameters:
            _apply_parameters_stochastic(model, parameters)
        
        # Storage for ensemble results
        all_A = []
        all_B = []
        time_points = None
        
        # Run ensemble
        for run_idx in range(n_runs):
            run_result = _run_single_stochastic(
                model, 
                sim_time, 
                num_points,
                run_seed=seed + run_idx if seed else None
            )
            
            if run_result['successful']:
                result['runs'].append(run_result)
                
                if run_result['A'] is not None:
                    all_A.append(run_result['A'])
                    result['final_A_distribution'].append(run_result['A'][-1])
                
                if run_result['B'] is not None:
                    all_B.append(run_result['B'])
                    result['final_B_distribution'].append(run_result['B'][-1])
                
                if time_points is None:
                    time_points = run_result['time']
                
                result['n_runs'] += 1
            
            # Reset model between runs
            model.reset()
            if parameters:
                _apply_parameters_stochastic(model, parameters)
        
        # Calculate ensemble statistics
        if all_A:
            A_array = np.array(all_A)
            result['mean_A'] = np.mean(A_array, axis=0)
            result['std_A'] = np.std(A_array, axis=0)
        
        if all_B:
            B_array = np.array(all_B)
            result['mean_B'] = np.mean(B_array, axis=0)
            result['std_B'] = np.std(B_array, axis=0)
        
        result['time'] = time_points
        result['simulation_successful'] = result['n_runs'] > 0
        
    except Exception as e:
        result['error_message'] = str(e)
        result['simulation_successful'] = False
    
    return result


def _load_stochastic_model(model_source: str):
    """Load a model for stochastic simulation."""
    if TELLURIUM_AVAILABLE:
        try:
            # Check if file path
            if os.path.exists(model_source):
                with open(model_source, 'r') as f:
                    model_str = f.read()
                model = te.loada(model_str)
            else:
                model = te.loada(model_source)
            
            # Configure for stochastic simulation
            model.integrator = 'gillespie'
            return model
            
        except Exception as e:
            print(f"Warning: Failed to load for stochastic simulation: {e}")
    
    return None


def _apply_parameters_stochastic(model, parameters: Dict[str, float]):
    """Apply parameter modifications for stochastic model."""
    for param, value in parameters.items():
        try:
            # Clamp to valid ranges
            if param.startswith('tx_') or param.startswith('tl_'):
                value = np.clip(value, 1e-4, 1e2)
            elif param == 'Kd':
                value = np.clip(value, 1e-3, 1e3)
            elif param == 'n':
                value = np.clip(value, 1.0, 6.0)
            
            model[param] = value
        except Exception as e:
            print(f"Warning: Could not set parameter {param}: {e}")


def _run_single_stochastic(
    model,
    sim_time: float,
    num_points: int,
    run_seed: Optional[int] = None
) -> Dict[str, Any]:
    """Run a single stochastic simulation."""
    result = {
        'time': None,
        'A': None,
        'B': None,
        'mRNA_A': None,
        'mRNA_B': None,
        'successful': False
    }
    
    try:
        # Set seed if provided
        if run_seed is not None:
            model.seed = run_seed
        
        # Set random initial condition (mimic cellular variability)
        model['A'] = np.random.uniform(0.1, 5.0)
        model['B'] = np.random.uniform(0.1, 5.0)
        model['mRNA_A'] = np.random.uniform(0.0, 1.0)
        model['mRNA_B'] = np.random.uniform(0.0, 1.0)
        
        # Run simulation
        sim_result = model.simulate(0, sim_time, num_points)
        
        result['time'] = sim_result[:, 0]
        
        # Extract species (handle different column naming)
        col_names = sim_result.colnames if hasattr(sim_result, 'colnames') else []
        
        for species in ['A', 'B', 'mRNA_A', 'mRNA_B']:
            try:
                for col_name in [f'[{species}]', species]:
                    if col_name in col_names:
                        idx = list(col_names).index(col_name)
                        result[species] = sim_result[:, idx]
                        break
                
                # Fallback by index
                if result[species] is None:
                    species_map = {'A': 3, 'B': 4, 'mRNA_A': 1, 'mRNA_B': 2}
                    if species in species_map and sim_result.shape[1] > species_map[species]:
                        result[species] = sim_result[:, species_map[species]]
            except:
                pass
        
        result['successful'] = True
        
    except Exception as e:
        print(f"Stochastic run error: {e}")
    
    return result


def analyze_bimodality(final_distribution: List[float]) -> Dict[str, float]:
    """
    Analyze a distribution of final values for bimodality.
    
    Useful for detecting bistability in stochastic simulations.
    
    Args:
        final_distribution: List of final concentration values
    
    Returns:
        Dict with bimodality statistics
    """
    if len(final_distribution) < 10:
        return {'bimodal': False, 'n_samples': len(final_distribution)}
    
    values = np.array(final_distribution)
    
    # Calculate basic statistics
    mean = np.mean(values)
    std = np.std(values)
    
    # Check coefficient of variation
    cv = std / mean if mean > 0 else 0
    
    # A high CV (>0.5) suggests bimodality or high noise
    bimodal = cv > 0.5
    
    # Try to identify two modes using simple threshold
    median = np.median(values)
    low_mode = np.mean(values[values < median])
    high_mode = np.mean(values[values >= median])
    
    return {
        'bimodal': bimodal,
        'mean': mean,
        'std': std,
        'cv': cv,
        'low_mode': low_mode,
        'high_mode': high_mode,
        'n_samples': len(final_distribution)
    }
