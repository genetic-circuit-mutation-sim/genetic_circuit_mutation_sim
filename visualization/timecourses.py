"""
Time Course Visualization Module

Provides functions for visualizing toggle switch simulation time series,
including multi-condition comparisons and trajectory overlays.

Public API:
    plot_timecourse(df, outpath) -> None
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for file output

from typing import Dict, List, Optional, Any, Union
import pandas as pd
from pathlib import Path


# Default style settings
STYLE_CONFIG = {
    'figure.figsize': (10, 6),
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
    'lines.linewidth': 2,
    'axes.grid': True,
    'grid.alpha': 0.3
}

# Color palette for species
COLORS = {
    'A': '#E63946',      # Red for protein A
    'B': '#457B9D',      # Blue for protein B
    'mRNA_A': '#F4A261', # Orange for mRNA A
    'mRNA_B': '#2A9D8F', # Teal for mRNA B
}


def plot_timecourse(
    data: Union[pd.DataFrame, Dict[str, Any]],
    outpath: Union[str, Path],
    title: str = "Toggle Switch Dynamics",
    show_mrna: bool = False,
    figsize: tuple = (10, 6),
    dpi: int = 150
) -> None:
    """
    Plot protein (and optionally mRNA) concentration time series.
    
    Generates a publication-quality figure showing toggle switch
    dynamics over time.
    
    Args:
        data: Either:
            - DataFrame with 'time', 'A', 'B', (optional 'mRNA_A', 'mRNA_B') columns
            - Dict from simulate_deterministic() result
        outpath: Output path for saved figure (PNG, PDF, etc.)
        title: Plot title
        show_mrna: Whether to show mRNA concentrations
        figsize: Figure size in inches
        dpi: Resolution in dots per inch
    
    Output:
        Saves figure to specified outpath
    
    Example:
        >>> plot_timecourse(sim_result, 'trajectory.png')
    """
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # Extract data
    if isinstance(data, pd.DataFrame):
        time = data['time'].values if 'time' in data else np.arange(len(data))
        A = data['A'].values if 'A' in data else None
        B = data['B'].values if 'B' in data else None
        mRNA_A = data.get('mRNA_A', pd.Series()).values if show_mrna else None
        mRNA_B = data.get('mRNA_B', pd.Series()).values if show_mrna else None
        two_start = False
    else:
        # Dict from simulation result
        result_highA = data.get('result_highA', {})
        result_highB = data.get('result_highB', {})
        time = data.get('time')
        two_start = result_highB is not None and 'A' in result_highB
        
        if isinstance(result_highA, dict):
            A = result_highA.get('A')
            B = result_highA.get('B')
            mRNA_A = result_highA.get('mRNA_A') if show_mrna else None
            mRNA_B = result_highA.get('mRNA_B') if show_mrna else None
            
            if two_start:
                A_highB = result_highB.get('A')
                B_highB = result_highB.get('B')
        else:
            A = B = mRNA_A = mRNA_B = None
    
    # Create figure
    if show_mrna:
        fig, axes = plt.subplots(2, 1, figsize=(figsize[0], figsize[1] * 1.5), 
                                  sharex=True)
        ax_protein = axes[0]
        ax_mrna = axes[1]
    else:
        fig, ax_protein = plt.subplots(figsize=figsize)
        ax_mrna = None
    
    # Plot proteins
    if A is not None and time is not None:
        ax_protein.plot(time, A, color=COLORS['A'], label='Protein A', 
                       linewidth=2)
    if B is not None and time is not None:
        ax_protein.plot(time, B, color=COLORS['B'], label='Protein B',
                       linewidth=2)
    
    # Plot second initial condition if available
    if two_start and 'A_highB' in dir():
        if A_highB is not None:
            ax_protein.plot(time, A_highB, color=COLORS['A'], linestyle='--',
                           alpha=0.7, label='Protein A (high-B start)')
        if B_highB is not None:
            ax_protein.plot(time, B_highB, color=COLORS['B'], linestyle='--',
                           alpha=0.7, label='Protein B (high-B start)')
    
    ax_protein.set_ylabel('Protein Concentration (AU)')
    ax_protein.set_title(title)
    ax_protein.legend(loc='best')
    ax_protein.set_ylim(bottom=0)
    
    # Plot mRNA if requested
    if show_mrna and ax_mrna is not None:
        if mRNA_A is not None and time is not None:
            ax_mrna.plot(time, mRNA_A, color=COLORS['mRNA_A'], 
                        label='mRNA A', linewidth=2)
        if mRNA_B is not None and time is not None:
            ax_mrna.plot(time, mRNA_B, color=COLORS['mRNA_B'],
                        label='mRNA B', linewidth=2)
        
        ax_mrna.set_xlabel('Time (AU)')
        ax_mrna.set_ylabel('mRNA Concentration (AU)')
        ax_mrna.legend(loc='best')
        ax_mrna.set_ylim(bottom=0)
    else:
        ax_protein.set_xlabel('Time (AU)')
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def plot_two_start_comparison(
    sim_result: Dict[str, Any],
    outpath: Union[str, Path],
    title: str = "Two-Start Bistability Test",
    figsize: tuple = (12, 5),
    dpi: int = 150
) -> None:
    """
    Plot side-by-side comparison of high-A and high-B starting conditions.
    
    This is the standard visualization for bistability testing.
    
    Args:
        sim_result: Result dict from simulate_deterministic()
        outpath: Output path for figure
        title: Plot title
        figsize: Figure size
        dpi: Resolution
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    
    time = sim_result.get('time')
    
    # High-A start
    result_highA = sim_result.get('result_highA', {})
    if isinstance(result_highA, dict):
        A_highA = result_highA.get('A')
        B_highA = result_highA.get('B')
        
        if A_highA is not None and time is not None:
            axes[0].plot(time, A_highA, color=COLORS['A'], label='Protein A', linewidth=2)
            axes[0].plot(time, B_highA, color=COLORS['B'], label='Protein B', linewidth=2)
        
        axes[0].set_title('High-A Start')
        axes[0].set_xlabel('Time (AU)')
        axes[0].set_ylabel('Concentration (AU)')
        axes[0].legend()
    
    # High-B start
    result_highB = sim_result.get('result_highB', {})
    if isinstance(result_highB, dict):
        A_highB = result_highB.get('A')
        B_highB = result_highB.get('B')
        
        if A_highB is not None and time is not None:
            axes[1].plot(time, A_highB, color=COLORS['A'], label='Protein A', linewidth=2)
            axes[1].plot(time, B_highB, color=COLORS['B'], label='Protein B', linewidth=2)
        
        axes[1].set_title('High-B Start')
        axes[1].set_xlabel('Time (AU)')
        axes[1].legend()
    
    fig.suptitle(title, fontsize=14)
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def plot_stochastic_ensemble(
    stoch_result: Dict[str, Any],
    outpath: Union[str, Path],
    title: str = "Stochastic Ensemble",
    show_individual: bool = True,
    max_individual: int = 20,
    figsize: tuple = (10, 6),
    dpi: int = 150
) -> None:
    """
    Plot stochastic simulation ensemble with mean and confidence interval.
    
    Args:
        stoch_result: Result from simulate_gillespie()
        outpath: Output path
        title: Plot title
        show_individual: Whether to show individual trajectories
        max_individual: Maximum number of individual trajectories to show
        figsize: Figure size
        dpi: Resolution
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    time = stoch_result.get('time')
    mean_A = stoch_result.get('mean_A')
    mean_B = stoch_result.get('mean_B')
    std_A = stoch_result.get('std_A')
    std_B = stoch_result.get('std_B')
    runs = stoch_result.get('runs', [])
    
    if time is None or mean_A is None:
        plt.close(fig)
        return
    
    # Plot individual trajectories (faint)
    if show_individual and runs:
        for i, run in enumerate(runs[:max_individual]):
            if run.get('A') is not None:
                ax.plot(time, run['A'], color=COLORS['A'], alpha=0.1, linewidth=0.5)
            if run.get('B') is not None:
                ax.plot(time, run['B'], color=COLORS['B'], alpha=0.1, linewidth=0.5)
    
    # Plot mean ± std
    if mean_A is not None and std_A is not None:
        ax.fill_between(time, mean_A - std_A, mean_A + std_A, 
                        color=COLORS['A'], alpha=0.3)
        ax.plot(time, mean_A, color=COLORS['A'], linewidth=2, 
                label=f'Protein A (n={stoch_result.get("n_runs", 0)})')
    
    if mean_B is not None and std_B is not None:
        ax.fill_between(time, mean_B - std_B, mean_B + std_B,
                        color=COLORS['B'], alpha=0.3)
        ax.plot(time, mean_B, color=COLORS['B'], linewidth=2,
                label=f'Protein B (n={stoch_result.get("n_runs", 0)})')
    
    ax.set_xlabel('Time (AU)')
    ax.set_ylabel('Concentration (AU)')
    ax.set_title(title)
    ax.legend()
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
