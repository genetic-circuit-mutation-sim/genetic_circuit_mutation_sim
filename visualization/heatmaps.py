"""
Heatmap Visualization Module

Provides functions for creating heatmaps of failure frequency
across sequence positions and genomic regions.

Public API:
    plot_failure_heatmap(matrix, outpath) -> None
    plot_region_failure_heatmap(results, regions, outpath) -> None
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

from typing import Dict, List, Optional, Any, Union
import pandas as pd
from pathlib import Path


# Color maps for different visualizations
FAILURE_CMAP = 'RdYlGn_r'  # Red = high failure, Green = low failure
REGION_CMAP = 'Blues'


def plot_failure_heatmap(
    matrix: Union[np.ndarray, pd.DataFrame],
    outpath: Union[str, Path],
    row_labels: Optional[List[str]] = None,
    col_labels: Optional[List[str]] = None,
    title: str = "Failure Frequency Heatmap",
    xlabel: str = "Position",
    ylabel: str = "Region",
    figsize: tuple = (14, 6),
    dpi: int = 150,
    cmap: str = FAILURE_CMAP,
    annot: bool = False
) -> None:
    """
    Create a heatmap showing failure frequency by position/region.
    
    Generates a PNG heatmap for visualizing which sequence positions
    or regions are most vulnerable to mutations.
    
    Args:
        matrix: 2D array or DataFrame of failure frequencies (0-1 scale)
                Rows = regions, Cols = positions or failure modes
        outpath: Output path for PNG file
        row_labels: Labels for rows (e.g., region names)
        col_labels: Labels for columns (e.g., positions or failure modes)
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size in inches
        dpi: Resolution
        cmap: Matplotlib colormap name
        annot: Whether to annotate cells with values
    
    Output:
        Saves PNG heatmap to outpath
    """
    # Convert to numpy if DataFrame
    if isinstance(matrix, pd.DataFrame):
        if row_labels is None:
            row_labels = list(matrix.index)
        if col_labels is None:
            col_labels = list(matrix.columns)
        matrix = matrix.values
    
    if matrix.size == 0:
        # Create placeholder for empty data
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center',
                fontsize=14, transform=ax.transAxes)
        ax.set_title(title)
        plt.savefig(outpath, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        return
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create heatmap
    im = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=0, vmax=1)
    
    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, label='Failure Rate', shrink=0.8)
    
    # Set tick labels
    if row_labels is not None and len(row_labels) == matrix.shape[0]:
        ax.set_yticks(np.arange(len(row_labels)))
        ax.set_yticklabels(row_labels)
    
    if col_labels is not None and len(col_labels) == matrix.shape[1]:
        # For many columns, show subset
        if len(col_labels) > 20:
            step = len(col_labels) // 10
            tick_positions = np.arange(0, len(col_labels), step)
            ax.set_xticks(tick_positions)
            ax.set_xticklabels([col_labels[i] for i in tick_positions], rotation=45)
        else:
            ax.set_xticks(np.arange(len(col_labels)))
            ax.set_xticklabels(col_labels, rotation=45, ha='right')
    
    # Annotate if requested (only for small matrices)
    if annot and matrix.size <= 100:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                text_color = 'white' if matrix[i, j] > 0.5 else 'black'
                ax.text(j, i, f'{matrix[i, j]:.2f}', ha='center', va='center',
                       color=text_color, fontsize=8)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def plot_region_failure_heatmap(
    results: pd.DataFrame,
    regions: List[Dict[str, Any]],
    outpath: Union[str, Path],
    label_column: str = 'failure_label',
    title: str = "Failure Modes by Genomic Region",
    figsize: tuple = (12, 6),
    dpi: int = 150
) -> None:
    """
    Create a heatmap of failure mode distribution across genomic regions.
    
    Args:
        results: DataFrame with simulation results
        regions: List of region dicts with 'name' key
        outpath: Output path for PNG
        label_column: Column containing failure labels
        title: Plot title
        figsize: Figure size
        dpi: Resolution
    
    Output:
        Saves heatmap showing which regions cause which failure modes
    """
    failure_modes = [
        'bistable', 'loss_of_bistability', 'leaky',
        'no_expression', 'oscillatory', 'simulation_failed'
    ]
    
    # Get region names
    region_names = []
    for r in regions:
        if isinstance(r, dict):
            region_names.append(r.get('name', str(r)))
        else:
            region_names.append(getattr(r, 'name', str(r)))
    
    # Build failure matrix
    n_regions = len(region_names)
    n_modes = len(failure_modes)
    matrix = np.zeros((n_regions, n_modes))
    
    if 'region' in results.columns and label_column in results.columns:
        for i, region_name in enumerate(region_names):
            region_data = results[results['region'] == region_name]
            total = len(region_data)
            
            if total > 0:
                for j, mode in enumerate(failure_modes):
                    count = len(region_data[region_data[label_column] == mode])
                    matrix[i, j] = count / total
    
    # Plot heatmap
    plot_failure_heatmap(
        matrix,
        outpath,
        row_labels=region_names,
        col_labels=[m.replace('_', '\n') for m in failure_modes],
        title=title,
        xlabel="Failure Mode",
        ylabel="Genomic Region",
        figsize=figsize,
        dpi=dpi,
        annot=True
    )


def plot_position_failure_profile(
    failure_rates: np.ndarray,
    regions: List[Dict[str, Any]],
    outpath: Union[str, Path],
    title: str = "Position-wise Failure Rate",
    figsize: tuple = (14, 4),
    dpi: int = 150
) -> None:
    """
    Plot failure rate as a function of sequence position.
    
    Shows which positions are most sensitive to mutations,
    with region boundaries marked.
    
    Args:
        failure_rates: Array of failure rates per position
        regions: List of region dicts with 'name', 'start', 'end'
        outpath: Output path
        title: Plot title
        figsize: Figure size
        dpi: Resolution
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    positions = np.arange(len(failure_rates))
    
    # Plot failure rate profile
    ax.fill_between(positions, 0, failure_rates, alpha=0.5, color='#E63946')
    ax.plot(positions, failure_rates, color='#E63946', linewidth=1)
    
    # Mark region boundaries
    colors = plt.cm.Set3(np.linspace(0, 1, len(regions)))
    
    for i, region in enumerate(regions):
        if isinstance(region, dict):
            start = region.get('start', 0)
            end = region.get('end', 0) 
            name = region.get('name', f'Region {i}')
        else:
            start = getattr(region, 'start', 0)
            end = getattr(region, 'end', 0)
            name = getattr(region, 'name', f'Region {i}')
        
        # Shade region
        ax.axvspan(start, end, alpha=0.2, color=colors[i], label=name)
        
        # Add region label
        mid = (start + end) / 2
        ax.text(mid, ax.get_ylim()[1] * 0.95, name, 
               ha='center', va='top', fontsize=8, rotation=45)
    
    ax.set_xlabel('Sequence Position')
    ax.set_ylabel('Failure Rate')
    ax.set_title(title)
    ax.set_xlim(0, len(failure_rates))
    ax.set_ylim(0, 1)
    
    # Add threshold line
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='50% threshold')
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def plot_failure_mode_pie(
    failure_counts: Dict[str, int],
    outpath: Union[str, Path],
    title: str = "Failure Mode Distribution",
    figsize: tuple = (8, 8),
    dpi: int = 150
) -> None:
    """
    Create a pie chart of failure mode distribution.
    
    Args:
        failure_counts: Dict mapping failure mode to count
        outpath: Output path
        title: Plot title
        figsize: Figure size
        dpi: Resolution
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    labels = list(failure_counts.keys())
    sizes = list(failure_counts.values())
    
    # Colors for each failure mode
    color_map = {
        'bistable': '#2E7D32',           # Green - success
        'loss_of_bistability': '#FFA726', # Orange
        'leaky': '#EF5350',               # Red
        'no_expression': '#78909C',       # Gray
        'oscillatory': '#AB47BC',         # Purple
        'simulation_failed': '#424242'     # Dark gray
    }
    
    colors = [color_map.get(label, '#757575') for label in labels]
    
    # Create pie chart
    wedges, texts, autotexts = ax.pie(
        sizes, 
        labels=labels, 
        colors=colors,
        autopct='%1.1f%%',
        startangle=90,
        pctdistance=0.75
    )
    
    # Style
    for text in texts:
        text.set_fontsize(10)
    for autotext in autotexts:
        autotext.set_fontsize(9)
        autotext.set_color('white')
        autotext.set_weight('bold')
    
    ax.set_title(title, fontsize=14)
    
    # Add total count
    total = sum(sizes)
    ax.text(0, 0, f'N={total}', ha='center', va='center', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def create_summary_figure(
    results: pd.DataFrame,
    regions: List[Dict[str, Any]],
    failure_rates: np.ndarray,
    outpath: Union[str, Path],
    dpi: int = 150
) -> None:
    """
    Create a comprehensive summary figure with multiple panels.
    
    Includes:
    - Pie chart of failure modes
    - Position-wise failure profile
    - Region failure heatmap
    
    Args:
        results: DataFrame with simulation results
        regions: List of region dicts
        failure_rates: Position-wise failure rates
        outpath: Output path
        dpi: Resolution
    """
    fig = plt.figure(figsize=(16, 10))
    
    # Create grid
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.25)
    
    # Panel 1: Pie chart (top left)
    ax1 = fig.add_subplot(gs[0, 0])
    if 'failure_label' in results.columns:
        failure_counts = results['failure_label'].value_counts().to_dict()
        labels = list(failure_counts.keys())
        sizes = list(failure_counts.values())
        
        color_map = {
            'bistable': '#2E7D32',
            'loss_of_bistability': '#FFA726',
            'leaky': '#EF5350',
            'no_expression': '#78909C',
            'oscillatory': '#AB47BC',
            'simulation_failed': '#424242'
        }
        colors = [color_map.get(label, '#757575') for label in labels]
        
        ax1.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
        ax1.set_title('Failure Mode Distribution')
    
    # Panel 2: Position profile (top right)
    ax2 = fig.add_subplot(gs[0, 1])
    if len(failure_rates) > 0:
        positions = np.arange(len(failure_rates))
        ax2.fill_between(positions, 0, failure_rates, alpha=0.5, color='#E63946')
        ax2.plot(positions, failure_rates, color='#E63946', linewidth=1)
        ax2.set_xlabel('Sequence Position')
        ax2.set_ylabel('Failure Rate')
        ax2.set_title('Position-wise Failure Rate')
        ax2.set_ylim(0, 1)
    
    # Panel 3: Region heatmap (bottom, spanning both columns)
    ax3 = fig.add_subplot(gs[1, :])
    
    failure_modes = ['bistable', 'loss_of_bistability', 'leaky',
                    'no_expression', 'oscillatory']
    
    region_names = [r.get('name', str(r)) if isinstance(r, dict) 
                   else getattr(r, 'name', str(r)) for r in regions]
    
    matrix = np.zeros((len(region_names), len(failure_modes)))
    
    if 'region' in results.columns and 'failure_label' in results.columns:
        for i, region_name in enumerate(region_names):
            region_data = results[results['region'] == region_name]
            total = len(region_data)
            if total > 0:
                for j, mode in enumerate(failure_modes):
                    count = len(region_data[region_data['failure_label'] == mode])
                    matrix[i, j] = count / total
    
    im = ax3.imshow(matrix, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=1)
    ax3.set_yticks(np.arange(len(region_names)))
    ax3.set_yticklabels(region_names)
    ax3.set_xticks(np.arange(len(failure_modes)))
    ax3.set_xticklabels([m.replace('_', '\n') for m in failure_modes])
    ax3.set_title('Failure Modes by Genomic Region')
    
    fig.colorbar(im, ax=ax3, label='Rate', shrink=0.6)
    
    fig.suptitle('Mutation Robustness Analysis Summary', fontsize=16, y=1.02)
    
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
