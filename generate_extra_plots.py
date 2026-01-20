
import os
import glob
import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def plot_parameter_distributions(df, out_dir):
    """Plot boxplots of parameters for bistable vs failed variants."""
    params = ['param_tx_A', 'param_tx_B', 'param_tl_A', 'param_tl_B', 'param_Kd', 'param_n']
    
    # Filter for valid data
    df_clean = df.dropna(subset=params + ['failure_label'])
    
    # Normalize failure labels to just 'Bistable' vs 'Failed' for cleaner plotting
    df_clean['Group'] = df_clean['failure_label'].apply(
        lambda x: 'Bistable' if x == 'bistable' else 'Failed'
    )
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    colors = {'Bistable': '#2E7D32', 'Failed': '#EF5350'}
    
    for i, param in enumerate(params):
        ax = axes[i]
        
        # Prepare data for boxplot
        data_to_plot = [
            df_clean[df_clean['Group'] == 'Bistable'][param],
            df_clean[df_clean['Group'] == 'Failed'][param]
        ]
        
        # Create boxplot
        bp = ax.boxplot(data_to_plot, patch_artist=True, labels=['Bistable', 'Failed'])
        
        # Color boxes
        for patch, label in zip(bp['boxes'], ['Bistable', 'Failed']):
            patch.set_facecolor(colors[label])
            patch.set_alpha(0.6)
            
        ax.set_title(param)
        ax.grid(True, linestyle='--', alpha=0.3)
        
        # Use log scale for Kd and rates if range is large
        if param in ['param_Kd', 'param_tx_A', 'param_tx_B']:
            ax.set_yscale('log')
            
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'parameter_distributions.png'), dpi=150)
    plt.close()
    print(f"Saved parameter distributions to {out_dir}")

def plot_state_space(df, out_dir):
    """Plot steady state A vs B."""
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Define colors
    color_map = {
        'bistable': '#2E7D32',           # Green
        'loss_of_bistability': '#FFA726', # Orange
        'leaky': '#EF5350',               # Red
        'no_expression': '#78909C',       # Gray
        'oscillatory': '#AB47BC',         # Purple
        'simulation_failed': '#424242'     # Dark gray
    }
    
    # Plot each category
    for label, color in color_map.items():
        subset = df[df['failure_label'] == label]
        if not subset.empty:
            ax.scatter(
                subset['ss_A_highA'], 
                subset['ss_B_highA'], 
                c=color, 
                label=label,
                alpha=0.6,
                edgecolors='none',
                s=30
            )
            
    ax.set_xlabel('[Repressor A] (nM)')
    ax.set_ylabel('[Repressor B] (nM)')
    ax.set_title('State Space Distribution (High A Initial Condition)')
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.3)
    
    plt.savefig(os.path.join(out_dir, 'state_space_scatter.png'), dpi=150)
    plt.close()
    print(f"Saved state space scatter to {out_dir}")

def plot_comparative_robustness(results_root):
    """Compare robustness across different scenarios."""
    summary_files = glob.glob(os.path.join(results_root, '*', 'robustness_summary.json'))
    
    data = []
    for f in summary_files:
        scenario_name = os.path.basename(os.path.dirname(f))
        try:
            with open(f, 'r') as json_file:
                stats = json.load(json_file)
                data.append({
                    'Scenario': scenario_name.replace('_', ' ').title(),
                    'Bistable (%)': stats.get('pct_bistable', 0),
                    'Robustness Score': stats.get('robustness_score', 0)
                })
        except Exception as e:
            print(f"Error reading {f}: {e}")
            
    if not data:
        return

    df = pd.DataFrame(data)
    df = df.sort_values('Bistable (%)', ascending=False)
    
    # Plotting
    fig, ax1 = plt.subplots(figsize=(12, 6))
    
    # Bar chart for Bistable %
    bars = ax1.bar(df['Scenario'], df['Bistable (%)'], color='#2196F3', alpha=0.7, label='Bistable %')
    ax1.set_ylabel('Percent Bistable Variants (%)', color='#2196F3', fontsize=12)
    ax1.set_ylim(0, 100)
    ax1.set_title('Robustness Comparison Across Scenarios', fontsize=14)
    
    # Add values on top of bars
    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{height:.1f}%', ha='center', va='bottom')
                
    plt.xticks(rotation=15)
    plt.tight_layout()
    
    out_path = os.path.join(results_root, 'comparative_robustness.png')
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved comparative chart to {out_path}")

def plot_parameter_sensitivity_heatmap(df, out_dir):
    """
    Generate a heatmap showing how much mutations in each region affect each parameter.
    X-axis: Parameters
    Y-axis: Regions
    Value: Average fold-change relative to WT (log scale or ratio)
    """
    # Define default params for normalization (approximate from main.py)
    defaults = {
        'param_tx_A': 1.0, 'param_tx_B': 1.0,
        'param_tl_A': 1.0, 'param_tl_B': 1.0,
        'param_Kd': 1.0, 'param_n': 2.0
    }
    
    # Clean data
    df_clean = df.dropna(subset=['region'] + list(defaults.keys()))
    df_clean = df_clean[df_clean['region'] != '']
    
    # Calculate fold change (absolute log ratio) helps visualize magnitude of change up or down
    # fold_change = |log(val / default)|
    # We add a small epsilon to avoid log(0)
    for param, default_val in defaults.items():
        # Handle zeros by setting to very small number
        vals = df_clean[param].replace(0, 1e-6)
        df_clean[f'change_{param}'] = np.abs(np.log10(vals / default_val))
        
    # Group by region and calc mean change
    change_cols = [f'change_{p}' for p in defaults.keys()]
    sensitivity = df_clean.groupby('region')[change_cols].mean()
    
    if sensitivity.empty:
        return

    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot heatmap
    im = ax.imshow(sensitivity, cmap='Reds', aspect='auto')
    
    # Labels
    ax.set_xticks(np.arange(len(defaults)))
    ax.set_xticklabels([p.replace('param_', '') for p in defaults.keys()], rotation=45)
    
    ax.set_yticks(np.arange(len(sensitivity.index)))
    ax.set_yticklabels(sensitivity.index)
    
    ax.set_title('Parameter Sensitivity by Region\n(Mean |log10 fold-change|)')
    
    plt.colorbar(im, label='Magnitude of Change')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'parameter_sensitivity_heatmap.png'), dpi=150)
    plt.close()
    print(f"Saved param sensitivity heatmap to {out_dir}")

def plot_mutation_type_region_heatmap(df, out_dir):
    """
    Generate a heatmap of Failure Rate for (Region x Mutation Type).
    X-axis: Mutation Types
    Y-axis: Regions
    Value: Failure Rate (0-1)
    """
    if 'mutation_type' not in df.columns or 'region' not in df.columns:
        return
        
    df_clean = df[df['region'] != '']
    df_clean = df_clean[df_clean['mutation_type'] != '']
    
    # Calculate failure rates
    # We group by Region and Type, then calc % failed
    # Failure is defined as failure_label != 'bistable'
    df_clean['is_failed'] = df_clean['failure_label'] != 'bistable'
    
    grouped = df_clean.groupby(['region', 'mutation_type'])['is_failed'].mean().unstack()
    
    # Fill NaN with 0 (or could be masked)
    grouped = grouped.fillna(0)
    
    if grouped.empty:
        return
        
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Plot
    im = ax.imshow(grouped, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=1)
    
    ax.set_xticks(np.arange(len(grouped.columns)))
    ax.set_xticklabels(grouped.columns, rotation=45)
    
    ax.set_yticks(np.arange(len(grouped.index)))
    ax.set_yticklabels(grouped.index)
    
    ax.set_title('Failure Rate by Region & Mutation Type')
    
    plt.colorbar(im, label='Failure Rate')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'mutation_type_region_heatmap.png'), dpi=150)
    plt.close()
    print(f"Saved mutation type heatmap to {out_dir}")

def main():
    results_root = 'results'
    
    # Process each subdirectory
    subdirs = [f.path for f in os.scandir(results_root) if f.is_dir()]
    
    for subdir in subdirs:
        print(f"Processing {subdir}...")
        csv_path = os.path.join(subdir, 'results.csv')
        
        if os.path.exists(csv_path):
            try:
                df = pd.read_csv(csv_path)
                plot_parameter_distributions(df, subdir)
                plot_state_space(df, subdir)
                plot_parameter_sensitivity_heatmap(df, subdir)  # New
                plot_mutation_type_region_heatmap(df, subdir)   # New
            except Exception as e:
                print(f"Failed to process {subdir}: {e}")
        else:
            print(f"No results.csv found in {subdir}")
            
    # Generate aggregate comparison
    print("Generating comparative analysis...")
    plot_comparative_robustness(results_root)

if __name__ == "__main__":
    main()
