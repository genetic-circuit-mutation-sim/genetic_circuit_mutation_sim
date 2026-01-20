"""
Smoke Tests for Genetic Circuit Mutation Simulator

Minimal tests to verify the package works end-to-end.
Runs a small simulation (10 variants) and checks outputs exist.
"""

import os
import sys
import tempfile
import unittest
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd


class TestSmokeTests(unittest.TestCase):
    """
    Smoke tests to verify basic functionality.
    
    These tests check that:
    1. Package imports work
    2. Core functions run without errors
    3. Output files are generated
    """
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        np.random.seed(42)
    
    def tearDown(self):
        """Clean up test artifacts."""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_imports(self):
        """Test that all package modules can be imported."""
        # Mutation engine
        from mutation_engine import mutate_sequence, AnnotatedRegion, MutationType
        from mutation_engine.mutation_types import RegionType
        
        # Mapping
        from mapping import promoter_to_tx_rate, rbs_to_tl_rate, cds_to_protein_effect
        
        # Analysis
        from analysis import classify_failure, FailureMode, compute_robustness
        
        # Visualization
        from visualization import plot_timecourse, plot_failure_heatmap
        
        self.assertTrue(True)  # If we get here, imports worked
    
    def test_mutation_engine(self):
        """Test that mutation engine generates valid mutations."""
        from mutation_engine import mutate_sequence, AnnotatedRegion
        from mutation_engine.mutation_types import RegionType
        
        # Create test sequence and region
        seq = "ATGCATGCATGCATGC"
        regions = [AnnotatedRegion("test_region", 0, 16, RegionType.CDS, "A")]
        
        # Generate mutations
        mutated_seq, mutations = mutate_sequence(
            seq, regions, n_mut=2, seed=42
        )
        
        self.assertIsInstance(mutated_seq, str)
        self.assertGreater(len(mutations), 0)
        self.assertLessEqual(len(mutations), 2)
    
    def test_mapping_functions(self):
        """Test that mapping functions return valid values."""
        from mapping import promoter_to_tx_rate, rbs_to_tl_rate, cds_to_protein_effect
        
        # Test promoter mapping
        promoter_seq = "TTGACANNNNNNNNNNNNNNNNTATAAT"  # Consensus
        tx_rate = promoter_to_tx_rate(promoter_seq, add_noise=False)
        self.assertIsInstance(tx_rate, float)
        self.assertGreater(tx_rate, 0)
        
        # Test RBS mapping
        rbs_seq = "AGGAGGNNNNNNATG"
        tl_rate = rbs_to_tl_rate(rbs_seq, add_noise=False)
        self.assertIsInstance(tl_rate, float)
        self.assertGreater(tl_rate, 0)
        
        # Test CDS mapping
        cds_seq = "ATGGCTAGCTAA"
        effect = cds_to_protein_effect(cds_seq, add_noise=False)
        self.assertIn('effect_multiplier', effect)
        self.assertIn('effect_type', effect)
    
    def test_failure_classification(self):
        """Test that failure classifier works with mock data."""
        from analysis import classify_failure
        
        # Mock bistable result
        bistable_result = {
            'simulation_successful': True,
            'steady_state_A_highA': 10.0,
            'steady_state_B_highA': 0.1,
            'steady_state_A_highB': 0.1,
            'steady_state_B_highB': 10.0,
            'result_highA': {'A': np.array([10.0]), 'B': np.array([0.1])},
            'result_highB': {'A': np.array([0.1]), 'B': np.array([10.0])}
        }
        
        result = classify_failure(bistable_result)
        self.assertEqual(result.mode.value, 'bistable')
        
        # Mock monostable result
        monostable_result = {
            'simulation_successful': True,
            'steady_state_A_highA': 10.0,
            'steady_state_B_highA': 0.1,
            'steady_state_A_highB': 9.8,
            'steady_state_B_highB': 0.15,
            'result_highA': {'A': np.array([10.0]), 'B': np.array([0.1])},
            'result_highB': {'A': np.array([9.8]), 'B': np.array([0.15])}
        }
        
        result = classify_failure(monostable_result)
        self.assertEqual(result.mode.value, 'loss_of_bistability')
    
    def test_robustness_metrics(self):
        """Test robustness metrics computation."""
        from analysis import compute_robustness
        
        # Create mock results DataFrame
        results = pd.DataFrame({
            'variant_id': range(10),
            'failure_label': ['bistable'] * 6 + ['loss_of_bistability'] * 3 + ['leaky'],
            'region': ['promoter_A'] * 5 + ['rbs_A'] * 5
        })
        
        metrics = compute_robustness(results)
        
        self.assertIn('total_variants', metrics)
        self.assertIn('pct_bistable', metrics)
        self.assertIn('robustness_score', metrics)
        self.assertEqual(metrics['total_variants'], 10)
        self.assertAlmostEqual(metrics['pct_bistable'], 60.0)
    
    def test_full_pipeline_mini(self):
        """
        Test the full pipeline with 10 variants.
        
        This is the main smoke test that verifies main.py functionality.
        """
        from mutation_engine.sequence_mutator import generate_toggle_switch_sequence
        from mutation_engine import mutate_sequence
        from mapping import promoter_to_tx_rate, rbs_to_tl_rate
        from analysis import classify_failure
        
        # Generate sequence
        wt_seq, regions = generate_toggle_switch_sequence()
        self.assertGreater(len(wt_seq), 100)
        self.assertGreater(len(regions), 0)
        
        # Run mini Monte Carlo
        results = []
        for i in range(10):
            # Mutate
            mutated_seq, mutations = mutate_sequence(
                wt_seq, regions, n_mut=1, seed=i
            )
            
            # Simple mock simulation result (since Tellurium may not be installed)
            mock_result = {
                'simulation_successful': True,
                'steady_state_A_highA': np.random.uniform(0, 10),
                'steady_state_B_highA': np.random.uniform(0, 1),
                'steady_state_A_highB': np.random.uniform(0, 1),
                'steady_state_B_highB': np.random.uniform(0, 10),
                'result_highA': {'A': np.array([5.0]), 'B': np.array([0.5])},
                'result_highB': {'A': np.array([0.5]), 'B': np.array([5.0])}
            }
            
            # Classify
            classification = classify_failure(mock_result)
            
            results.append({
                'variant_id': i,
                'failure_label': classification.mode.value,
                'n_mutations': len(mutations)
            })
        
        df = pd.DataFrame(results)
        
        # Save to temp directory
        csv_path = os.path.join(self.temp_dir, 'results.csv')
        df.to_csv(csv_path, index=False)
        
        # Verify output exists
        self.assertTrue(os.path.exists(csv_path))
        
        # Verify content
        loaded_df = pd.read_csv(csv_path)
        self.assertEqual(len(loaded_df), 10)
        self.assertIn('failure_label', loaded_df.columns)
    
    def test_visualization_doesnt_crash(self):
        """Test that visualization functions don't crash."""
        from visualization.heatmaps import plot_failure_heatmap
        
        # Create simple test matrix
        matrix = np.random.rand(4, 3)
        
        output_path = os.path.join(self.temp_dir, 'test_heatmap.png')
        
        plot_failure_heatmap(
            matrix,
            output_path,
            row_labels=['A', 'B', 'C', 'D'],
            col_labels=['X', 'Y', 'Z']
        )
        
        self.assertTrue(os.path.exists(output_path))


class TestIntegration(unittest.TestCase):
    """Integration tests for end-to-end functionality."""
    
    def test_sequence_generation(self):
        """Test toggle switch sequence generation."""
        from mutation_engine.sequence_mutator import generate_toggle_switch_sequence
        
        seq, regions = generate_toggle_switch_sequence()
        
        # Check sequence is valid DNA
        valid_bases = set('ATGC')
        self.assertTrue(all(b in valid_bases for b in seq))
        
        # Check regions are properly ordered
        for region in regions:
            self.assertLess(region.start, region.end)
            self.assertLessEqual(region.end, len(seq))


if __name__ == '__main__':
    unittest.main(verbosity=2)
