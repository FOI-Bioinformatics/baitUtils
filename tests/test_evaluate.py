#!/usr/bin/env python3

"""
test_evaluate.py

Comprehensive test suite for the evaluate command and its associated modules.
Tests coverage analysis, gap analysis, and visualization functionality.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
import numpy as np
import pandas as pd

from baitUtils.evaluate import main, validate_inputs, generate_recommendations, print_summary
from baitUtils.coverage_stats import CoverageAnalyzer
from baitUtils.gap_analysis import GapAnalyzer
from baitUtils.coverage_viz import CoverageVisualizer


class TestEvaluateCommand(unittest.TestCase):
    """Test suite for the main evaluate command."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
        
        # Create mock files
        self.mock_oligos_file = self.temp_path / "oligos.fasta"
        self.mock_reference_file = self.temp_path / "reference.fasta"
        self.mock_output_dir = self.temp_path / "output"
        
        # Create mock FASTA files
        with open(self.mock_oligos_file, 'w') as f:
            f.write(">oligo1\nATGCGCGCGCGCATGCGCGCGCGC\n")
            f.write(">oligo2\nGCGCATGCATGCATGCGCGCGCGC\n")
        
        with open(self.mock_reference_file, 'w') as f:
            f.write(">chr1\nATGCGCGCGCGCATGCGCGCGCGCGCGCATGCATGCATGCGCGCGCGCGCGCGCGCATGC\n")
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_validate_inputs_valid(self):
        """Test input validation with valid inputs."""
        args = MagicMock()
        args.input = str(self.mock_oligos_file)
        args.reference = str(self.mock_reference_file)
        args.min_identity = 90.0
        args.min_coverage = 1.0
        args.target_coverage = 10.0
        
        # Should not raise exception
        with patch('subprocess.run'):
            validate_inputs(args)
    
    def test_validate_inputs_missing_files(self):
        """Test input validation with missing files."""
        args = MagicMock()
        args.input = str(self.temp_path / "nonexistent.fasta")
        args.reference = str(self.mock_reference_file)
        args.min_identity = 90.0
        args.min_coverage = 1.0
        args.target_coverage = 10.0
        
        with self.assertRaises(SystemExit):
            validate_inputs(args)
    
    def test_validate_inputs_invalid_parameters(self):
        """Test input validation with invalid parameters."""
        args = MagicMock()
        args.input = str(self.mock_oligos_file)
        args.reference = str(self.mock_reference_file)
        args.min_identity = 150.0  # Invalid
        args.min_coverage = 1.0
        args.target_coverage = 10.0
        
        with self.assertRaises(SystemExit):
            validate_inputs(args)
    
    def test_generate_recommendations(self):
        """Test recommendation generation."""
        coverage_stats = {
            'coverage_breadth': 75.0,
            'mean_depth': 8.5,
            'depth_std': 12.0,
            'mapping_efficiency': 85.0
        }
        
        gap_analysis = {
            'total_gaps': 150,
            'max_gap_size': 15000
        }
        
        recommendations = generate_recommendations(coverage_stats, gap_analysis)
        
        self.assertIsInstance(recommendations, list)
        self.assertTrue(len(recommendations) > 0)
        self.assertIn('coverage breadth', recommendations[0].lower())
    
    def test_print_summary(self):
        """Test summary printing."""
        coverage_stats = {
            'coverage_breadth': 85.5,
            'mean_depth': 12.3,
            'mapping_efficiency': 90.0
        }
        
        gap_analysis = {
            'total_gaps': 25,
            'max_gap_size': 5000
        }
        
        with patch('builtins.print') as mock_print:
            print_summary(coverage_stats, gap_analysis)
            mock_print.assert_called()


class TestCoverageAnalyzer(unittest.TestCase):
    """Test suite for the CoverageAnalyzer class."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
        
        # Create mock PSL file
        self.psl_file = self.temp_path / "test.psl"
        with open(self.psl_file, 'w') as f:
            f.write("psLayout version 3\n")
            f.write("\n")
            f.write("match\tmis-\trep.\tN's\tQ gap\tQ gap\tT gap\tT gap\tstrand\tQ\tQ\tQ\tQ\tT\tT\tT\tT\tblock\tblockSizes\tqStarts\ttStarts\n")
            f.write("\tmatch\tmatch\tcount\tcount\tbases\tcount\tbases\t\tname\tsize\tstart\tend\tname\tsize\tstart\tend\tcount\n")
            f.write("-----------------------------------------------------------------------------------------------\n")
            f.write("120\t0\t0\t0\t0\t0\t0\t0\t+\toligo1\t120\t0\t120\tchr1\t1000\t0\t120\t1\t120,\t0,\t0,\n")
            f.write("115\t5\t0\t0\t0\t0\t0\t0\t+\toligo2\t120\t0\t120\tchr1\t1000\t200\t320\t1\t120,\t0,\t200,\n")
        
        # Create mock reference file
        self.ref_file = self.temp_path / "reference.fasta"
        with open(self.ref_file, 'w') as f:
            f.write(">chr1\n")
            f.write("A" * 1000 + "\n")
        
        self.analyzer = CoverageAnalyzer(
            psl_file=self.psl_file,
            reference_file=self.ref_file,
            min_coverage=1.0,
            target_coverage=10.0
        )
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_load_reference_sequences(self):
        """Test reference sequence loading."""
        self.analyzer._load_reference_sequences()
        
        self.assertIn('chr1', self.analyzer.reference_sequences)
        self.assertEqual(self.analyzer.reference_sequences['chr1']['length'], 1000)
    
    def test_parse_psl_file(self):
        """Test PSL file parsing."""
        self.analyzer._parse_psl_file()
        
        self.assertEqual(len(self.analyzer.mappings), 2)
        self.assertEqual(self.analyzer.mappings[0]['query_name'], 'oligo1')
        self.assertEqual(self.analyzer.mappings[0]['target_start'], 0)
        self.assertEqual(self.analyzer.mappings[0]['target_end'], 120)
    
    def test_compute_coverage_arrays(self):
        """Test coverage array computation."""
        self.analyzer._load_reference_sequences()
        self.analyzer._parse_psl_file()
        self.analyzer._compute_coverage_arrays()
        
        self.assertIn('chr1', self.analyzer.coverage_arrays)
        cov_array = self.analyzer.coverage_arrays['chr1']
        
        # Check that covered regions have coverage > 0
        self.assertTrue(np.all(cov_array[0:120] >= 1))
        self.assertTrue(np.all(cov_array[200:320] >= 1))
    
    def test_analyze_complete(self):
        """Test complete analysis workflow."""
        stats = self.analyzer.analyze()
        
        self.assertIn('total_mappings', stats)
        self.assertIn('coverage_breadth', stats)
        self.assertIn('mean_depth', stats)
        self.assertIn('per_reference', stats)
        
        # Check basic statistics
        self.assertEqual(stats['total_mappings'], 2)
        self.assertEqual(stats['mapped_oligos'], 2)
        self.assertTrue(stats['coverage_breadth'] > 0)
    
    def test_export_coverage_data(self):
        """Test coverage data export."""
        self.analyzer.analyze()
        
        output_dir = self.temp_path / "export_test"
        output_dir.mkdir()
        
        coverage_df, gap_regions = self.analyzer.export_coverage_data(output_dir)
        
        self.assertIsInstance(coverage_df, pd.DataFrame)
        self.assertTrue((output_dir / "data" / "coverage_per_position.csv").exists())


class TestGapAnalyzer(unittest.TestCase):
    """Test suite for the GapAnalyzer class."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
        
        # Create mock reference file
        self.ref_file = self.temp_path / "reference.fasta"
        with open(self.ref_file, 'w') as f:
            f.write(">chr1\n")
            f.write("ATGCGCGCGCGC" + "N" * 20 + "GCGCATGCATGC" + "A" * 100 + "\n")
        
        # Mock coverage data
        self.coverage_data = {
            'reference_length': 1000,
            'per_reference': {
                'chr1': {
                    'length': 132,
                    'coverage_breadth': 60.0,
                    'gaps': 3
                }
            }
        }
        
        self.gap_analyzer = GapAnalyzer(
            coverage_data=self.coverage_data,
            reference_file=self.ref_file,
            min_gap_size=50
        )
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_load_reference_sequences(self):
        """Test reference sequence loading."""
        self.gap_analyzer._load_reference_sequences()
        
        self.assertIn('chr1', self.gap_analyzer.reference_sequences)
        self.assertEqual(len(self.gap_analyzer.reference_sequences['chr1']), 132)
    
    def test_calculate_sequence_features(self):
        """Test sequence feature calculation."""
        sequence = "ATGCGCGCGCGCNNNNNNGCGCATGCATGC"
        features = self.gap_analyzer._calculate_sequence_features(sequence)
        
        self.assertIn('gc_content', features)
        self.assertIn('complexity', features)
        self.assertIn('n_content', features)
        self.assertIn('max_homopolymer', features)
        
        # Check N content calculation
        expected_n_content = (6 / 30) * 100  # 6 N's in 30 bases
        self.assertAlmostEqual(features['n_content'], expected_n_content, places=1)
    
    def test_calculate_shannon_entropy(self):
        """Test Shannon entropy calculation."""
        # Equal distribution of 4 nucleotides
        sequence = "ATGC" * 10
        entropy = self.gap_analyzer._calculate_shannon_entropy(sequence)
        self.assertAlmostEqual(entropy, 2.0, places=1)  # log2(4) = 2
        
        # Single nucleotide (no entropy)
        sequence = "AAAA"
        entropy = self.gap_analyzer._calculate_shannon_entropy(sequence)
        self.assertEqual(entropy, 0.0)
    
    def test_find_max_homopolymer(self):
        """Test homopolymer detection."""
        sequence = "ATGCCCCCGAT"
        max_run = self.gap_analyzer._find_max_homopolymer(sequence)
        self.assertEqual(max_run, 5)  # 5 C's in a row
        
        sequence = "ATGCATGC"
        max_run = self.gap_analyzer._find_max_homopolymer(sequence)
        self.assertEqual(max_run, 1)  # No homopolymers
    
    def test_analyze_complete(self):
        """Test complete gap analysis workflow."""
        results = self.gap_analyzer.analyze()
        
        self.assertIn('total_gaps', results)
        self.assertIn('size_distribution', results)
        self.assertIn('suggestions', results)
        
        # Check that suggestions are generated
        self.assertIsInstance(results['suggestions'], list)
        self.assertTrue(len(results['suggestions']) > 0)


class TestCoverageVisualizer(unittest.TestCase):
    """Test suite for the CoverageVisualizer class."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
        
        # Mock coverage statistics
        self.coverage_stats = {
            'coverage_breadth': 85.5,
            'mean_depth': 12.3,
            'depth_std': 5.2,
            'coverage_cv': 0.42,
            'uniformity_score': 0.58,
            'depth_distribution': {1: 85.5, 5: 70.2, 10: 55.8, 20: 30.1},
            'depth_histogram': {
                'counts': [100, 150, 200, 180, 120, 80, 50, 30],
                'bin_edges': [0, 5, 10, 15, 20, 25, 30, 35, 40]
            },
            'per_reference': {
                'chr1': {'coverage_breadth': 88.0, 'mean_depth': 15.2, 'max_depth': 45, 'gaps': 5, 'length': 10000},
                'chr2': {'coverage_breadth': 82.0, 'mean_depth': 9.8, 'max_depth': 32, 'gaps': 8, 'length': 8000}
            },
            'parameters': {'target_coverage': 10}
        }
        
        # Mock gap analysis
        self.gap_analysis = {
            'total_gaps': 25,
            'mean_gap_size': 250.5,
            'max_gap_size': 5000,
            'gap_percentage': 12.5,
            'size_distribution': {'100-500': 15, '500-1000': 6, '1000-5000': 3, '≥5000': 1},
            'largest_gaps': [
                {'chromosome': 'chr1', 'start': 1000, 'end': 6000, 'length': 5000},
                {'chromosome': 'chr2', 'start': 2000, 'end': 4500, 'length': 2500}
            ]
        }
        
        self.visualizer = CoverageVisualizer(
            coverage_stats=self.coverage_stats,
            gap_analysis=self.gap_analysis,
            output_dir=self.temp_path,
            plot_format='png'
        )
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    @patch('matplotlib.pyplot.savefig')
    @patch('matplotlib.pyplot.close')
    def test_plot_coverage_overview(self, mock_close, mock_savefig):
        """Test coverage overview plot generation."""
        self.visualizer.plot_coverage_overview()
        
        mock_savefig.assert_called_once()
        mock_close.assert_called_once()
    
    @patch('matplotlib.pyplot.savefig')
    @patch('matplotlib.pyplot.close')
    def test_plot_depth_distribution(self, mock_close, mock_savefig):
        """Test depth distribution plot generation."""
        self.visualizer.plot_depth_distribution()
        
        mock_savefig.assert_called_once()
        mock_close.assert_called_once()
    
    @patch('matplotlib.pyplot.savefig')
    @patch('matplotlib.pyplot.close')
    def test_plot_gap_analysis(self, mock_close, mock_savefig):
        """Test gap analysis plot generation."""
        self.visualizer.plot_gap_analysis()
        
        mock_savefig.assert_called_once()
        mock_close.assert_called_once()
    
    @patch('matplotlib.pyplot.savefig')
    @patch('matplotlib.pyplot.close')
    def test_plot_coverage_heatmap(self, mock_close, mock_savefig):
        """Test coverage heatmap generation."""
        self.visualizer.plot_coverage_heatmap()
        
        mock_savefig.assert_called_once()
        mock_close.assert_called_once()
    
    @patch('matplotlib.pyplot.savefig')
    @patch('matplotlib.pyplot.close')
    def test_plot_per_reference_coverage(self, mock_close, mock_savefig):
        """Test per-reference coverage plots."""
        self.visualizer.plot_per_reference_coverage()
        
        mock_savefig.assert_called_once()
        mock_close.assert_called_once()
    
    def test_plots_directory_creation(self):
        """Test that plots directory is created."""
        plots_dir = self.visualizer.plots_dir
        self.assertTrue(plots_dir.exists())
        self.assertTrue(plots_dir.is_dir())


class TestIntegrationScenarios(unittest.TestCase):
    """Integration tests for realistic usage scenarios."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_high_coverage_scenario(self):
        """Test scenario with good coverage."""
        coverage_stats = {
            'coverage_breadth': 95.0,
            'mean_depth': 15.0,
            'depth_std': 3.0,
            'mapping_efficiency': 92.0
        }
        
        gap_analysis = {
            'total_gaps': 5,
            'max_gap_size': 200
        }
        
        recommendations = generate_recommendations(coverage_stats, gap_analysis)
        
        # Should indicate good coverage
        self.assertTrue(any('good' in rec.lower() for rec in recommendations))
    
    def test_low_coverage_scenario(self):
        """Test scenario with poor coverage."""
        coverage_stats = {
            'coverage_breadth': 45.0,
            'mean_depth': 3.0,
            'depth_std': 8.0,
            'mapping_efficiency': 60.0
        }
        
        gap_analysis = {
            'total_gaps': 200,
            'max_gap_size': 25000
        }
        
        recommendations = generate_recommendations(coverage_stats, gap_analysis)
        
        # Should suggest improvements
        self.assertTrue(len(recommendations) > 2)
        self.assertTrue(any('breadth' in rec.lower() for rec in recommendations))
        self.assertTrue(any('gaps' in rec.lower() for rec in recommendations))


if __name__ == '__main__':
    # Set up test suite
    unittest.main(verbosity=2)