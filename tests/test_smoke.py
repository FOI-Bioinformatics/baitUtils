#!/usr/bin/env python3
"""
Smoke tests for baitUtils core functionality.

These tests verify basic functionality works without external dependencies.
They should always pass and serve as a minimum functional verification.
"""

import unittest
import tempfile
import os
from pathlib import Path

# Core imports - these should never fail
from baitUtils.sequence_analysis import SequenceAnalyzer
from baitUtils.sequence_statistics import SequenceStatsCalculator, SequenceFilter
from baitUtils.plotting_utils import validate_columns


class TestSmokeBasicImports(unittest.TestCase):
    """Test that core modules can be imported."""
    
    def test_sequence_analysis_import(self):
        """Test sequence analysis module imports."""
        from baitUtils.sequence_analysis import SequenceAnalyzer, reverse_complement, clean_sequence
        self.assertTrue(True)  # If we got here, import succeeded
    
    def test_sequence_statistics_import(self):
        """Test sequence statistics module imports."""
        from baitUtils.sequence_statistics import SequenceStatsCalculator, SequenceFilter
        self.assertTrue(True)  # If we got here, import succeeded
    
    def test_plotting_utils_import(self):
        """Test plotting utilities import."""
        from baitUtils.plotting_utils import PlotGenerator, StatisticalPlotter
        self.assertTrue(True)  # If we got here, import succeeded

    def test_main_package_import(self):
        """Test main package imports properly."""
        import baitUtils
        self.assertTrue(hasattr(baitUtils, '__version__'))


class TestSmokeSequenceAnalysis(unittest.TestCase):
    """Test basic sequence analysis functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = SequenceAnalyzer()
    
    def test_gc_content_calculation(self):
        """Test GC content calculation with known sequences."""
        # 100% AT
        gc_at = self.analyzer.calculate_gc_content('AAATTT')
        self.assertEqual(gc_at, 0.0)
        
        # 100% GC
        gc_gc = self.analyzer.calculate_gc_content('GGGCCC')
        self.assertEqual(gc_gc, 100.0)
        
        # 50% GC
        gc_mixed = self.analyzer.calculate_gc_content('ATGC')
        self.assertEqual(gc_mixed, 50.0)
    
    def test_entropy_calculation(self):
        """Test entropy calculation."""
        # Uniform distribution should have high entropy
        entropy_uniform = self.analyzer.calculate_entropy('ATGC')
        self.assertEqual(entropy_uniform, 2.0)  # log2(4) = 2
        
        # Single character should have zero entropy
        entropy_single = self.analyzer.calculate_entropy('AAAA')
        self.assertEqual(entropy_single, 0.0)
        
        # Empty sequence should have zero entropy
        entropy_empty = self.analyzer.calculate_entropy('')
        self.assertEqual(entropy_empty, 0.0)
    
    def test_n_count(self):
        """Test N base counting."""
        n_count = self.analyzer.count_n_bases('ATGCNNN')
        self.assertEqual(n_count, 3)
        
        n_count_none = self.analyzer.count_n_bases('ATGC')
        self.assertEqual(n_count_none, 0)
    
    def test_masked_count(self):
        """Test masked base counting."""
        masked_count = self.analyzer.count_masked_bases('ATGCatgc')
        self.assertEqual(masked_count, 4)
        
        masked_count_none = self.analyzer.count_masked_bases('ATGC')
        self.assertEqual(masked_count_none, 0)

    def test_complexity_calculation(self):
        """Test sequence complexity calculation."""
        # Perfect complexity (all unique 2-mers)
        complexity_high = self.analyzer.calculate_complexity('ATGC', k=2)
        self.assertEqual(complexity_high, 1.0)  # 3 unique 2-mers out of 3 total
        
        # Low complexity (repeated 2-mers) 
        complexity_low = self.analyzer.calculate_complexity('AAAA', k=2)
        expected = 1/3  # 1 unique 2-mer (AA) out of 3 total
        self.assertAlmostEqual(complexity_low, expected, places=5)

    def test_homopolymer_runs(self):
        """Test homopolymer run detection."""
        runs = self.analyzer.count_homopolymer_runs('AAACCCGGGTTT', min_length=3)
        self.assertEqual(len(runs), 4)  # AAA, CCC, GGG, TTT
        
        # Check specific run details
        bases = [run[0] for run in runs]
        lengths = [run[2] for run in runs]
        self.assertEqual(bases, ['A', 'C', 'G', 'T'])
        self.assertEqual(lengths, [3, 3, 3, 3])


class TestSmokeSequenceStatistics(unittest.TestCase):
    """Test basic sequence statistics functionality."""
    
    def test_sequence_filter_initialization(self):
        """Test sequence filter can be initialized."""
        filter_obj = SequenceFilter(length=120, complete=True, no_ns=True)
        self.assertEqual(filter_obj.length, 120)
        self.assertTrue(filter_obj.complete)
        self.assertTrue(filter_obj.no_ns)
    
    def test_sequence_filter_logic(self):
        """Test basic filtering logic."""
        filter_obj = SequenceFilter(length=10, complete=True)
        
        # Mock stats for a sequence that should pass
        good_stats = {
            'n_count': 0,
            'gc_content': 50.0,
            'melting_temperature': 55.0,
            'masked_percentage': 0.0
        }
        
        # Test with correct length sequence
        result = filter_obj.passes_filters('A' * 10, good_stats)
        self.assertTrue(result)
        
        # Test with wrong length sequence  
        result = filter_obj.passes_filters('A' * 5, good_stats)
        self.assertFalse(result)
    
    def test_stats_calculator_initialization(self):
        """Test stats calculator can be initialized."""
        calculator = SequenceStatsCalculator(num_processes=1)
        self.assertEqual(calculator.num_processes, 1)
        self.assertIsInstance(calculator.analyzer, SequenceAnalyzer)


class TestSmokeUtilityFunctions(unittest.TestCase):
    """Test utility functions."""
    
    def test_reverse_complement(self):
        """Test reverse complement function."""
        from baitUtils.sequence_analysis import reverse_complement
        
        result = reverse_complement('ATGC')
        self.assertEqual(result, 'GCAT')
        
        result = reverse_complement('AAATTT')
        self.assertEqual(result, 'AAATTT')
    
    def test_clean_sequence(self):
        """Test sequence cleaning function."""
        from baitUtils.sequence_analysis import clean_sequence
        
        # Remove invalid characters
        result = clean_sequence('ATGCXYZ123')
        self.assertEqual(result, 'ATGC')
        
        # Remove N characters when requested
        result = clean_sequence('ATGCNNN', remove_ns=True)
        self.assertEqual(result, 'ATGC')
        
        # Keep N characters by default
        result = clean_sequence('ATGCNNN', remove_ns=False)
        self.assertEqual(result, 'ATGCNNN')
    
    def test_validate_columns_function(self):
        """Test column validation utility."""
        import pandas as pd
        
        # Create test DataFrame
        df = pd.DataFrame({
            'col1': [1, 2, 3],
            'col2': [4.0, 5.0, 6.0],
            'col3': ['a', 'b', 'c']  # Non-numeric
        })
        
        # Test with valid numeric columns
        result = validate_columns(df, ['col1', 'col2'], ['histogram'])
        self.assertEqual(result, ['col1', 'col2'])
        
        # Test with mixed columns (should filter out non-numeric)
        result = validate_columns(df, ['col1', 'col3'], ['histogram'])
        self.assertEqual(result, ['col1'])


class TestSmokeEndToEnd(unittest.TestCase):
    """End-to-end smoke tests with real data."""
    
    def test_basic_sequence_analysis_workflow(self):
        """Test basic workflow with small dataset."""
        # Create test sequences
        sequences = [
            ('seq1', 'ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC'),
            ('seq2', 'GCTATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'),
            ('seq3', 'CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT')
        ]
        
        # Initialize calculator
        calculator = SequenceStatsCalculator(num_processes=1)
        
        # Process sequences (minimal simulation)
        results = []
        for seq_id, sequence in sequences:
            stats = calculator.analyzer.analyze_sequence(sequence, seq_id)
            results.append(stats)
        
        # Verify we got results for all sequences
        self.assertEqual(len(results), 3)
        
        # Verify basic stats are present
        for result in results:
            self.assertIn('sequence_id', result)
            self.assertIn('length', result)
            self.assertIn('gc_content', result)
            self.assertIsInstance(result['length'], int)
            self.assertIsInstance(result['gc_content'], float)


if __name__ == '__main__':
    unittest.main()