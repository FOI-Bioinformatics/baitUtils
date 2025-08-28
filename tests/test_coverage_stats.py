#!/usr/bin/env python3

"""
test_coverage_stats.py

Unit tests for the coverage_stats module.
Tests all statistical calculations and data processing functions.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock
import numpy as np
import pandas as pd

from baitUtils.coverage_stats import CoverageAnalyzer


class TestCoverageStatsCalculations(unittest.TestCase):
    """Test statistical calculation functions."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
        
        # Create comprehensive test PSL file
        self.psl_file = self.temp_path / "comprehensive_test.psl"
        with open(self.psl_file, 'w') as f:
            f.write("psLayout version 3\n")
            f.write("\n")
            f.write("match\tmis-\trep.\tN's\tQ gap\tQ gap\tT gap\tT gap\tstrand\tQ\tQ\tQ\tQ\tT\tT\tT\tT\tblock\tblockSizes\tqStarts\ttStarts\n")
            f.write("\tmatch\tmatch\tcount\tcount\tbases\tcount\tbases\t\tname\tsize\tstart\tend\tname\tsize\tstart\tend\tcount\n")
            f.write("-----------------------------------------------------------------------------------------------\n")
            # High identity mapping
            f.write("120\t0\t0\t0\t0\t0\t0\t0\t+\toligo1\t120\t0\t120\tchr1\t2000\t0\t120\t1\t120,\t0,\t0,\n")
            # Medium identity mapping
            f.write("100\t15\t5\t0\t0\t0\t0\t0\t+\toligo2\t120\t0\t120\tchr1\t2000\t500\t620\t1\t120,\t0,\t500,\n")
            # Low identity mapping (should be filtered out)
            f.write("80\t30\t10\t0\t0\t0\t0\t0\t+\toligo3\t120\t0\t120\tchr1\t2000\t1000\t1120\t1\t120,\t0,\t1000,\n")
            # Another chromosome
            f.write("110\t10\t0\t0\t0\t0\t0\t0\t+\toligo4\t120\t0\t120\tchr2\t1500\t100\t220\t1\t120,\t0,\t100,\n")
            # Overlapping mapping
            f.write("115\t5\t0\t0\t0\t0\t0\t0\t+\toligo5\t120\t0\t120\tchr1\t2000\t50\t170\t1\t120,\t0,\t50,\n")
        
        # Create comprehensive reference file
        self.ref_file = self.temp_path / "reference.fasta"
        with open(self.ref_file, 'w') as f:
            f.write(">chr1\n")
            f.write("A" * 2000 + "\n")
            f.write(">chr2\n")
            f.write("G" * 1500 + "\n")
        
        self.analyzer = CoverageAnalyzer(
            psl_file=self.psl_file,
            reference_file=self.ref_file,
            min_coverage=1.0,
            target_coverage=5.0,
            min_identity=85.0,
            min_length=100
        )
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_psl_filtering(self):
        """Test PSL file parsing and filtering."""
        self.analyzer._parse_psl_file()
        
        # Should have 4 mappings (oligo3 filtered out due to low identity)
        self.assertEqual(len(self.analyzer.mappings), 4)
        
        # Check identity filtering
        identities = [m['identity'] for m in self.analyzer.mappings]
        self.assertTrue(all(identity >= 85.0 for identity in identities))
        
        # Check oligo3 is filtered out
        oligo_names = [m['query_name'] for m in self.analyzer.mappings]
        self.assertNotIn('oligo3', oligo_names)
    
    def test_coverage_array_computation(self):
        """Test coverage depth array computation."""
        self.analyzer._load_reference_sequences()
        self.analyzer._parse_psl_file()
        self.analyzer._compute_coverage_arrays()
        
        # Check chr1 coverage
        chr1_cov = self.analyzer.coverage_arrays['chr1']
        
        # Position 0-50: oligo1 only (depth 1)
        self.assertTrue(np.all(chr1_cov[0:50] == 1))
        
        # Position 50-120: oligo1 + oligo5 overlap (depth 2) 
        self.assertTrue(np.all(chr1_cov[50:120] == 2))
        
        # Position 120-170: only oligo5 (depth 1)
        self.assertTrue(np.all(chr1_cov[120:170] == 1))
        
        # Position 500-620: oligo2 (depth 1)
        self.assertTrue(np.all(chr1_cov[500:620] == 1))
        
        # Uncovered regions should have depth 0
        self.assertTrue(np.all(chr1_cov[170:500] == 0))
        self.assertTrue(np.all(chr1_cov[620:] == 0))
        
        # Check chr2 coverage
        chr2_cov = self.analyzer.coverage_arrays['chr2']
        self.assertTrue(np.all(chr2_cov[100:220] == 1))
        self.assertTrue(np.all(chr2_cov[0:100] == 0))
        self.assertTrue(np.all(chr2_cov[220:] == 0))
    
    def test_depth_statistics(self):
        """Test coverage depth statistics calculation."""
        self.analyzer._load_reference_sequences()
        self.analyzer._parse_psl_file()
        self.analyzer._compute_coverage_arrays()
        self.analyzer._calculate_depth_statistics()
        
        stats = self.analyzer.stats
        
        # Check basic depth statistics exist
        self.assertIn('mean_depth', stats)
        self.assertIn('median_depth', stats)
        self.assertIn('depth_std', stats)
        self.assertIn('max_depth', stats)
        self.assertIn('min_depth', stats)
        
        # Check depth distribution
        self.assertIn('depth_distribution', stats)
        depth_dist = stats['depth_distribution']
        
        # Should have entries for different thresholds
        self.assertIn(1, depth_dist)
        self.assertIn(5, depth_dist)
        
        # Percentage at depth ≥1 should be > 0
        self.assertGreater(depth_dist[1], 0)
        
        # Max depth should be 2 (overlap region)
        self.assertEqual(stats['max_depth'], 2)
    
    def test_breadth_statistics(self):
        """Test coverage breadth statistics calculation."""
        self.analyzer._load_reference_sequences()
        self.analyzer._parse_psl_file()
        self.analyzer._compute_coverage_arrays()
        self.analyzer._calculate_breadth_statistics()
        
        stats = self.analyzer.stats
        
        # Check breadth statistics exist
        self.assertIn('total_bases', stats)
        self.assertIn('covered_bases', stats)
        self.assertIn('coverage_breadth', stats)
        self.assertIn('target_coverage_breadth', stats)
        
        # Total bases should be 2000 + 1500 = 3500
        self.assertEqual(stats['total_bases'], 3500)
        
        # Coverage breadth should be reasonable
        self.assertGreater(stats['coverage_breadth'], 0)
        self.assertLessEqual(stats['coverage_breadth'], 100)
        
        # Target coverage breadth should be lower than general breadth
        self.assertLessEqual(stats['target_coverage_breadth'], stats['coverage_breadth'])
    
    def test_uniformity_statistics(self):
        """Test coverage uniformity statistics."""
        self.analyzer.analyze()
        
        stats = self.analyzer.stats
        
        # Check uniformity metrics exist
        self.assertIn('coverage_cv', stats)
        self.assertIn('coverage_gini', stats)
        self.assertIn('uniformity_score', stats)
        
        # Check reasonable ranges
        self.assertGreaterEqual(stats['coverage_gini'], 0)
        self.assertLessEqual(stats['coverage_gini'], 1)
        
        self.assertGreaterEqual(stats['uniformity_score'], 0)
        self.assertLessEqual(stats['uniformity_score'], 1)
    
    def test_per_reference_statistics(self):
        """Test per-reference sequence statistics."""
        self.analyzer.analyze()
        
        stats = self.analyzer.stats
        per_ref = stats['per_reference']
        
        # Check both chromosomes are present
        self.assertIn('chr1', per_ref)
        self.assertIn('chr2', per_ref)
        
        # Check chr1 statistics
        chr1_stats = per_ref['chr1']
        self.assertEqual(chr1_stats['length'], 2000)
        self.assertGreater(chr1_stats['mean_depth'], 0)
        self.assertEqual(chr1_stats['max_depth'], 2)
        self.assertGreater(chr1_stats['covered_bases'], 0)
        
        # Check chr2 statistics
        chr2_stats = per_ref['chr2']
        self.assertEqual(chr2_stats['length'], 1500)
        self.assertGreater(chr2_stats['mean_depth'], 0)
        self.assertEqual(chr2_stats['max_depth'], 1)
    
    def test_gap_counting(self):
        """Test gap counting functionality."""
        # Create coverage array with known gaps
        coverage_array = np.array([1, 1, 0, 0, 0, 1, 1, 0, 1, 1])
        
        gaps = self.analyzer._count_gaps_in_sequence(coverage_array)
        
        # Should detect 2 gaps: positions 2-4 and position 7
        self.assertEqual(gaps, 2)
        
        # Test edge cases
        no_gaps = np.array([1, 1, 1, 1, 1])
        self.assertEqual(self.analyzer._count_gaps_in_sequence(no_gaps), 0)
        
        all_gaps = np.array([0, 0, 0, 0, 0])
        self.assertEqual(self.analyzer._count_gaps_in_sequence(all_gaps), 1)
    
    def test_gini_coefficient(self):
        """Test Gini coefficient calculation."""
        # Perfect equality (all values equal)
        equal_values = np.array([10, 10, 10, 10, 10])
        gini_equal = self.analyzer._calculate_gini_coefficient(equal_values)
        self.assertAlmostEqual(gini_equal, 0.0, places=2)
        
        # Perfect inequality (one has everything)
        unequal_values = np.array([0, 0, 0, 0, 100])
        gini_unequal = self.analyzer._calculate_gini_coefficient(unequal_values)
        self.assertGreater(gini_unequal, 0.5)
        
        # Empty array
        empty_values = np.array([])
        gini_empty = self.analyzer._calculate_gini_coefficient(empty_values)
        self.assertEqual(gini_empty, 1.0)
    
    def test_export_functionality(self):
        """Test data export functionality."""
        self.analyzer.analyze()
        
        export_dir = self.temp_path / "export_test"
        export_dir.mkdir()
        
        coverage_df, gap_regions = self.analyzer.export_coverage_data(export_dir)
        
        # Check coverage data export
        self.assertIsInstance(coverage_df, pd.DataFrame)
        self.assertIn('chromosome', coverage_df.columns)
        self.assertIn('position', coverage_df.columns)
        self.assertIn('coverage', coverage_df.columns)
        
        # Check files were created
        data_dir = export_dir / "data"
        self.assertTrue(data_dir.exists())
        self.assertTrue((data_dir / "coverage_per_position.csv").exists())
        
        # Check gap regions
        if len(gap_regions) > 0:
            self.assertTrue((data_dir / "gap_regions.bed").exists())
    
    def test_complete_analysis_workflow(self):
        """Test complete analysis workflow."""
        stats = self.analyzer.analyze()
        
        # Check all major statistics are present
        required_keys = [
            'analysis_date', 'parameters', 'total_mappings', 'mapped_oligos',
            'reference_length', 'mean_depth', 'coverage_breadth', 'per_reference'
        ]
        
        for key in required_keys:
            self.assertIn(key, stats)
        
        # Check parameter preservation
        params = stats['parameters']
        self.assertEqual(params['min_coverage'], 1.0)
        self.assertEqual(params['target_coverage'], 5.0)
        self.assertEqual(params['min_identity'], 85.0)
        self.assertEqual(params['min_length'], 100)
        
        # Check reasonable values
        self.assertGreater(stats['total_mappings'], 0)
        self.assertGreater(stats['mapped_oligos'], 0)
        self.assertEqual(stats['reference_length'], 3500)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error conditions."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_empty_psl_file(self):
        """Test handling of empty PSL file."""
        empty_psl = self.temp_path / "empty.psl"
        empty_psl.touch()  # Create empty file
        
        ref_file = self.temp_path / "ref.fasta"
        with open(ref_file, 'w') as f:
            f.write(">chr1\nATGC\n")
        
        analyzer = CoverageAnalyzer(empty_psl, ref_file)
        stats = analyzer.analyze()
        
        # Should handle empty input gracefully
        self.assertEqual(stats['total_mappings'], 0)
        self.assertEqual(stats['mapped_oligos'], 0)
    
    def test_malformed_psl_lines(self):
        """Test handling of malformed PSL lines."""
        malformed_psl = self.temp_path / "malformed.psl"
        with open(malformed_psl, 'w') as f:
            f.write("psLayout version 3\n")
            f.write("------\n")
            f.write("incomplete line\n")  # Too few fields
            f.write("120\tNOT_A_NUMBER\t0\t0\t0\t0\t0\t0\t+\toligo1\t120\t0\t120\tchr1\t1000\t0\t120\t1\t120,\t0,\t0,\n")  # Invalid number
            f.write("120\t0\t0\t0\t0\t0\t0\t0\t+\toligo2\t120\t0\t120\tchr1\t1000\t0\t120\t1\t120,\t0,\t0,\n")  # Valid line
        
        ref_file = self.temp_path / "ref.fasta"
        with open(ref_file, 'w') as f:
            f.write(">chr1\n" + "A" * 1000 + "\n")
        
        analyzer = CoverageAnalyzer(malformed_psl, ref_file)
        stats = analyzer.analyze()
        
        # Should only process valid lines
        self.assertEqual(stats['total_mappings'], 1)
        self.assertEqual(stats['mapped_oligos'], 1)
    
    def test_missing_reference_sequences(self):
        """Test handling when PSL references sequences not in FASTA."""
        psl_file = self.temp_path / "test.psl"
        with open(psl_file, 'w') as f:
            f.write("psLayout version 3\n")
            f.write("------\n")
            f.write("120\t0\t0\t0\t0\t0\t0\t0\t+\toligo1\t120\t0\t120\tmissing_chr\t1000\t0\t120\t1\t120,\t0,\t0,\n")
        
        ref_file = self.temp_path / "ref.fasta"
        with open(ref_file, 'w') as f:
            f.write(">chr1\n" + "A" * 1000 + "\n")
        
        analyzer = CoverageAnalyzer(psl_file, ref_file)
        stats = analyzer.analyze()
        
        # Should handle missing references gracefully
        # Mapping will be recorded but won't contribute to coverage arrays
        self.assertEqual(stats['total_mappings'], 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)