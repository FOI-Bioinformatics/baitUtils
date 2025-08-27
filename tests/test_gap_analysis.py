#!/usr/bin/env python3

"""
test_gap_analysis.py

Unit tests for the gap_analysis module.
Tests sequence feature analysis and gap characterization functions.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import numpy as np

from baitUtils.gap_analysis import GapAnalyzer


class TestSequenceFeatureAnalysis(unittest.TestCase):
    """Test sequence feature calculation functions."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
        
        # Create test reference file
        self.ref_file = self.temp_path / "reference.fasta"
        with open(self.ref_file, 'w') as f:
            f.write(">chr1\n")
            # Mix of different sequence features
            f.write("ATGCGCGCGC" +         # Normal GC content
                   "AAAAAAAAAA" +         # AT-rich homopolymer
                   "NNNNNNNNNN" +         # N's
                   "GCGCGCGCGC" +         # GC-rich
                   "ATATATATATAT" +       # Low complexity
                   "CGTACGTACGTA" +       # Higher complexity
                   "\n")
        
        # Mock coverage data
        self.coverage_data = {
            'reference_length': 500,
            'per_reference': {
                'chr1': {
                    'length': 72,  # Length of sequence above
                    'coverage_breadth': 50.0,
                    'gaps': 2
                }
            }
        }
        
        self.analyzer = GapAnalyzer(
            coverage_data=self.coverage_data,
            reference_file=self.ref_file,
            min_gap_size=20
        )
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_gc_content_calculation(self):
        """Test GC content calculation in sequence features."""
        # Test different GC contents
        at_rich = "AAATTTAAATTT"  # 0% GC
        balanced = "ATGCATGCATGC"  # 50% GC  
        gc_rich = "GCGCGCGCGCGC"  # 100% GC
        
        features_at = self.analyzer._calculate_sequence_features(at_rich)
        features_balanced = self.analyzer._calculate_sequence_features(balanced)
        features_gc = self.analyzer._calculate_sequence_features(gc_rich)
        
        self.assertAlmostEqual(features_at['gc_content'], 0.0, places=1)
        self.assertAlmostEqual(features_balanced['gc_content'], 50.0, places=1)
        self.assertAlmostEqual(features_gc['gc_content'], 100.0, places=1)
        
        # Test boolean flags
        self.assertTrue(features_at['at_rich'])
        self.assertFalse(features_balanced['at_rich'])
        self.assertFalse(features_gc['at_rich'])
        
        self.assertFalse(features_at['gc_rich'])
        self.assertFalse(features_balanced['gc_rich'])
        self.assertTrue(features_gc['gc_rich'])
    
    def test_shannon_entropy_calculation(self):
        """Test Shannon entropy calculation."""
        # Maximum entropy: equal distribution of all 4 bases
        max_entropy_seq = "ATGC" * 25  # 100 bases
        entropy_max = self.analyzer._calculate_shannon_entropy(max_entropy_seq)
        self.assertAlmostEqual(entropy_max, 2.0, places=1)  # log2(4) = 2
        
        # Minimum entropy: single nucleotide
        min_entropy_seq = "AAAA"
        entropy_min = self.analyzer._calculate_shannon_entropy(min_entropy_seq)
        self.assertEqual(entropy_min, 0.0)
        
        # Medium entropy: two nucleotides
        medium_entropy_seq = "ATAT" * 25
        entropy_medium = self.analyzer._calculate_shannon_entropy(medium_entropy_seq)
        self.assertAlmostEqual(entropy_medium, 1.0, places=1)  # log2(2) = 1
        
        # Empty sequence
        entropy_empty = self.analyzer._calculate_shannon_entropy("")
        self.assertEqual(entropy_empty, 0.0)
    
    def test_homopolymer_detection(self):
        """Test homopolymer run detection."""
        # No homopolymers
        no_homo = "ATGCATGC"
        max_run = self.analyzer._find_max_homopolymer(no_homo)
        self.assertEqual(max_run, 1)
        
        # Short homopolymer
        short_homo = "ATGCCCCGAT"
        max_run = self.analyzer._find_max_homopolymer(short_homo)
        self.assertEqual(max_run, 4)
        
        # Long homopolymer
        long_homo = "ATGCCCCCCCCCGAT"
        max_run = self.analyzer._find_max_homopolymer(long_homo)
        self.assertEqual(max_run, 10)
        
        # Homopolymer at start
        start_homo = "AAAAATGC"
        max_run = self.analyzer._find_max_homopolymer(start_homo)
        self.assertEqual(max_run, 5)
        
        # Homopolymer at end
        end_homo = "ATGCCCCC"
        max_run = self.analyzer._find_max_homopolymer(end_homo)
        self.assertEqual(max_run, 5)
        
        # Empty sequence
        empty_seq = ""
        max_run = self.analyzer._find_max_homopolymer(empty_seq)
        self.assertEqual(max_run, 0)
    
    def test_n_content_calculation(self):
        """Test N content calculation."""
        # No N's
        no_n = "ATGCATGC"
        features = self.analyzer._calculate_sequence_features(no_n)
        self.assertEqual(features['n_content'], 0.0)
        
        # Some N's
        some_n = "ATGCNNNATGC"  # 3 N's out of 11 bases
        features = self.analyzer._calculate_sequence_features(some_n)
        expected_n_content = (3 / 11) * 100
        self.assertAlmostEqual(features['n_content'], expected_n_content, places=1)
        
        # All N's
        all_n = "NNNNNNNN"
        features = self.analyzer._calculate_sequence_features(all_n)
        self.assertEqual(features['n_content'], 100.0)
    
    def test_repeat_content_estimation(self):
        """Test repetitive content estimation."""
        # No repeats
        no_repeats = "ATGCTAGCTA"
        repeat_content = self.analyzer._calculate_repeat_content(no_repeats)
        self.assertEqual(repeat_content, 0.0)
        
        # Simple tandem repeat
        tandem_repeat = "ATGATGATGATG"
        repeat_content = self.analyzer._calculate_repeat_content(tandem_repeat)
        self.assertGreater(repeat_content, 0.0)
        
        # Sequence too short for analysis
        too_short = "ATGC"
        repeat_content = self.analyzer._calculate_repeat_content(too_short)
        self.assertEqual(repeat_content, 0.0)
    
    def test_low_complexity_regions(self):
        """Test low complexity region detection."""
        # High complexity sequence
        high_complexity = "ATGCTAGCATGCTAGCATGC"
        low_complex_fraction = self.analyzer._calculate_low_complexity(high_complexity)
        self.assertLessEqual(low_complex_fraction, 0.5)
        
        # Low complexity sequence (homopolymer)
        low_complexity = "AAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        low_complex_fraction = self.analyzer._calculate_low_complexity(low_complexity)
        self.assertGreater(low_complex_fraction, 0.5)
        
        # Sequence too short for analysis
        too_short = "ATGC"
        low_complex_fraction = self.analyzer._calculate_low_complexity(too_short)
        self.assertEqual(low_complex_fraction, 0.0)
    
    def test_dinucleotide_bias(self):
        """Test dinucleotide composition bias calculation."""
        # Uniform dinucleotide composition
        uniform_seq = "ATGCATGCATGC"
        bias = self.analyzer._calculate_dinucleotide_bias(uniform_seq)
        self.assertIsInstance(bias, float)
        
        # Biased composition (only AT dinucleotides)
        biased_seq = "ATATATATATATAT"
        bias_biased = self.analyzer._calculate_dinucleotide_bias(biased_seq)
        self.assertGreater(bias_biased, bias)
        
        # Sequence with N's (should be ignored)
        n_seq = "ATNGCNATNGC"
        bias_n = self.analyzer._calculate_dinucleotide_bias(n_seq)
        self.assertIsInstance(bias_n, float)
        
        # Too short sequence
        short_seq = "A"
        bias_short = self.analyzer._calculate_dinucleotide_bias(short_seq)
        self.assertEqual(bias_short, 0.0)
    
    def test_comprehensive_feature_analysis(self):
        """Test comprehensive sequence feature calculation."""
        test_sequence = (
            "ATGCGCGCGC" +      # Normal sequence
            "AAAAAAAAAA" +      # Homopolymer
            "NNNNN" +           # N's
            "GCGCGCGCGC" +      # GC-rich
            "ATATATATAT"        # Low complexity
        )
        
        features = self.analyzer._calculate_sequence_features(test_sequence)
        
        # Check all expected features are present
        expected_features = [
            'gc_content', 'complexity', 'repeat_content', 'max_homopolymer',
            'n_content', 'low_complexity_fraction', 'at_rich', 'gc_rich',
            'extreme_composition', 'dinuc_bias'
        ]
        
        for feature in expected_features:
            self.assertIn(feature, features)
        
        # Check reasonable ranges
        self.assertGreaterEqual(features['gc_content'], 0)
        self.assertLessEqual(features['gc_content'], 100)
        
        self.assertGreaterEqual(features['complexity'], 0)
        self.assertLessEqual(features['complexity'], 2)  # Max for DNA
        
        self.assertGreaterEqual(features['n_content'], 0)
        self.assertLessEqual(features['n_content'], 100)
        
        # Should detect homopolymer
        self.assertGreaterEqual(features['max_homopolymer'], 10)
        
        # Should detect N content
        self.assertGreater(features['n_content'], 0)


class TestGapAnalysisWorkflow(unittest.TestCase):
    """Test complete gap analysis workflow."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
        
        # Create complex reference file
        self.ref_file = self.temp_path / "complex_reference.fasta"
        with open(self.ref_file, 'w') as f:
            f.write(">chr1\n")
            # Create sequence with different problematic regions
            f.write("ATGCATGCATGC" * 10 +           # Normal region (120bp)
                   "AAAAAAAAAA" * 10 +            # AT-rich homopolymers (100bp)
                   "NNNNNNNNNN" * 5 +             # N region (50bp)
                   "GCGCGCGCGC" * 15 +            # GC-rich region (150bp)
                   "CGTACGTACGTA" * 20 +          # Complex region (240bp)
                   "\n")
            
            f.write(">chr2\n")
            f.write("ATGC" * 200 + "\n")  # Simple sequence (800bp)
        
        # Mock coverage data with multiple gaps
        self.coverage_data = {
            'reference_length': 1460,  # 660 + 800
            'per_reference': {
                'chr1': {
                    'length': 660,
                    'coverage_breadth': 60.0,
                    'gaps': 4
                },
                'chr2': {
                    'length': 800,
                    'coverage_breadth': 85.0,
                    'gaps': 2
                }
            }
        }
        
        self.analyzer = GapAnalyzer(
            coverage_data=self.coverage_data,
            reference_file=self.ref_file,
            min_gap_size=50,
            extend_bp=10
        )
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_reference_sequence_loading(self):
        """Test reference sequence loading."""
        self.analyzer._load_reference_sequences()
        
        self.assertIn('chr1', self.analyzer.reference_sequences)
        self.assertIn('chr2', self.analyzer.reference_sequences)
        
        # Check sequence lengths
        self.assertEqual(len(self.analyzer.reference_sequences['chr1']), 660)
        self.assertEqual(len(self.analyzer.reference_sequences['chr2']), 800)
    
    def test_gap_extraction(self):
        """Test gap extraction from coverage data."""
        self.analyzer._load_reference_sequences()
        self.analyzer._extract_gaps()
        
        # Should extract gaps based on mock data
        self.assertGreater(len(self.analyzer.gaps), 0)
        
        # Check gap structure
        for gap in self.analyzer.gaps:
            self.assertIn('id', gap)
            self.assertIn('chromosome', gap)
            self.assertIn('start', gap)
            self.assertIn('end', gap)
            self.assertIn('length', gap)
            self.assertIn('extended_start', gap)
            self.assertIn('extended_end', gap)
            
            # Check gap size filter
            self.assertGreaterEqual(gap['length'], self.analyzer.min_gap_size)
            
            # Check extension
            self.assertLessEqual(gap['extended_start'], gap['start'])
            self.assertGreaterEqual(gap['extended_end'], gap['end'])
    
    def test_gap_feature_analysis(self):
        """Test gap sequence feature analysis."""
        self.analyzer._load_reference_sequences()
        self.analyzer._extract_gaps()
        self.analyzer._analyze_gap_features()
        
        # Should have analyzed features for each gap
        self.assertEqual(len(self.analyzer.gap_features), len(self.analyzer.gaps))
        
        # Check feature structure
        if self.analyzer.gap_features:
            features = self.analyzer.gap_features[0]
            expected_keys = [
                'gap_id', 'chromosome', 'start', 'end', 'length',
                'gc_content', 'complexity', 'max_homopolymer', 'n_content'
            ]
            
            for key in expected_keys:
                self.assertIn(key, features)
    
    def test_gap_statistics_computation(self):
        """Test gap statistics computation."""
        self.analyzer._load_reference_sequences()
        self.analyzer._extract_gaps()
        self.analyzer._compute_gap_statistics()
        
        stats = self.analyzer.analysis_results
        
        # Check basic gap statistics
        expected_keys = [
            'total_gaps', 'total_gap_length', 'gap_percentage',
            'mean_gap_size', 'median_gap_size', 'max_gap_size',
            'size_distribution', 'largest_gaps'
        ]
        
        for key in expected_keys:
            self.assertIn(key, stats)
        
        # Check reasonable values
        self.assertGreaterEqual(stats['total_gaps'], 0)
        self.assertGreaterEqual(stats['total_gap_length'], 0)
        self.assertGreaterEqual(stats['gap_percentage'], 0)
        self.assertLessEqual(stats['gap_percentage'], 100)
        
        # Check size distribution structure
        size_dist = stats['size_distribution']
        self.assertIsInstance(size_dist, dict)
        
        # Check largest gaps structure
        largest_gaps = stats['largest_gaps']
        self.assertIsInstance(largest_gaps, list)
    
    def test_feature_correlation_analysis(self):
        """Test gap feature correlation analysis."""
        self.analyzer.analyze()  # Full analysis
        
        stats = self.analyzer.analysis_results
        
        if 'feature_analysis' in stats:
            feature_analysis = stats['feature_analysis']
            
            # Check expected feature categories
            expected_categories = ['gc_content', 'complexity', 'repeats', 'homopolymers']
            
            for category in expected_categories:
                if category in feature_analysis:
                    category_stats = feature_analysis[category]
                    self.assertIsInstance(category_stats, dict)
                    
                    # Check for mean values
                    mean_key = f'mean_{category}' if category != 'gc_content' else 'mean_gc'
                    if mean_key in category_stats:
                        self.assertIsInstance(category_stats[mean_key], (int, float))
    
    def test_suggestion_generation(self):
        """Test improvement suggestion generation."""
        self.analyzer.analyze()
        
        stats = self.analyzer.analysis_results
        
        self.assertIn('suggestions', stats)
        suggestions = stats['suggestions']
        
        self.assertIsInstance(suggestions, list)
        self.assertGreater(len(suggestions), 0)
        
        # Each suggestion should be a string
        for suggestion in suggestions:
            self.assertIsInstance(suggestion, str)
            self.assertGreater(len(suggestion), 0)
    
    def test_complete_analysis_workflow(self):
        """Test complete gap analysis workflow."""
        results = self.analyzer.analyze()
        
        # Check all major components are present
        expected_keys = [
            'total_gaps', 'gap_percentage', 'size_distribution',
            'largest_gaps', 'suggestions'
        ]
        
        for key in expected_keys:
            self.assertIn(key, results)
        
        # Check data types
        self.assertIsInstance(results['total_gaps'], int)
        self.assertIsInstance(results['gap_percentage'], (int, float))
        self.assertIsInstance(results['size_distribution'], dict)
        self.assertIsInstance(results['largest_gaps'], list)
        self.assertIsInstance(results['suggestions'], list)


class TestEdgeCasesGapAnalysis(unittest.TestCase):
    """Test edge cases in gap analysis."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
    
    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir)
    
    def test_no_gaps_scenario(self):
        """Test analysis when there are no gaps."""
        ref_file = self.temp_path / "ref.fasta"
        with open(ref_file, 'w') as f:
            f.write(">chr1\nATGCATGC\n")
        
        # Perfect coverage (no gaps)
        coverage_data = {
            'reference_length': 8,
            'per_reference': {
                'chr1': {
                    'length': 8,
                    'coverage_breadth': 100.0,
                    'gaps': 0
                }
            }
        }
        
        analyzer = GapAnalyzer(coverage_data, ref_file, min_gap_size=50)
        results = analyzer.analyze()
        
        # Should handle no gaps gracefully
        self.assertEqual(results['total_gaps'], 0)
        self.assertEqual(results['total_gap_length'], 0)
        self.assertEqual(results['gap_percentage'], 0.0)
        self.assertIn('suggestions', results)
    
    def test_all_gaps_scenario(self):
        """Test analysis when coverage is very poor."""
        ref_file = self.temp_path / "ref.fasta"
        with open(ref_file, 'w') as f:
            f.write(">chr1\n" + "A" * 1000 + "\n")
        
        # Very poor coverage
        coverage_data = {
            'reference_length': 1000,
            'per_reference': {
                'chr1': {
                    'length': 1000,
                    'coverage_breadth': 5.0,
                    'gaps': 1
                }
            }
        }
        
        analyzer = GapAnalyzer(coverage_data, ref_file, min_gap_size=100)
        results = analyzer.analyze()
        
        # Should handle poor coverage
        self.assertGreater(results['gap_percentage'], 80.0)
        self.assertGreater(len(results['suggestions']), 1)
    
    def test_empty_sequence_features(self):
        """Test feature calculation on empty sequences."""
        ref_file = self.temp_path / "empty_ref.fasta"
        with open(ref_file, 'w') as f:
            f.write(">chr1\n\n")  # Empty sequence
        
        coverage_data = {
            'reference_length': 0,
            'per_reference': {}
        }
        
        analyzer = GapAnalyzer(coverage_data, ref_file)
        
        # Test empty sequence feature calculation
        features = analyzer._calculate_sequence_features("")
        
        expected_defaults = {
            'gc_content': 0, 'complexity': 0, 'repeats': 0,
            'homopolymers': 0, 'n_content': 0
        }
        
        for key, expected_value in expected_defaults.items():
            if key in features:
                self.assertEqual(features[key], expected_value)


if __name__ == '__main__':
    unittest.main(verbosity=2)