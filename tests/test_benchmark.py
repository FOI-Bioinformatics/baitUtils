#!/usr/bin/env python3

"""
test_benchmark.py

Unit tests for the benchmarking system.
"""

import unittest
import sys
import os
from unittest.mock import patch, MagicMock

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from baitUtils.benchmark import BenchmarkAnalyzer, TheoreticalOptimal, BenchmarkResult
from baitUtils.quality_scorer import QualityScore, QualityCategory


class TestBenchmarkAnalyzer(unittest.TestCase):
    """Test BenchmarkAnalyzer class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.coverage_stats = {
            'coverage_breadth': 75.5,
            'mean_depth': 8.2,
            'gini_coefficient': 0.35,
            'mapping_efficiency': 88.3,
            'reference_length': 12000,
            'covered_bases': 9060
        }
        
        self.gap_analysis = {
            'total_gaps': 35,
            'total_gap_length': 2940,
            'mean_gap_size': 84,
            'max_gap_size': 650
        }
        
        self.reference_analysis = {
            'total_length': 12000,
            'challenging_regions': [
                {'start': 1500, 'end': 1800, 'length': 300, 'type': 'low_complexity'},
                {'start': 6000, 'end': 6200, 'length': 200, 'type': 'homopolymer'},
                {'start': 9500, 'end': 9750, 'length': 250, 'type': 'repeat'}
            ]
        }
        
        self.quality_score = QualityScore(
            overall_score=7.1,
            category=QualityCategory.GOOD,
            component_scores={
                'coverage_score': 7.5,
                'depth_uniformity_score': 6.5,
                'gap_score': 6.8,
                'mapping_efficiency_score': 8.8
            },
            weighted_scores={
                'coverage_score': 7.5,
                'depth_uniformity_score': 6.5,
                'gap_score': 6.8,
                'mapping_efficiency_score': 8.8
            },
            benchmarks={
                'excellent_threshold': 9.0,
                'good_threshold': 7.0,
                'fair_threshold': 5.0
            },
            recommendations=["Consider optimizing coverage distribution", "Review gap filling strategy"]
        )
        
        self.analyzer = BenchmarkAnalyzer(
            self.coverage_stats,
            self.gap_analysis,
            self.reference_analysis,
            self.quality_score
        )
    
    def test_init(self):
        """Test analyzer initialization."""
        self.assertEqual(self.analyzer.coverage_stats, self.coverage_stats)
        self.assertEqual(self.analyzer.gap_analysis, self.gap_analysis)
        self.assertEqual(self.analyzer.reference_analysis, self.reference_analysis)
        self.assertEqual(self.analyzer.quality_score, self.quality_score)
        self.assertEqual(self.analyzer.oligo_length, 120)
        self.assertEqual(self.analyzer.target_spacing, 500)
    
    def test_calculate_theoretical_optimal(self):
        """Test theoretical optimal calculation."""
        theoretical = self.analyzer.calculate_theoretical_optimal()
        
        # Check return type
        self.assertIsInstance(theoretical, TheoreticalOptimal)
        
        # Check that theoretical metrics are reasonable
        self.assertTrue(0 <= theoretical.max_coverage_breadth <= 100)
        self.assertTrue(0 <= theoretical.optimal_depth_uniformity <= 100)
        self.assertGreaterEqual(theoretical.min_possible_gaps, 0)
        self.assertTrue(0 <= theoretical.optimal_mapping_efficiency <= 100)
        self.assertTrue(0 <= theoretical.theoretical_quality_score <= 10)
        
        # Check that theoretical values account for challenging regions
        challenging_bp = sum(region['length'] for region in self.reference_analysis['challenging_regions'])
        expected_max_breadth = (self.reference_analysis['total_length'] - challenging_bp) / self.reference_analysis['total_length'] * 100
        self.assertAlmostEqual(theoretical.max_coverage_breadth, expected_max_breadth, places=1)
        
        # Check that minimum gaps equals number of challenging regions
        self.assertEqual(theoretical.min_possible_gaps, len(self.reference_analysis['challenging_regions']))
    
    def test_benchmark_coverage_breadth(self):
        """Test coverage breadth benchmarking."""
        theoretical = self.analyzer.calculate_theoretical_optimal()
        result = self.analyzer.benchmark_coverage_breadth(theoretical)
        
        # Check return type
        self.assertIsInstance(result, BenchmarkResult)
        
        # Check result values
        self.assertEqual(result.actual_score, self.coverage_stats['coverage_breadth'])
        self.assertEqual(result.theoretical_optimal, theoretical.max_coverage_breadth)
        
        # Check efficiency ratio calculation
        expected_ratio = self.coverage_stats['coverage_breadth'] / theoretical.max_coverage_breadth
        self.assertAlmostEqual(result.efficiency_ratio, expected_ratio, places=3)
        
        # Check improvement potential
        expected_improvement = theoretical.max_coverage_breadth - self.coverage_stats['coverage_breadth']
        self.assertAlmostEqual(result.improvement_potential, expected_improvement, places=1)
        
        # Check category assignment
        self.assertIn(result.category, ['Excellent', 'Good', 'Fair', 'Poor'])
        
        # Check recommendations
        self.assertIsInstance(result.recommendations, list)
    
    def test_benchmark_depth_uniformity(self):
        """Test depth uniformity benchmarking."""
        theoretical = self.analyzer.calculate_theoretical_optimal()
        result = self.analyzer.benchmark_depth_uniformity(theoretical)
        
        # Check return type
        self.assertIsInstance(result, BenchmarkResult)
        
        # Check actual score calculation (conversion from Gini)
        expected_uniformity = (1 - self.coverage_stats['gini_coefficient']) * 100
        self.assertAlmostEqual(result.actual_score, expected_uniformity, places=1)
        
        # Check efficiency ratio
        self.assertTrue(0 <= result.efficiency_ratio <= 1)
        
        # Check category and recommendations
        self.assertIn(result.category, ['Excellent', 'Good', 'Fair', 'Poor'])
        self.assertIsInstance(result.recommendations, list)
    
    def test_benchmark_gap_reduction(self):
        """Test gap reduction benchmarking."""
        theoretical = self.analyzer.calculate_theoretical_optimal()
        result = self.analyzer.benchmark_gap_reduction(theoretical)
        
        # Check return type
        self.assertIsInstance(result, BenchmarkResult)
        
        # Check actual score
        self.assertEqual(result.actual_score, float(self.gap_analysis['total_gaps']))
        self.assertEqual(result.theoretical_optimal, float(theoretical.min_possible_gaps))
        
        # Check efficiency ratio (lower gaps = higher efficiency)
        expected_ratio = max(theoretical.min_possible_gaps / self.gap_analysis['total_gaps'], 0)
        self.assertAlmostEqual(result.efficiency_ratio, expected_ratio, places=3)
        
        # Check improvement potential
        expected_improvement = self.gap_analysis['total_gaps'] - theoretical.min_possible_gaps
        self.assertEqual(result.improvement_potential, expected_improvement)
    
    def test_benchmark_overall_quality(self):
        """Test overall quality benchmarking."""
        theoretical = self.analyzer.calculate_theoretical_optimal()
        result = self.analyzer.benchmark_overall_quality(theoretical)
        
        # Check return type
        self.assertIsInstance(result, BenchmarkResult)
        
        # Check actual score
        self.assertEqual(result.actual_score, self.quality_score.overall_score)
        self.assertEqual(result.theoretical_optimal, theoretical.theoretical_quality_score)
        
        # Check efficiency ratio
        expected_ratio = self.quality_score.overall_score / theoretical.theoretical_quality_score
        self.assertAlmostEqual(result.efficiency_ratio, expected_ratio, places=3)
        
        # Check recommendations include component-specific advice
        self.assertIsInstance(result.recommendations, list)
    
    def test_run_full_benchmark(self):
        """Test complete benchmark analysis."""
        benchmarks = self.analyzer.run_full_benchmark()
        
        # Check return type
        self.assertIsInstance(benchmarks, dict)
        
        # Check all expected metrics are present
        expected_metrics = ['coverage_breadth', 'depth_uniformity', 'gap_reduction', 'overall_quality']
        for metric in expected_metrics:
            self.assertIn(metric, benchmarks)
            self.assertIsInstance(benchmarks[metric], BenchmarkResult)
    
    def test_generate_benchmark_report(self):
        """Test benchmark report generation."""
        benchmarks = self.analyzer.run_full_benchmark()
        report = self.analyzer.generate_benchmark_report(benchmarks)
        
        # Check return type
        self.assertIsInstance(report, str)
        
        # Check report content
        self.assertIn("PERFORMANCE BENCHMARK ANALYSIS", report)
        self.assertIn("Performance vs Theoretical Optimal", report)
        self.assertIn("Improvement Recommendations", report)
        
        # Check all metrics are included in report
        self.assertIn("Coverage Breadth", report)
        self.assertIn("Depth Uniformity", report)
        self.assertIn("Gap Reduction", report)
        self.assertIn("Overall Quality", report)
        
        # Check efficiency percentages are formatted correctly
        for metric_result in benchmarks.values():
            efficiency_pct = f"{metric_result.efficiency_ratio*100:8.1f}%"
            self.assertIn(efficiency_pct.strip(), report)
    
    def test_component_recommendations(self):
        """Test component-specific recommendations."""
        recommendations = self.analyzer._generate_component_recommendations()
        
        # Check return type
        self.assertIsInstance(recommendations, list)
        
        # Test with low component scores
        low_score_quality = QualityScore(
            overall_score=5.0,
            category=QualityCategory.POOR,
            component_scores={
                'coverage_score': 6.0,
                'depth_uniformity_score': 5.0,
                'gap_score': 4.0,
                'mapping_efficiency_score': 6.0
            },
            weighted_scores={
                'coverage_score': 6.0,
                'depth_uniformity_score': 5.0,
                'gap_score': 4.0,
                'mapping_efficiency_score': 6.0
            },
            benchmarks={
                'excellent_threshold': 9.0,
                'good_threshold': 7.0,
                'fair_threshold': 5.0
            },
            recommendations=["Improve overall coverage", "Optimize depth distribution"]
        )
        
        low_analyzer = BenchmarkAnalyzer(
            self.coverage_stats,
            self.gap_analysis,
            self.reference_analysis,
            low_score_quality
        )
        
        low_recommendations = low_analyzer._generate_component_recommendations()
        
        # Should have recommendations for low-scoring components
        self.assertTrue(len(low_recommendations) > 0)
        
        # Check for specific recommendation types
        rec_text = ' '.join(low_recommendations).lower()
        self.assertIn('gap', rec_text)  # Gap score is lowest (4.0)
    
    def test_edge_cases(self):
        """Test edge cases and boundary conditions."""
        # Test with zero gaps
        zero_gap_analysis = self.gap_analysis.copy()
        zero_gap_analysis['total_gaps'] = 0
        
        zero_analyzer = BenchmarkAnalyzer(
            self.coverage_stats,
            zero_gap_analysis,
            self.reference_analysis,
            self.quality_score
        )
        
        theoretical = zero_analyzer.calculate_theoretical_optimal()
        result = zero_analyzer.benchmark_gap_reduction(theoretical)
        
        self.assertEqual(result.efficiency_ratio, 1.0)
        self.assertEqual(result.improvement_potential, 0)
        
        # Test with perfect coverage
        perfect_coverage = self.coverage_stats.copy()
        perfect_coverage['coverage_breadth'] = 100.0
        perfect_coverage['gini_coefficient'] = 0.0
        perfect_coverage['mapping_efficiency'] = 100.0
        
        perfect_analyzer = BenchmarkAnalyzer(
            perfect_coverage,
            zero_gap_analysis,
            self.reference_analysis,
            QualityScore(
                overall_score=10.0,
                category=QualityCategory.EXCELLENT,
                component_scores={},
                weighted_scores={},
                benchmarks={'excellent_threshold': 9.0},
                recommendations=[]
            )
        )
        
        benchmarks = perfect_analyzer.run_full_benchmark()
        
        # All efficiency ratios should be high
        for result in benchmarks.values():
            self.assertGreaterEqual(result.efficiency_ratio, 0.8)


class TestTheoreticalOptimal(unittest.TestCase):
    """Test TheoreticalOptimal NamedTuple."""
    
    def test_theoretical_optimal_creation(self):
        """Test TheoreticalOptimal creation and access."""
        theoretical = TheoreticalOptimal(
            max_coverage_breadth=95.5,
            optimal_depth_uniformity=92.0,
            min_possible_gaps=3,
            optimal_mapping_efficiency=96.0,
            theoretical_quality_score=9.2
        )
        
        self.assertEqual(theoretical.max_coverage_breadth, 95.5)
        self.assertEqual(theoretical.optimal_depth_uniformity, 92.0)
        self.assertEqual(theoretical.min_possible_gaps, 3)
        self.assertEqual(theoretical.optimal_mapping_efficiency, 96.0)
        self.assertEqual(theoretical.theoretical_quality_score, 9.2)


class TestBenchmarkResult(unittest.TestCase):
    """Test BenchmarkResult NamedTuple."""
    
    def test_benchmark_result_creation(self):
        """Test BenchmarkResult creation and access."""
        result = BenchmarkResult(
            actual_score=75.5,
            theoretical_optimal=90.0,
            efficiency_ratio=0.839,
            improvement_potential=14.5,
            category='Good',
            recommendations=['Add more oligos', 'Improve uniformity']
        )
        
        self.assertEqual(result.actual_score, 75.5)
        self.assertEqual(result.theoretical_optimal, 90.0)
        self.assertAlmostEqual(result.efficiency_ratio, 0.839, places=3)
        self.assertEqual(result.improvement_potential, 14.5)
        self.assertEqual(result.category, 'Good')
        self.assertEqual(len(result.recommendations), 2)


if __name__ == '__main__':
    unittest.main()