#!/usr/bin/env python3

"""
test_quality_scorer.py

Unit tests for the quality scoring system.
"""

import unittest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from baitUtils.quality_scorer import QualityScorer, QualityScore


class TestQualityScorer(unittest.TestCase):
    """Test QualityScorer class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.coverage_stats = {
            'coverage_breadth': 82.5,
            'mean_depth': 9.7,
            'gini_coefficient': 0.28,
            'mapping_efficiency': 91.2,
            'reference_length': 15000,
            'covered_bases': 12375,
            'depth_distribution': {1: 82.5, 5: 68.2, 10: 45.8, 20: 23.1}
        }
        
        self.gap_analysis = {
            'total_gaps': 28,
            'total_gap_length': 2625,
            'mean_gap_size': 93.8,
            'max_gap_size': 750,
            'gap_percentage': 17.5
        }
        
        self.reference_analysis = {
            'total_length': 15000,
            'num_sequences': 3,
            'challenging_regions': [
                {'start': 2000, 'end': 2250, 'length': 250, 'type': 'low_complexity'},
                {'start': 7500, 'end': 7750, 'length': 250, 'type': 'homopolymer'}
            ],
            'sequence_features': {
                'mean_gc_content': 0.42,
                'mean_complexity': 0.85
            }
        }
        
        self.scorer = QualityScorer(
            self.coverage_stats,
            self.gap_analysis,
            self.reference_analysis
        )
    
    def test_init(self):
        """Test scorer initialization."""
        self.assertEqual(self.scorer.coverage_stats, self.coverage_stats)
        self.assertEqual(self.scorer.gap_analysis, self.gap_analysis)
        self.assertEqual(self.scorer.reference_analysis, self.reference_analysis)
        
        # Check default weights
        expected_weights = {'coverage': 0.3, 'depth': 0.25, 'gaps': 0.25, 'mapping': 0.2}
        self.assertEqual(self.scorer.weights, expected_weights)
    
    def test_init_custom_weights(self):
        """Test scorer initialization with custom weights."""
        custom_weights = {'coverage': 0.4, 'depth': 0.3, 'gaps': 0.2, 'mapping': 0.1}
        scorer = QualityScorer(
            self.coverage_stats,
            self.gap_analysis,
            self.reference_analysis,
            weights=custom_weights
        )
        
        self.assertEqual(scorer.weights, custom_weights)
    
    def test_score_coverage_breadth(self):
        """Test coverage breadth scoring."""
        score = self.scorer._score_coverage_breadth()
        
        # Check score is in valid range
        self.assertTrue(0 <= score <= 10)
        
        # Check score calculation
        # 82.5% should give a good score
        self.assertGreater(score, 7.0)
        
        # Test edge cases
        perfect_stats = self.coverage_stats.copy()
        perfect_stats['coverage_breadth'] = 100.0
        perfect_scorer = QualityScorer(perfect_stats, self.gap_analysis, self.reference_analysis)
        perfect_score = perfect_scorer._score_coverage_breadth()
        self.assertAlmostEqual(perfect_score, 10.0, places=1)
        
        poor_stats = self.coverage_stats.copy()
        poor_stats['coverage_breadth'] = 30.0
        poor_scorer = QualityScorer(poor_stats, self.gap_analysis, self.reference_analysis)
        poor_score = poor_scorer._score_coverage_breadth()
        self.assertLess(poor_score, 5.0)
    
    def test_score_depth_uniformity(self):
        """Test depth uniformity scoring."""
        score = self.scorer._score_depth_uniformity()
        
        # Check score is in valid range
        self.assertTrue(0 <= score <= 10)
        
        # Check that lower Gini coefficient gives higher score
        uniform_stats = self.coverage_stats.copy()
        uniform_stats['gini_coefficient'] = 0.1  # More uniform
        uniform_scorer = QualityScorer(uniform_stats, self.gap_analysis, self.reference_analysis)
        uniform_score = uniform_scorer._score_depth_uniformity()
        
        self.assertGreater(uniform_score, score)
        
        # Test edge cases
        variable_stats = self.coverage_stats.copy()
        variable_stats['gini_coefficient'] = 0.8  # Highly variable
        variable_scorer = QualityScorer(variable_stats, self.gap_analysis, self.reference_analysis)
        variable_score = variable_scorer._score_depth_uniformity()
        
        self.assertLess(variable_score, 3.0)
    
    def test_score_gap_analysis(self):
        """Test gap analysis scoring."""
        score = self.scorer._score_gap_analysis()
        
        # Check score is in valid range
        self.assertTrue(0 <= score <= 10)
        
        # Test with fewer gaps (should score higher)
        good_gaps = self.gap_analysis.copy()
        good_gaps['total_gaps'] = 5
        good_gaps['gap_percentage'] = 3.3
        good_scorer = QualityScorer(self.coverage_stats, good_gaps, self.reference_analysis)
        good_score = good_scorer._score_gap_analysis()
        
        self.assertGreater(good_score, score)
        
        # Test with many gaps (should score lower)
        bad_gaps = self.gap_analysis.copy()
        bad_gaps['total_gaps'] = 200
        bad_gaps['gap_percentage'] = 45.0
        bad_scorer = QualityScorer(self.coverage_stats, bad_gaps, self.reference_analysis)
        bad_score = bad_scorer._score_gap_analysis()
        
        self.assertLess(bad_score, score)
    
    def test_score_mapping_efficiency(self):
        """Test mapping efficiency scoring."""
        score = self.scorer._score_mapping_efficiency()
        
        # Check score is in valid range
        self.assertTrue(0 <= score <= 10)
        
        # 91.2% efficiency should give a good score
        self.assertGreater(score, 8.0)
        
        # Test edge cases
        perfect_stats = self.coverage_stats.copy()
        perfect_stats['mapping_efficiency'] = 100.0
        perfect_scorer = QualityScorer(perfect_stats, self.gap_analysis, self.reference_analysis)
        perfect_score = perfect_scorer._score_mapping_efficiency()
        self.assertAlmostEqual(perfect_score, 10.0, places=1)
        
        poor_stats = self.coverage_stats.copy()
        poor_stats['mapping_efficiency'] = 45.0
        poor_scorer = QualityScorer(poor_stats, self.gap_analysis, self.reference_analysis)
        poor_score = poor_scorer._score_mapping_efficiency()
        self.assertLess(poor_score, 5.0)
    
    def test_calculate_score(self):
        """Test overall score calculation."""
        quality_score = self.scorer.calculate_score()
        
        # Check return type
        self.assertIsInstance(quality_score, QualityScore)
        
        # Check overall score is in valid range
        self.assertTrue(0 <= quality_score.overall_score <= 10)
        
        # Check grade assignment
        self.assertIn(quality_score.grade, ['A', 'B', 'C', 'D', 'F'])
        
        # Check component scores exist
        self.assertIsInstance(quality_score.component_scores, dict)
        expected_components = ['coverage_score', 'depth_uniformity_score', 
                              'gap_score', 'mapping_efficiency_score']
        for component in expected_components:
            self.assertIn(component, quality_score.component_scores)
            self.assertTrue(0 <= quality_score.component_scores[component] <= 10)
        
        # Check grade assignment logic
        if quality_score.overall_score >= 9:
            self.assertEqual(quality_score.grade, 'A')
        elif quality_score.overall_score >= 8:
            self.assertEqual(quality_score.grade, 'B')
        elif quality_score.overall_score >= 7:
            self.assertEqual(quality_score.grade, 'C')
        elif quality_score.overall_score >= 6:
            self.assertEqual(quality_score.grade, 'D')
        else:
            self.assertEqual(quality_score.grade, 'F')
    
    def test_weighted_average_calculation(self):
        """Test that weighted average is calculated correctly."""
        quality_score = self.scorer.calculate_score()
        
        # Calculate expected weighted average
        component_scores = quality_score.component_scores
        weights = self.scorer.weights
        
        expected_score = (
            component_scores['coverage_score'] * weights['coverage'] +
            component_scores['depth_uniformity_score'] * weights['depth'] +
            component_scores['gap_score'] * weights['gaps'] +
            component_scores['mapping_efficiency_score'] * weights['mapping']
        )
        
        self.assertAlmostEqual(quality_score.overall_score, expected_score, places=2)
    
    def test_grade_boundaries(self):
        """Test grade boundary assignments."""
        # Test each grade boundary
        test_cases = [
            (9.5, 'A'),
            (8.5, 'B'),
            (7.5, 'C'),
            (6.5, 'D'),
            (5.0, 'F'),
            (2.0, 'F')
        ]
        
        for score, expected_grade in test_cases:
            grade = self.scorer._assign_grade(score)
            self.assertEqual(grade, expected_grade)
    
    def test_edge_cases(self):
        """Test edge cases and boundary conditions."""
        # Test with zero coverage
        zero_stats = {
            'coverage_breadth': 0.0,
            'mean_depth': 0.0,
            'gini_coefficient': 1.0,
            'mapping_efficiency': 0.0,
            'reference_length': 15000,
            'covered_bases': 0
        }
        
        zero_gaps = {
            'total_gaps': 1000,
            'total_gap_length': 15000,
            'mean_gap_size': 15.0,
            'max_gap_size': 15000,
            'gap_percentage': 100.0
        }
        
        zero_scorer = QualityScorer(zero_stats, zero_gaps, self.reference_analysis)
        zero_quality = zero_scorer.calculate_score()
        
        # Should get very low score
        self.assertLess(zero_quality.overall_score, 3.0)
        self.assertEqual(zero_quality.grade, 'F')
        
        # Test with perfect metrics
        perfect_stats = {
            'coverage_breadth': 100.0,
            'mean_depth': 20.0,
            'gini_coefficient': 0.0,
            'mapping_efficiency': 100.0,
            'reference_length': 15000,
            'covered_bases': 15000
        }
        
        perfect_gaps = {
            'total_gaps': 0,
            'total_gap_length': 0,
            'mean_gap_size': 0,
            'max_gap_size': 0,
            'gap_percentage': 0.0
        }
        
        perfect_scorer = QualityScorer(perfect_stats, perfect_gaps, self.reference_analysis)
        perfect_quality = perfect_scorer.calculate_score()
        
        # Should get high score
        self.assertGreater(perfect_quality.overall_score, 8.5)
        self.assertIn(perfect_quality.grade, ['A', 'B'])
    
    def test_missing_data_handling(self):
        """Test handling of missing data."""
        # Test with missing keys
        incomplete_stats = {
            'coverage_breadth': 80.0,
            'mapping_efficiency': 85.0
            # Missing mean_depth, gini_coefficient, etc.
        }
        
        incomplete_scorer = QualityScorer(incomplete_stats, self.gap_analysis, self.reference_analysis)
        
        # Should handle missing data gracefully with defaults
        score = incomplete_scorer._score_depth_uniformity()
        self.assertTrue(0 <= score <= 10)
        
        coverage_score = incomplete_scorer._score_coverage_breadth()
        self.assertTrue(0 <= coverage_score <= 10)


class TestQualityScore(unittest.TestCase):
    """Test QualityScore dataclass."""
    
    def test_quality_score_creation(self):
        """Test QualityScore creation and access."""
        component_scores = {
            'coverage_score': 8.2,
            'depth_uniformity_score': 7.1,
            'gap_score': 6.8,
            'mapping_efficiency_score': 9.1
        }
        
        quality_score = QualityScore(
            overall_score=7.8,
            grade='B',
            component_scores=component_scores
        )
        
        self.assertEqual(quality_score.overall_score, 7.8)
        self.assertEqual(quality_score.grade, 'B')
        self.assertEqual(quality_score.component_scores, component_scores)
    
    def test_quality_score_string_representation(self):
        """Test QualityScore string representation."""
        quality_score = QualityScore(
            overall_score=8.5,
            grade='A',
            component_scores={'coverage_score': 9.0}
        )
        
        str_repr = str(quality_score)
        self.assertIn('8.5', str_repr)
        self.assertIn('A', str_repr)


if __name__ == '__main__':
    unittest.main()