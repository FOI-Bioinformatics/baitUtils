#!/usr/bin/env python3

"""
test_phase2_integration.py

Integration tests for Phase 2 features of the baitUtils evaluate command.
Tests interactive reporting, reference analysis, quality scoring, and benchmarking.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from baitUtils.reference_analyzer import ReferenceAnalyzer
from baitUtils.quality_scorer import QualityScorer, QualityScore, QualityCategory
from baitUtils.report_generator import InteractiveReportGenerator
from baitUtils.interactive_plots import InteractivePlotter
from baitUtils.benchmark import BenchmarkAnalyzer, TheoreticalOptimal, BenchmarkResult


def create_quality_score(score, grade_letter):
    """Helper function to create QualityScore instances with proper API."""
    category_map = {
        'A': QualityCategory.EXCELLENT,
        'B': QualityCategory.GOOD,
        'C': QualityCategory.FAIR,
        'D': QualityCategory.POOR
    }
    
    return QualityScore(
        overall_score=score,
        category=category_map.get(grade_letter, QualityCategory.FAIR),
        component_scores={
            'coverage_breadth': score,
            'coverage_depth': score,
            'gap_characteristics': score,
            'mapping_efficiency': score,
            'reference_difficulty': score
        },
        weighted_scores={
            'coverage_breadth': score,
            'coverage_depth': score,
            'gap_characteristics': score,
            'mapping_efficiency': score,
            'reference_difficulty': score
        },
        benchmarks={
            'excellent_threshold': 9.0,
            'good_threshold': 7.0,
            'fair_threshold': 5.0
        },
        recommendations=[]
    )


class TestReferenceAnalyzer(unittest.TestCase):
    """Test reference sequence analysis functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.ref_file = self.test_dir / "test_reference.fasta"
        
        # Create mock reference file
        with open(self.ref_file, 'w') as f:
            f.write(">chr1\n")
            f.write("ATCGATCGATCG" * 100 + "\n")  # 1200 bp
            f.write(">chr2\n") 
            f.write("GCTAGCTAGCTA" * 50 + "\n")   # 600 bp
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    @patch('baitUtils.reference_analyzer.SeqIO')
    def test_analyze_reference(self, mock_seqio):
        """Test reference sequence analysis."""
        # Mock SeqIO.parse
        mock_record1 = MagicMock()
        mock_record1.id = "chr1"
        mock_record1.seq = "ATCGATCGATCG" * 100
        mock_record2 = MagicMock()
        mock_record2.id = "chr2" 
        mock_record2.seq = "GCTAGCTAGCTA" * 50
        mock_seqio.parse.return_value = [mock_record1, mock_record2]
        
        coverage_data = {'coverage_breadth': 85.0, 'mean_depth': 10.5}
        analyzer = ReferenceAnalyzer(self.ref_file, coverage_data, window_size=100)
        results = analyzer.analyze()
        
        # Check results structure
        self.assertIn('analysis_summary', results)
        self.assertIn('challenging_regions', results)
        self.assertIn('sequence_features', results)
        
        # Check basic statistics
        self.assertEqual(results['analysis_summary']['total_sequences'], 2)
    
    def test_calculate_sequence_features(self):
        """Test sequence feature calculations."""
        coverage_data = {'coverage_breadth': 85.0, 'mean_depth': 10.5}
        analyzer = ReferenceAnalyzer(self.ref_file, coverage_data)
        
        test_seq = "ATCGATCGATCGAAAAAATTTTTCCCCCGGGGG"
        # Test the actual analyze method instead since _calculate_sequence_features is not public
        results = analyzer.analyze()
        features = results['sequence_features']
        
        # Check feature structure - features are organized by chromosome
        self.assertIn('chr1', features)
        self.assertIn('chr2', features)
        
        # Check feature keys for chr1
        chr1_features = features['chr1']
        expected_keys = ['gc_content', 'at_content', 'n_content', 'length']
        for key in expected_keys:
            self.assertIn(key, chr1_features)
        
        # Check GC content calculation for chr1
        self.assertAlmostEqual(chr1_features['gc_content'], 50.0, places=1)
        
        # Check length
        self.assertEqual(chr1_features['length'], 1200)


class TestQualityScorer(unittest.TestCase):
    """Test quality scoring system."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.coverage_stats = {
            'coverage_breadth': 85.5,
            'mean_depth': 12.3,
            'gini_coefficient': 0.3,
            'mapping_efficiency': 92.1,
            'reference_length': 10000,
            'covered_bases': 8550
        }
        
        self.gap_analysis = {
            'total_gaps': 25,
            'total_gap_length': 1450,
            'mean_gap_size': 58,
            'max_gap_size': 500
        }
        
        self.reference_analysis = {
            'total_length': 10000,
            'challenging_regions': [
                {'start': 1000, 'end': 1200, 'length': 200},
                {'start': 5000, 'end': 5100, 'length': 100}
            ]
        }
    
    def test_calculate_quality_score(self):
        """Test overall quality score calculation."""
        scorer = QualityScorer(
            self.coverage_stats,
            self.gap_analysis, 
            self.reference_analysis
        )
        
        quality_score = scorer.calculate_quality_score()
        
        # Check quality score structure
        self.assertIsInstance(quality_score, QualityScore)
        self.assertTrue(0 <= quality_score.overall_score <= 10)
        self.assertIn(quality_score.category, [QualityCategory.EXCELLENT, QualityCategory.GOOD, QualityCategory.FAIR, QualityCategory.POOR])
        
        # Check component scores exist
        self.assertIsInstance(quality_score.component_scores, dict)
        self.assertIn('coverage_breadth', quality_score.component_scores)
    
    def test_component_scoring(self):
        """Test individual component scoring."""
        scorer = QualityScorer(
            self.coverage_stats,
            self.gap_analysis,
            self.reference_analysis
        )
        
        # Test coverage scoring
        coverage_score = scorer._score_coverage_breadth()
        self.assertTrue(0 <= coverage_score <= 1)
        
        # Test depth uniformity scoring
        depth_score = scorer._score_coverage_depth()
        self.assertTrue(0 <= depth_score <= 1)
        
        # Test gap scoring
        gap_score = scorer._score_gap_characteristics()
        self.assertTrue(0 <= gap_score <= 1)


class TestBenchmarkAnalyzer(unittest.TestCase):
    """Test benchmarking system."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.coverage_stats = {
            'coverage_breadth': 78.5,
            'mean_depth': 8.2,
            'gini_coefficient': 0.4,
            'mapping_efficiency': 85.3,
            'reference_length': 15000
        }
        
        self.gap_analysis = {
            'total_gaps': 45,
            'total_gap_length': 3225,
            'mean_gap_size': 71.7,
            'max_gap_size': 800
        }
        
        self.reference_analysis = {
            'total_length': 15000,
            'challenging_regions': [
                {'start': 2000, 'end': 2300, 'length': 300},
                {'start': 8000, 'end': 8150, 'length': 150}
            ]
        }
        
        self.quality_score = create_quality_score(7.2, 'B')
    
    def test_theoretical_optimal_calculation(self):
        """Test theoretical optimal metrics calculation."""
        analyzer = BenchmarkAnalyzer(
            self.coverage_stats,
            self.gap_analysis,
            self.reference_analysis,
            self.quality_score
        )
        
        theoretical = analyzer.calculate_theoretical_optimal()
        
        # Check theoretical optimal structure
        self.assertIsInstance(theoretical, TheoreticalOptimal)
        self.assertTrue(0 <= theoretical.max_coverage_breadth <= 100)
        self.assertTrue(0 <= theoretical.theoretical_quality_score <= 10)
        self.assertGreaterEqual(theoretical.min_possible_gaps, 0)
    
    def test_benchmark_coverage_breadth(self):
        """Test coverage breadth benchmarking."""
        analyzer = BenchmarkAnalyzer(
            self.coverage_stats,
            self.gap_analysis,
            self.reference_analysis,
            self.quality_score
        )
        
        theoretical = analyzer.calculate_theoretical_optimal()
        result = analyzer.benchmark_coverage_breadth(theoretical)
        
        # Check benchmark result structure
        self.assertIsInstance(result, BenchmarkResult)
        self.assertTrue(0 <= result.efficiency_ratio <= 1)
        self.assertIn(result.category, ['Excellent', 'Good', 'Fair', 'Poor'])
        self.assertIsInstance(result.recommendations, list)
    
    def test_full_benchmark_analysis(self):
        """Test complete benchmark analysis."""
        analyzer = BenchmarkAnalyzer(
            self.coverage_stats,
            self.gap_analysis,
            self.reference_analysis,
            self.quality_score
        )
        
        benchmarks = analyzer.run_full_benchmark()
        
        # Check benchmark results structure
        expected_metrics = ['coverage_breadth', 'depth_uniformity', 
                           'gap_reduction', 'overall_quality']
        for metric in expected_metrics:
            self.assertIn(metric, benchmarks)
            self.assertIsInstance(benchmarks[metric], BenchmarkResult)
    
    def test_benchmark_report_generation(self):
        """Test benchmark report generation."""
        analyzer = BenchmarkAnalyzer(
            self.coverage_stats,
            self.gap_analysis,
            self.reference_analysis,
            self.quality_score
        )
        
        benchmarks = analyzer.run_full_benchmark()
        report = analyzer.generate_benchmark_report(benchmarks)
        
        # Check report content
        self.assertIsInstance(report, str)
        self.assertIn("PERFORMANCE BENCHMARK ANALYSIS", report)
        self.assertIn("Improvement Recommendations", report)


class TestInteractiveReportGenerator(unittest.TestCase):
    """Test interactive HTML report generation."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.output_dir = self.test_dir / "output"
        self.output_dir.mkdir()
        
        # Mock data
        self.coverage_stats = {'coverage_breadth': 85.0, 'mean_depth': 10.5}
        self.gap_analysis = {'total_gaps': 20}
        self.reference_analysis = {'total_length': 10000}
        self.quality_score = create_quality_score(8.5, 'A')
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    @patch('builtins.open', new_callable=mock_open)
    def test_report_generation(self, mock_file):
        """Test HTML report generation."""
        quality_scores = {
            'overall_score': self.quality_score.overall_score,
            'component_scores': self.quality_score.component_scores
        }
        generator = InteractiveReportGenerator(
            self.coverage_stats,
            self.gap_analysis,
            quality_scores,
            self.reference_analysis,
            self.test_dir
        )
        
        # Test report generation doesn't crash
        try:
            generator.generate_report()
        except Exception as e:
            self.fail(f"Report generation failed: {e}")
    
    def test_html_structure_validation(self):
        """Test HTML report structure."""
        quality_scores = {
            'overall_score': self.quality_score.overall_score,
            'component_scores': self.quality_score.component_scores
        }
        generator = InteractiveReportGenerator(
            self.coverage_stats,
            self.gap_analysis,
            quality_scores,
            self.reference_analysis,
            self.test_dir
        )
        
        html_content = generator._assemble_html()
        
        # Check for essential HTML elements
        self.assertIn('<html', html_content)  # Check for html tag (with or without attributes)
        self.assertIn('<head>', html_content)
        self.assertIn('<body>', html_content)
        self.assertIn('</html>', html_content)
        
        # The _assemble_html method only returns the template structure
        # Content is populated by the generate_report method
        # Just check that we have a valid HTML structure
        self.assertIn('container', html_content)
        self.assertIn('baitUtils Coverage Evaluation Report', html_content)


class TestInteractivePlotter(unittest.TestCase):
    """Test interactive plotting functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.output_dir = self.test_dir / "plots"
        self.output_dir.mkdir()
        
        # Mock data
        self.coverage_stats = {
            'coverage_breadth': 82.3,
            'mean_depth': 9.7,
            'depth_distribution': {1: 82.3, 5: 65.2, 10: 45.1}
        }
        self.gap_analysis = {
            'gaps': [
                {'chromosome': 'chr1', 'start': 1000, 'end': 1200, 'length': 200},
                {'chromosome': 'chr1', 'start': 3000, 'end': 3150, 'length': 150}
            ]
        }
        self.reference_analysis = {'total_length': 10000}
        self.quality_score = create_quality_score(7.8, 'B')
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    @patch('plotly.graph_objects.Figure')
    def test_plot_generation(self, mock_figure):
        """Test interactive plot generation."""
        # Mock plotly figure
        mock_fig = MagicMock()
        mock_figure.return_value = mock_fig
        
        quality_scores = {
            'overall_score': self.quality_score.overall_score,
            'component_scores': self.quality_score.component_scores
        }
        plotter = InteractivePlotter(
            self.coverage_stats,
            self.gap_analysis,
            quality_scores,
            self.reference_analysis,
            self.test_dir
        )
        
        # Test individual plot generation
        try:
            plotter.create_coverage_dashboard()
            plotter.create_reference_analysis_plot()
            plotter.create_quality_assessment_plot()
        except Exception as e:
            self.fail(f"Plot generation failed: {e}")
    
    @patch('plotly.graph_objects.Figure.write_html')
    def test_plot_saving(self, mock_write_html):
        """Test plot file saving."""
        quality_scores = {
            'overall_score': self.quality_score.overall_score,
            'component_scores': self.quality_score.component_scores
        }
        plotter = InteractivePlotter(
            self.coverage_stats,
            self.gap_analysis,
            quality_scores,
            self.reference_analysis,
            self.test_dir
        )
        
        # Test that plots are saved
        plotter.create_all_interactive_plots()
        
        # Check that write_html was called
        self.assertTrue(mock_write_html.called)


class TestPhase2Integration(unittest.TestCase):
    """Test Phase 2 integration with main evaluate command."""
    
    @patch('baitUtils.evaluate.ReferenceAnalyzer')
    @patch('baitUtils.evaluate.QualityScorer')
    @patch('baitUtils.evaluate.BenchmarkAnalyzer')
    @patch('baitUtils.evaluate.InteractiveReportGenerator')
    @patch('baitUtils.evaluate.InteractivePlotter')
    def test_phase2_workflow_integration(self, mock_plotter, mock_report,
                                       mock_benchmark, mock_quality, mock_ref):
        """Test Phase 2 workflow integration."""
        # This would test the integration in the main evaluate.py
        # For now, just ensure the imports work correctly
        
        from baitUtils import evaluate
        
        # Check that all Phase 2 classes are importable
        self.assertTrue(hasattr(evaluate, 'ReferenceAnalyzer'))
        self.assertTrue(hasattr(evaluate, 'QualityScorer'))
        self.assertTrue(hasattr(evaluate, 'BenchmarkAnalyzer'))
        self.assertTrue(hasattr(evaluate, 'InteractiveReportGenerator'))
        self.assertTrue(hasattr(evaluate, 'InteractivePlotter'))


if __name__ == '__main__':
    unittest.main()