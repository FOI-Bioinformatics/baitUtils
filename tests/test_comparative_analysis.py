#!/usr/bin/env python3

"""
test_comparative_analysis.py

Unit tests for Phase 3 comparative analysis functionality.
Tests comparative analyzer, differential analysis, visualizations, and reporting.
"""

import unittest
import tempfile
import shutil
import sys
import os
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
import pandas as pd
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from baitUtils.comparative_analyzer import ComparativeAnalyzer, OligoSetResult, ComparisonMetrics
from baitUtils.differential_analysis import DifferentialAnalyzer, StatisticalTest, CoverageDistributionComparison
from baitUtils.comparative_visualizations import ComparativeVisualizer
from baitUtils.comparative_report_generator import ComparativeReportGenerator
from baitUtils.quality_scorer import QualityScore, QualityCategory


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
            'coverage_score': score,
            'depth_score': score,
            'gap_score': score,
            'mapping_score': score
        },
        weighted_scores={
            'coverage_score': score,
            'depth_score': score,
            'gap_score': score,
            'mapping_score': score
        },
        benchmarks={
            'excellent_threshold': 9.0,
            'good_threshold': 7.0,
            'fair_threshold': 5.0
        },
        recommendations=[]
    )


class TestComparativeAnalyzer(unittest.TestCase):
    """Test ComparativeAnalyzer class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.reference_file = self.test_dir / "reference.fasta"
        self.output_dir = self.test_dir / "output"
        
        # Create mock reference file
        with open(self.reference_file, 'w') as f:
            f.write(">chr1\n")
            f.write("ATCGATCGATCG" * 100 + "\n")
        
        # Create mock oligo files
        self.oligo_files = {}
        for i, name in enumerate(['set1', 'set2', 'set3']):
            oligo_file = self.test_dir / f"{name}.fasta"
            with open(oligo_file, 'w') as f:
                f.write(f">{name}_oligo_1\n")
                f.write("ATCGATCGATCG" * 10 + "\n")
                f.write(f">{name}_oligo_2\n")
                f.write("GCTAGCTAGCTA" * 10 + "\n")
            self.oligo_files[name] = str(oligo_file)
        
        self.analyzer = ComparativeAnalyzer(
            reference_file=str(self.reference_file),
            output_dir=self.output_dir
        )
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    def test_init(self):
        """Test analyzer initialization."""
        self.assertEqual(str(self.analyzer.reference_file), str(self.reference_file))
        self.assertEqual(self.analyzer.output_dir, self.output_dir)
        self.assertEqual(self.analyzer.min_identity, 90.0)
        self.assertEqual(self.analyzer.min_length, 100)
        self.assertEqual(len(self.analyzer.oligo_sets), 0)
    
    @patch('baitUtils.comparative_analyzer.CoverageAnalyzer')
    @patch('baitUtils.comparative_analyzer.GapAnalyzer')  
    @patch('baitUtils.comparative_analyzer.ReferenceAnalyzer')
    @patch('baitUtils.comparative_analyzer.QualityScorer')
    @patch('baitUtils.comparative_analyzer.BenchmarkAnalyzer')
    @patch('subprocess.run')
    def test_add_oligo_set(self, mock_subprocess, mock_benchmark, mock_quality, 
                          mock_ref, mock_gap, mock_coverage):
        """Test adding oligo sets for comparison."""
        # Mock subprocess for pblat
        mock_subprocess.return_value = MagicMock(returncode=0)
        
        # Mock analysis results
        mock_coverage.return_value.analyze.return_value = {
            'coverage_breadth': 85.0,
            'mean_depth': 10.5,
            'mapping_efficiency': 90.0
        }
        
        mock_gap.return_value.analyze.return_value = {
            'total_gaps': 25,
            'mean_gap_size': 150
        }
        
        mock_ref.return_value.analyze.return_value = {
            'total_length': 10000,
            'challenging_regions': []
        }
        
        mock_quality.return_value.calculate_score.return_value = create_quality_score(8.5, 'A')
        
        mock_benchmark.return_value.run_full_benchmark.return_value = {}
        
        # Add oligo set
        self.analyzer.add_oligo_set('test_set', self.oligo_files['set1'])
        
        # Check that set was added
        self.assertEqual(len(self.analyzer.oligo_sets), 1)
        self.assertEqual(self.analyzer.oligo_sets[0].name, 'test_set')
    
    def test_generate_comparison_matrix(self):
        """Test comparison matrix generation."""
        # Add mock oligo sets
        mock_results = [
            OligoSetResult(
                name='Set1',
                file_path='set1.fasta',
                coverage_stats={'coverage_breadth': 80.0, 'mean_depth': 8.0, 'mapping_efficiency': 85.0},
                gap_analysis={'total_gaps': 30, 'max_gap_size': 500},
                quality_score=create_quality_score(7.5, 'B'),
                benchmark_results={}
            ),
            OligoSetResult(
                name='Set2', 
                file_path='set2.fasta',
                coverage_stats={'coverage_breadth': 90.0, 'mean_depth': 12.0, 'mapping_efficiency': 92.0},
                gap_analysis={'total_gaps': 15, 'max_gap_size': 300},
                quality_score=create_quality_score(8.8, 'A'),
                benchmark_results={}
            )
        ]
        
        self.analyzer.oligo_sets = mock_results
        
        # Generate comparison matrix
        matrix = self.analyzer.generate_comparison_matrix()
        
        # Check matrix structure
        self.assertIsInstance(matrix, pd.DataFrame)
        self.assertEqual(len(matrix), 2)
        self.assertIn('Name', matrix.columns)
        self.assertIn('Quality_Score', matrix.columns)
        self.assertIn('Coverage_Breadth_%', matrix.columns)
    
    def test_calculate_pairwise_comparisons(self):
        """Test pairwise comparison calculation."""
        # Add mock oligo sets
        mock_results = [
            OligoSetResult(
                name='Set1',
                file_path='set1.fasta',
                coverage_stats={'coverage_breadth': 80.0, 'mean_depth': 8.0},
                gap_analysis={'total_gaps': 30},
                quality_score=create_quality_score(7.5, 'B')
            ),
            OligoSetResult(
                name='Set2',
                file_path='set2.fasta', 
                coverage_stats={'coverage_breadth': 90.0, 'mean_depth': 12.0},
                gap_analysis={'total_gaps': 15},
                quality_score=create_quality_score(8.8, 'A')
            )
        ]
        
        self.analyzer.oligo_sets = mock_results
        
        # Calculate pairwise comparisons
        comparisons = self.analyzer.calculate_pairwise_comparisons()
        
        # Check comparisons
        self.assertEqual(len(comparisons), 1)  # Only one pair
        comparison_key = ('Set1', 'Set2')
        self.assertIn(comparison_key, comparisons)
        
        metrics = comparisons[comparison_key]
        self.assertIsInstance(metrics, ComparisonMetrics)
        self.assertEqual(metrics.coverage_breadth_diff, -10.0)  # Set1 - Set2
        self.assertEqual(metrics.gap_count_diff, 15)  # Set1 - Set2
    
    def test_identify_best_performer(self):
        """Test best performer identification."""
        # Add mock oligo sets
        mock_results = [
            OligoSetResult(
                name='Set1',
                file_path='set1.fasta',
                coverage_stats={'coverage_breadth': 80.0},
                gap_analysis={'total_gaps': 30},
                quality_score=create_quality_score(7.5, 'B')
            ),
            OligoSetResult(
                name='Set2',
                file_path='set2.fasta',
                coverage_stats={'coverage_breadth': 90.0},
                gap_analysis={'total_gaps': 15},
                quality_score=create_quality_score(8.8, 'A')
            )
        ]
        
        self.analyzer.oligo_sets = mock_results
        
        # Test different metrics
        best_quality = self.analyzer.identify_best_performer('quality_score')
        self.assertEqual(best_quality.name, 'Set2')
        
        best_coverage = self.analyzer.identify_best_performer('coverage_breadth')
        self.assertEqual(best_coverage.name, 'Set2')
        
        best_gaps = self.analyzer.identify_best_performer('gap_count')
        self.assertEqual(best_gaps.name, 'Set2')  # Fewer gaps is better
    
    def test_generate_ranking(self):
        """Test composite ranking generation."""
        # Add mock oligo sets
        mock_results = [
            OligoSetResult(
                name='Set1',
                file_path='set1.fasta',
                coverage_stats={'coverage_breadth': 80.0, 'mapping_efficiency': 85.0},
                gap_analysis={'total_gaps': 30},
                quality_score=create_quality_score(7.5, 'B')
            ),
            OligoSetResult(
                name='Set2',
                file_path='set2.fasta',
                coverage_stats={'coverage_breadth': 90.0, 'mapping_efficiency': 92.0},
                gap_analysis={'total_gaps': 15},
                quality_score=create_quality_score(8.8, 'A')
            )
        ]
        
        self.analyzer.oligo_sets = mock_results
        
        # Generate ranking
        ranking = self.analyzer.generate_ranking()
        
        # Check ranking
        self.assertEqual(len(ranking), 2)
        self.assertEqual(ranking[0][0], 'Set2')  # Set2 should be first
        self.assertEqual(ranking[1][0], 'Set1')  # Set1 should be second
        self.assertGreater(ranking[0][1], ranking[1][1])  # First should have higher score


class TestDifferentialAnalyzer(unittest.TestCase):
    """Test DifferentialAnalyzer class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = DifferentialAnalyzer(significance_level=0.05)
        
        # Mock oligo set results
        self.mock_set1 = OligoSetResult(
            name='Set1',
            file_path='set1.fasta',
            coverage_stats={
                'coverage_breadth': 80.0,
                'mean_depth': 8.0,
                'gini_coefficient': 0.3,
                'mapping_efficiency': 85.0,
                'reference_length': 10000
            },
            gap_analysis={'total_gaps': 30, 'gaps': []},
            quality_score=create_quality_score(7.5, 'B')
        )
        
        self.mock_set2 = OligoSetResult(
            name='Set2',
            file_path='set2.fasta',
            coverage_stats={
                'coverage_breadth': 90.0,
                'mean_depth': 12.0,
                'gini_coefficient': 0.2,
                'mapping_efficiency': 92.0,
                'reference_length': 10000
            },
            gap_analysis={'total_gaps': 15, 'gaps': []},
            quality_score=create_quality_score(8.8, 'A')
        )
    
    def test_init(self):
        """Test analyzer initialization."""
        self.assertEqual(self.analyzer.significance_level, 0.05)
        self.assertIn(0.05, self.analyzer.alpha_levels)
    
    @patch('scipy.stats.kstest')
    @patch('scipy.stats.mannwhitneyu')
    @patch('scipy.stats.levene')
    def test_compare_coverage_distributions(self, mock_levene, mock_mannwhitney, mock_kstest):
        """Test coverage distribution comparison."""
        # Mock statistical test results
        mock_kstest.return_value = (0.2, 0.01)
        mock_mannwhitney.return_value = (150, 0.03)
        mock_levene.return_value = (2.5, 0.12)
        
        comparison = self.analyzer.compare_coverage_distributions(self.mock_set1, self.mock_set2)
        
        # Check comparison structure
        self.assertIsInstance(comparison, CoverageDistributionComparison)
        self.assertEqual(comparison.set1_name, 'Set1')
        self.assertEqual(comparison.set2_name, 'Set2')
        
        # Check test results
        self.assertIsInstance(comparison.ks_test, StatisticalTest)
        self.assertIsInstance(comparison.mann_whitney_test, StatisticalTest)
        self.assertIsInstance(comparison.levene_test, StatisticalTest)
        
        # Check statistical significance
        self.assertEqual(comparison.ks_test.significance_level, '**')
        self.assertEqual(comparison.mann_whitney_test.significance_level, '*')
        self.assertEqual(comparison.levene_test.significance_level, 'ns')
    
    def test_compare_quality_metrics(self):
        """Test quality metrics comparison."""
        oligo_sets = [self.mock_set1, self.mock_set2]
        
        # Compare metrics
        results = self.analyzer.compare_quality_metrics(oligo_sets)
        
        # Check results structure
        expected_metrics = ['quality_score', 'coverage_breadth', 'mean_depth', 
                           'mapping_efficiency', 'gap_count', 'gini_coefficient']
        
        for metric in expected_metrics:
            self.assertIn(metric, results)
            self.assertIsInstance(results[metric], StatisticalTest)
    
    def test_multiple_comparison_correction(self):
        """Test multiple comparison corrections."""
        p_values = [0.01, 0.03, 0.08, 0.15, 0.25]
        
        # Test Bonferroni correction
        bonferroni = self.analyzer.multiple_comparison_correction(p_values, 'bonferroni')
        self.assertEqual(len(bonferroni), len(p_values))
        self.assertTrue(all(corr >= orig for corr, orig in zip(bonferroni, p_values)))
        
        # Test FDR correction
        fdr = self.analyzer.multiple_comparison_correction(p_values, 'fdr')
        self.assertEqual(len(fdr), len(p_values))
        
        # Test Holm correction
        holm = self.analyzer.multiple_comparison_correction(p_values, 'holm')
        self.assertEqual(len(holm), len(p_values))
    
    def test_significance_level_determination(self):
        """Test significance level assignment."""
        test_cases = [
            (0.0005, '***'),
            (0.005, '**'),
            (0.03, '*'),
            (0.1, 'ns')
        ]
        
        for p_value, expected in test_cases:
            result = self.analyzer._get_significance_level(p_value)
            self.assertEqual(result, expected)


class TestComparativeVisualizer(unittest.TestCase):
    """Test ComparativeVisualizer class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.visualizer = ComparativeVisualizer(
            output_dir=self.test_dir,
            plot_format='png'
        )
        
        # Mock analyzer with oligo sets
        self.mock_analyzer = MagicMock()
        self.mock_analyzer.oligo_sets = [
            OligoSetResult(
                name='Set1',
                file_path='set1.fasta',
                coverage_stats={'coverage_breadth': 80.0, 'mean_depth': 8.0, 'mapping_efficiency': 85.0},
                gap_analysis={'total_gaps': 30, 'max_gap_size': 500, 'gap_percentage': 20.0},
                quality_score=create_quality_score(7.5, 'B')
            ),
            OligoSetResult(
                name='Set2',
                file_path='set2.fasta',
                coverage_stats={'coverage_breadth': 90.0, 'mean_depth': 12.0, 'mapping_efficiency': 92.0},
                gap_analysis={'total_gaps': 15, 'max_gap_size': 300, 'gap_percentage': 10.0},
                quality_score=create_quality_score(8.8, 'A')
            )
        ]
        
        # Mock comparison matrix
        self.mock_analyzer.generate_comparison_matrix.return_value = pd.DataFrame({
            'Name': ['Set1', 'Set2'],
            'Quality_Score': [7.5, 8.8],
            'Quality_Grade': ['B', 'A'],
            'Coverage_Breadth_%': [80.0, 90.0],
            'Mean_Depth_x': [8.0, 12.0],
            'Total_Gaps': [30, 15],
            'Mapping_Efficiency_%': [85.0, 92.0],
            'Gini_Coefficient': [0.3, 0.2]
        })
        
        self.mock_analyzer.generate_ranking.return_value = [('Set2', 8.5), ('Set1', 7.2)]
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    def test_init(self):
        """Test visualizer initialization."""
        self.assertEqual(self.visualizer.output_dir, self.test_dir)
        self.assertEqual(self.visualizer.plot_format, 'png')
        self.assertTrue(self.visualizer.plots_dir.exists())
        self.assertTrue(self.visualizer.interactive_dir.exists())
    
    @patch('plotly.graph_objects.Figure')
    def test_create_comparison_dashboard(self, mock_figure):
        """Test comparison dashboard creation."""
        mock_fig = MagicMock()
        mock_figure.return_value = mock_fig
        
        result = self.visualizer.create_comparison_dashboard(self.mock_analyzer)
        
        # Check that figure was created and saved
        mock_figure.assert_called()
        mock_fig.write_html.assert_called_once()
        
        # Check return value
        self.assertIsInstance(result, str)
        self.assertTrue(result.endswith('.html'))
    
    @patch('matplotlib.pyplot.savefig')
    @patch('matplotlib.pyplot.subplots')
    def test_create_coverage_distribution_plots(self, mock_subplots, mock_savefig):
        """Test coverage distribution plots creation."""
        # Mock matplotlib
        mock_fig = MagicMock()
        mock_axes = np.array([[MagicMock(), MagicMock()], [MagicMock(), MagicMock()]])
        mock_subplots.return_value = (mock_fig, mock_axes)
        
        result = self.visualizer.create_coverage_distribution_plots(self.mock_analyzer.oligo_sets)
        
        # Check that plot was created
        mock_subplots.assert_called_once()
        mock_savefig.assert_called_once()
        
        # Check return value
        self.assertIsInstance(result, str)
    
    @patch('matplotlib.pyplot.savefig')
    @patch('matplotlib.pyplot.subplots')
    def test_create_gap_analysis_comparison(self, mock_subplots, mock_savefig):
        """Test gap analysis comparison plots."""
        mock_fig = MagicMock()
        mock_axes = np.array([[MagicMock(), MagicMock()], [MagicMock(), MagicMock()]])
        mock_subplots.return_value = (mock_fig, mock_axes)
        
        result = self.visualizer.create_gap_analysis_comparison(self.mock_analyzer.oligo_sets)
        
        # Check that plot was created
        mock_subplots.assert_called_once()
        mock_savefig.assert_called_once()
        
        # Check return value
        self.assertIsInstance(result, str)
    
    @patch('matplotlib.pyplot.savefig')
    @patch('seaborn.heatmap')
    def test_create_performance_heatmap(self, mock_heatmap, mock_savefig):
        """Test performance heatmap creation."""
        result = self.visualizer.create_performance_heatmap(self.mock_analyzer)
        
        # Check that heatmap was created
        mock_heatmap.assert_called_once()
        mock_savefig.assert_called_once()
        
        # Check return value
        self.assertIsInstance(result, str)
    
    @patch('plotly.graph_objects.Figure')
    def test_create_interactive_ranking_plot(self, mock_figure):
        """Test interactive ranking plot creation."""
        mock_fig = MagicMock()
        mock_figure.return_value = mock_fig
        
        result = self.visualizer.create_interactive_ranking_plot(self.mock_analyzer)
        
        # Check that figure was created and saved
        mock_figure.assert_called()
        mock_fig.write_html.assert_called_once()
        
        # Check return value
        self.assertIsInstance(result, str)
        self.assertTrue(result.endswith('.html'))


class TestComparativeReportGenerator(unittest.TestCase):
    """Test ComparativeReportGenerator class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        
        # Mock analyzer
        self.mock_analyzer = MagicMock()
        self.mock_analyzer.oligo_sets = [
            OligoSetResult(
                name='Set1',
                file_path='set1.fasta',
                coverage_stats={'coverage_breadth': 80.0, 'mean_depth': 8.0},
                gap_analysis={'total_gaps': 30},
                quality_score=create_quality_score(7.5, 'B')
            ),
            OligoSetResult(
                name='Set2',
                file_path='set2.fasta',
                coverage_stats={'coverage_breadth': 90.0, 'mean_depth': 12.0},
                gap_analysis={'total_gaps': 15},
                quality_score=create_quality_score(8.8, 'A')
            )
        ]
        
        self.mock_analyzer.generate_comparison_matrix.return_value = pd.DataFrame({
            'Name': ['Set1', 'Set2'],
            'Quality_Score': [7.5, 8.8],
            'Quality_Grade': ['B', 'A']
        })
        
        self.mock_analyzer.generate_ranking.return_value = [('Set2', 8.5), ('Set1', 7.2)]
        self.mock_analyzer.identify_best_performer.return_value = self.mock_analyzer.oligo_sets[1]
        self.mock_analyzer.analyze_gap_overlap.return_value = {}
        
        self.generator = ComparativeReportGenerator(
            analyzer=self.mock_analyzer,
            output_dir=self.test_dir
        )
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    def test_init(self):
        """Test report generator initialization."""
        self.assertEqual(self.generator.analyzer, self.mock_analyzer)
        self.assertEqual(self.generator.output_dir, self.test_dir)
        self.assertTrue(self.generator.report_file.parent.exists())
    
    @patch.object(ComparativeVisualizer, 'generate_all_comparative_plots')
    def test_generate_report(self, mock_generate_plots):
        """Test report generation."""
        mock_generate_plots.return_value = {
            'dashboard': str(self.test_dir / 'dashboard.html'),
            'coverage_distributions': str(self.test_dir / 'coverage.png')
        }
        
        result = self.generator.generate_report()
        
        # Check that report file was created
        self.assertTrue(Path(result).exists())
        self.assertTrue(result.endswith('.html'))
        
        # Check that plots were generated
        mock_generate_plots.assert_called_once()
    
    def test_html_content_structure(self):
        """Test HTML content structure."""
        html_header = self.generator._get_html_header()
        html_footer = self.generator._get_html_footer()
        
        # Check header contains essential elements
        self.assertIn('<!DOCTYPE html>', html_header)
        self.assertIn('<head>', html_header)
        self.assertIn('<title>', html_header)
        self.assertIn('<style>', html_header)
        
        # Check footer contains closing tags
        self.assertIn('</body>', html_footer)
        self.assertIn('</html>', html_footer)
    
    def test_executive_summary_generation(self):
        """Test executive summary section."""
        summary = self.generator._generate_executive_summary()
        
        # Check content
        self.assertIn('Executive Summary', summary)
        self.assertIn('Set2', summary)  # Best performer
        self.assertIn('8.8/10', summary)  # Quality score
        self.assertIn('Grade A', summary)
    
    def test_comparison_overview_generation(self):
        """Test comparison overview table."""
        comparison_matrix = self.mock_analyzer.generate_comparison_matrix.return_value
        overview = self.generator._generate_comparison_overview(comparison_matrix)
        
        # Check content
        self.assertIn('Comparison Overview', overview)
        self.assertIn('<table', overview)
        self.assertIn('Set1', overview)
        self.assertIn('Set2', overview)


if __name__ == '__main__':
    unittest.main()