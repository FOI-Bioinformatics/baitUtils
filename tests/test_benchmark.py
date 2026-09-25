"""
Tests for BenchmarkAnalyzer: observed metrics against stated design targets.
"""

import pytest

from baitUtils.benchmark import BenchmarkAnalyzer, BenchmarkResult
from baitUtils.quality_scorer import QualityScorer, QualityScore, QualityCategory


@pytest.fixture
def coverage_stats():
    return {
        'coverage_breadth': 75.0,
        'mean_depth': 8.2,
        'coverage_gini': 0.35,
        'mapping_efficiency': 88.3,
        'reference_length': 12000,
    }


@pytest.fixture
def gap_analysis():
    return {'total_gaps': 35, 'total_gap_length': 2940, 'gap_percentage': 24.5, 'max_gap_size': 650}


@pytest.fixture
def quality_score():
    return QualityScore(
        overall_score=0.71, category=QualityCategory.GOOD, component_scores={},
        weighted_scores={}, benchmarks={}, recommendations=[])


class TestBenchmarkValues:
    def test_targets_default_to_quality_scorer(self, coverage_stats, gap_analysis):
        analyzer = BenchmarkAnalyzer(coverage_stats, gap_analysis)
        assert analyzer.targets == QualityScorer.DEFAULT_BENCHMARKS

    def test_coverage_breadth_ratio(self, coverage_stats, gap_analysis):
        result = BenchmarkAnalyzer(coverage_stats, gap_analysis).benchmark_coverage_breadth()
        assert result.actual_score == 75.0
        assert result.target == 95.0
        assert result.efficiency_ratio == pytest.approx(75.0 / 95.0)
        assert result.improvement_potential == pytest.approx(20.0)
        assert result.category == "Good"
        assert result.recommendations != []

    def test_depth_uniformity_uses_one_minus_gini(self, coverage_stats, gap_analysis):
        result = BenchmarkAnalyzer(coverage_stats, gap_analysis).benchmark_depth_uniformity()
        assert result.actual_score == pytest.approx(0.65)
        assert result.target == 0.8
        assert result.efficiency_ratio == pytest.approx(0.65 / 0.8)

    def test_gap_reduction_ratio_is_target_over_actual(self, coverage_stats, gap_analysis):
        result = BenchmarkAnalyzer(coverage_stats, gap_analysis).benchmark_gap_reduction()
        assert result.actual_score == 24.5
        assert result.efficiency_ratio == pytest.approx(2.0 / 24.5)
        assert result.improvement_potential == pytest.approx(22.5)
        assert result.category == "Poor"
        assert any("fill" in r for r in result.recommendations)

    def test_metrics_at_target_are_excellent(self, gap_analysis):
        stats = {'coverage_breadth': 100.0, 'coverage_gini': 0.0}
        gaps = {'gap_percentage': 0.0}
        benchmarks = BenchmarkAnalyzer(stats, gaps).run_full_benchmark()
        assert set(benchmarks) == {'coverage_breadth', 'depth_uniformity', 'gap_reduction'}
        for result in benchmarks.values():
            assert isinstance(result, BenchmarkResult)
            assert result.efficiency_ratio == 1.0
            assert result.improvement_potential == 0.0
            assert result.category == "Excellent"
            assert result.recommendations == []

    def test_target_override(self, coverage_stats, gap_analysis):
        analyzer = BenchmarkAnalyzer(coverage_stats, gap_analysis, targets={'coverage_breadth': 75.0})
        assert analyzer.benchmark_coverage_breadth().efficiency_ratio == 1.0


class TestBenchmarkReport:
    def test_report_lists_each_metric_and_recommendations(self, coverage_stats, gap_analysis, quality_score):
        analyzer = BenchmarkAnalyzer(coverage_stats, gap_analysis, quality_score=quality_score)
        report = analyzer.generate_benchmark_report(analyzer.run_full_benchmark())
        assert "Overall quality score: 0.71 (Good)" in report
        for heading in ("Coverage Breadth:", "Depth Uniformity:", "Gap Reduction:"):
            assert heading in report
        assert "Target:" in report
        assert "1. " in report

    def test_report_without_recommendations(self):
        analyzer = BenchmarkAnalyzer({'coverage_breadth': 100.0, 'coverage_gini': 0.0}, {'gap_percentage': 0.0})
        report = analyzer.generate_benchmark_report(analyzer.run_full_benchmark())
        assert "at or near their targets" in report
