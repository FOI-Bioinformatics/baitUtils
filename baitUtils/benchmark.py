#!/usr/bin/env python3

"""
benchmark.py

Compares observed coverage metrics with the design targets used by the
quality scorer and reports how far each metric is from its target.

The targets are those in QualityScorer.DEFAULT_BENCHMARKS unless overridden.
No attempt is made to derive a "theoretical optimum" from the reference; the
comparison is against stated targets, which keeps the result interpretable.
"""

import logging
from typing import Dict, List, NamedTuple, Optional

from baitUtils.quality_scorer import QualityScorer, QualityScore


class BenchmarkResult(NamedTuple):
    """Comparison of one observed metric with its target."""
    actual_score: float
    target: float
    efficiency_ratio: float          # 0-1, fraction of the target achieved
    improvement_potential: float     # distance to target in the metric's unit
    category: str
    recommendations: List[str]


class BenchmarkAnalyzer:
    """Benchmark observed coverage metrics against design targets."""

    def __init__(
        self,
        coverage_stats: Dict,
        gap_analysis: Dict,
        reference_analysis: Optional[Dict] = None,
        quality_score: Optional[QualityScore] = None,
        targets: Optional[Dict[str, float]] = None,
    ):
        """
        Args:
            coverage_stats: Output of CoverageAnalyzer.analyze
            gap_analysis: Output of GapAnalyzer.analyze
            reference_analysis: Output of ReferenceAnalyzer.analyze (unused
                by the benchmarks but kept for report context)
            quality_score: Optional QualityScore for the report header
            targets: Overrides for QualityScorer.DEFAULT_BENCHMARKS
        """
        self.coverage_stats = coverage_stats
        self.gap_analysis = gap_analysis
        self.reference_analysis = reference_analysis or {}
        self.quality_score = quality_score
        self.targets = dict(QualityScorer.DEFAULT_BENCHMARKS)
        if targets:
            self.targets.update(targets)

    @staticmethod
    def _categorize(efficiency_ratio: float) -> str:
        thresholds = QualityScorer.DEFAULT_THRESHOLDS
        if efficiency_ratio >= thresholds['excellent']:
            return "Excellent"
        if efficiency_ratio >= thresholds['good']:
            return "Good"
        if efficiency_ratio >= thresholds['fair']:
            return "Fair"
        return "Poor"

    def benchmark_coverage_breadth(self) -> BenchmarkResult:
        """Observed coverage breadth (%) against the target breadth."""
        actual = float(self.coverage_stats.get('coverage_breadth', 0.0))
        target = float(self.targets['coverage_breadth'])
        ratio = min(1.0, actual / target) if target > 0 else 1.0
        potential = max(0.0, target - actual)

        recommendations = []
        if ratio < 0.8:
            recommendations.append(
                f"Coverage breadth is {actual:.1f}% against a target of {target:.1f}%; "
                "add oligos in uncovered regions or relax mapping thresholds"
            )
        return BenchmarkResult(actual, target, ratio, potential, self._categorize(ratio), recommendations)

    def benchmark_depth_uniformity(self) -> BenchmarkResult:
        """Observed uniformity (1 - Gini) against the target uniformity."""
        gini = float(self.coverage_stats.get('coverage_gini', 0.5))
        actual = 1.0 - gini
        target = float(self.targets['depth_uniformity'])
        ratio = min(1.0, actual / target) if target > 0 else 1.0
        potential = max(0.0, target - actual)

        recommendations = []
        if ratio < 0.7:
            recommendations.append(
                f"Depth uniformity is {actual:.2f} (1 - Gini) against a target of {target:.2f}; "
                "consider redistributing oligo density"
            )
        return BenchmarkResult(actual, target, ratio, potential, self._categorize(ratio), recommendations)

    def benchmark_gap_reduction(self) -> BenchmarkResult:
        """Observed gap percentage against the target gap percentage."""
        actual = float(self.gap_analysis.get('gap_percentage', 0.0))
        target = float(self.targets['gap_percentage'])
        if actual <= target:
            ratio = 1.0
        else:
            ratio = target / actual if actual > 0 else 1.0
        potential = max(0.0, actual - target)

        recommendations = []
        if ratio < 0.6:
            recommendations.append(
                f"Gaps cover {actual:.1f}% of the reference against a target of {target:.1f}%; "
                "use 'baitUtils fill' to close gaps with additional oligos"
            )
        return BenchmarkResult(actual, target, ratio, potential, self._categorize(ratio), recommendations)

    def run_full_benchmark(self) -> Dict[str, BenchmarkResult]:
        """Run all benchmarks."""
        logging.info("Benchmarking coverage metrics against design targets...")
        return {
            'coverage_breadth': self.benchmark_coverage_breadth(),
            'depth_uniformity': self.benchmark_depth_uniformity(),
            'gap_reduction': self.benchmark_gap_reduction(),
        }

    def generate_benchmark_report(self, benchmarks: Dict[str, BenchmarkResult]) -> str:
        """Format benchmark results as plain text."""
        lines = ["=" * 80, "BENCHMARK AGAINST DESIGN TARGETS", "=" * 80, ""]
        if self.quality_score is not None:
            lines.append(f"Overall quality score: {self.quality_score.overall_score:.2f} "
                         f"({self.quality_score.category.value})")
            lines.append("")

        for metric, result in benchmarks.items():
            lines.append(f"{metric.replace('_', ' ').title()}:")
            lines.append(f"  Observed:    {result.actual_score:8.2f}")
            lines.append(f"  Target:      {result.target:8.2f}")
            lines.append(f"  Achieved:    {result.efficiency_ratio * 100:8.1f}%")
            lines.append(f"  Category:    {result.category:>8s}")
            lines.append("")

        lines.append("Recommendations:")
        lines.append("-" * 40)
        recommendations = [r for result in benchmarks.values() for r in result.recommendations]
        if recommendations:
            lines.extend(f"{i}. {rec}" for i, rec in enumerate(recommendations, 1))
        else:
            lines.append("All metrics are at or near their targets.")
        lines.extend(["", "=" * 80])
        return "\n".join(lines)
