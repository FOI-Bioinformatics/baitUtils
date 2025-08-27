#!/usr/bin/env python3

"""
quality_scorer.py

Comprehensive coverage quality scoring system for oligo set evaluation.
Provides standardized metrics and benchmarking against theoretical optimal coverage.
"""

import logging
from typing import Dict, List, Tuple, Any, Optional
import numpy as np
import pandas as pd
from dataclasses import dataclass
from enum import Enum


class QualityCategory(Enum):
    """Quality categories for coverage assessment."""
    EXCELLENT = "Excellent"
    GOOD = "Good"
    FAIR = "Fair"
    POOR = "Poor"


@dataclass
class QualityScore:
    """Quality score with breakdown and recommendations."""
    overall_score: float
    category: QualityCategory
    component_scores: Dict[str, float]
    weighted_scores: Dict[str, float]
    benchmarks: Dict[str, float]
    recommendations: List[str]


class QualityScorer:
    """Comprehensive quality scoring system for coverage evaluation."""
    
    def __init__(
        self,
        coverage_stats: Dict[str, Any],
        gap_analysis: Dict[str, Any],
        reference_analysis: Dict[str, Any] = None
    ):
        """
        Initialize the quality scorer.
        
        Args:
            coverage_stats: Coverage statistics from CoverageAnalyzer
            gap_analysis: Gap analysis results from GapAnalyzer
            reference_analysis: Reference sequence analysis (optional)
        """
        self.coverage_stats = coverage_stats
        self.gap_analysis = gap_analysis
        self.reference_analysis = reference_analysis or {}
        
        # Scoring weights (can be customized)
        self.weights = {
            'coverage_breadth': 0.30,      # Most important metric
            'coverage_depth': 0.20,        # Depth adequacy and uniformity
            'mapping_efficiency': 0.15,    # Oligo utilization
            'gap_characteristics': 0.20,   # Gap size and distribution
            'reference_difficulty': 0.15   # Reference sequence challenges
        }
        
        # Quality thresholds
        self.thresholds = {
            'excellent': 0.85,
            'good': 0.70,
            'fair': 0.50,
            'poor': 0.0
        }
        
        # Benchmark values (ideal targets)
        self.benchmarks = {
            'coverage_breadth': 95.0,      # 95% coverage
            'mean_depth': 20.0,            # 20x mean depth
            'depth_uniformity': 0.8,       # Low CV (high uniformity)
            'mapping_efficiency': 90.0,    # 90% oligos mapping
            'gap_percentage': 2.0,         # <2% gaps
            'largest_gap': 1000            # <1kb largest gap
        }
    
    def calculate_quality_score(self) -> QualityScore:
        """
        Calculate comprehensive quality score.
        
        Returns:
            QualityScore object with overall score and breakdown
        """
        logging.info("Calculating comprehensive quality score...")
        
        # Calculate component scores
        component_scores = {
            'coverage_breadth': self._score_coverage_breadth(),
            'coverage_depth': self._score_coverage_depth(),
            'mapping_efficiency': self._score_mapping_efficiency(),
            'gap_characteristics': self._score_gap_characteristics(),
            'reference_difficulty': self._score_reference_difficulty()
        }
        
        # Calculate weighted scores
        weighted_scores = {
            component: score * self.weights[component]
            for component, score in component_scores.items()
        }
        
        # Overall score
        overall_score = sum(weighted_scores.values())
        
        # Determine category
        category = self._determine_category(overall_score)
        
        # Generate recommendations
        recommendations = self._generate_quality_recommendations(component_scores)
        
        # Calculate benchmarks
        benchmarks = self._calculate_benchmarks()
        
        return QualityScore(
            overall_score=overall_score,
            category=category,
            component_scores=component_scores,
            weighted_scores=weighted_scores,
            benchmarks=benchmarks,
            recommendations=recommendations
        )
    
    def _score_coverage_breadth(self) -> float:
        """Score coverage breadth component."""
        breadth = self.coverage_stats.get('coverage_breadth', 0)
        target_breadth = self.benchmarks['coverage_breadth']
        
        # Score based on how close to target
        if breadth >= target_breadth:
            return 1.0
        elif breadth >= target_breadth * 0.8:  # 80% of target
            return 0.8 + (breadth - target_breadth * 0.8) / (target_breadth * 0.2) * 0.2
        elif breadth >= target_breadth * 0.5:  # 50% of target
            return 0.5 + (breadth - target_breadth * 0.5) / (target_breadth * 0.3) * 0.3
        else:
            return max(0.0, breadth / (target_breadth * 0.5) * 0.5)
    
    def _score_coverage_depth(self) -> float:
        """Score coverage depth adequacy and uniformity."""
        mean_depth = self.coverage_stats.get('mean_depth', 0)
        depth_cv = self.coverage_stats.get('coverage_cv', float('inf'))
        target_depth = self.benchmarks['mean_depth']
        target_uniformity = self.benchmarks['depth_uniformity']
        
        # Depth adequacy score
        if mean_depth >= target_depth:
            depth_score = 1.0
        elif mean_depth >= target_depth * 0.5:
            depth_score = 0.5 + (mean_depth - target_depth * 0.5) / (target_depth * 0.5) * 0.5
        else:
            depth_score = max(0.0, mean_depth / (target_depth * 0.5) * 0.5)
        
        # Uniformity score (lower CV is better)
        if depth_cv == float('inf') or depth_cv > 2.0:
            uniformity_score = 0.0
        elif depth_cv <= 0.5:  # Very uniform
            uniformity_score = 1.0
        else:
            uniformity_score = max(0.0, 1.0 - (depth_cv - 0.5) / 1.5)
        
        # Combined score (equal weights)
        return (depth_score + uniformity_score) / 2
    
    def _score_mapping_efficiency(self) -> float:
        """Score mapping efficiency component."""
        efficiency = self.coverage_stats.get('mapping_efficiency', 0)
        target_efficiency = self.benchmarks['mapping_efficiency']
        
        if efficiency >= target_efficiency:
            return 1.0
        elif efficiency >= target_efficiency * 0.7:
            return 0.7 + (efficiency - target_efficiency * 0.7) / (target_efficiency * 0.3) * 0.3
        else:
            return max(0.0, efficiency / (target_efficiency * 0.7) * 0.7)
    
    def _score_gap_characteristics(self) -> float:
        """Score gap characteristics component."""
        gap_percentage = self.gap_analysis.get('gap_percentage', 100)
        largest_gap = self.gap_analysis.get('max_gap_size', float('inf'))
        total_gaps = self.gap_analysis.get('total_gaps', float('inf'))
        
        target_gap_pct = self.benchmarks['gap_percentage']
        target_largest_gap = self.benchmarks['largest_gap']
        
        # Gap percentage score
        if gap_percentage <= target_gap_pct:
            gap_pct_score = 1.0
        elif gap_percentage <= target_gap_pct * 3:
            gap_pct_score = 1.0 - (gap_percentage - target_gap_pct) / (target_gap_pct * 2) * 0.5
        else:
            gap_pct_score = max(0.0, 0.5 - (gap_percentage - target_gap_pct * 3) / 50 * 0.5)
        
        # Largest gap score
        if largest_gap <= target_largest_gap:
            largest_gap_score = 1.0
        elif largest_gap <= target_largest_gap * 5:
            largest_gap_score = 1.0 - (largest_gap - target_largest_gap) / (target_largest_gap * 4) * 0.6
        else:
            largest_gap_score = max(0.0, 0.4 - (largest_gap - target_largest_gap * 5) / 50000 * 0.4)
        
        # Gap count score (fewer gaps is better)
        if total_gaps <= 10:
            gap_count_score = 1.0
        elif total_gaps <= 50:
            gap_count_score = 1.0 - (total_gaps - 10) / 40 * 0.4
        else:
            gap_count_score = max(0.0, 0.6 - (total_gaps - 50) / 200 * 0.6)
        
        # Combined score
        return (gap_pct_score * 0.5 + largest_gap_score * 0.3 + gap_count_score * 0.2)
    
    def _score_reference_difficulty(self) -> float:
        """Score based on reference sequence difficulty."""
        if not self.reference_analysis:
            return 0.7  # Neutral score if no reference analysis
        
        summary = self.reference_analysis.get('analysis_summary', {})
        challenging_regions = self.reference_analysis.get('challenging_regions', {})
        
        total_sequences = summary.get('total_sequences', 1)
        challenging_count = summary.get('challenging_sequences', 0)
        
        # Calculate difficulty factors
        extreme_gc_fraction = summary.get('sequences_with_extreme_gc', 0) / total_sequences
        high_repeat_fraction = summary.get('sequences_with_high_repeats', 0) / total_sequences
        challenging_fraction = challenging_count / total_sequences
        
        # Score based on reference difficulty (inverse relationship)
        difficulty_score = 1.0
        
        # Penalize extreme GC content
        difficulty_score -= extreme_gc_fraction * 0.3
        
        # Penalize high repeat content
        difficulty_score -= high_repeat_fraction * 0.2
        
        # Penalize challenging sequences
        difficulty_score -= challenging_fraction * 0.3
        
        # Adjust based on actual coverage performance relative to difficulty
        expected_difficulty = extreme_gc_fraction + high_repeat_fraction + challenging_fraction
        actual_breadth = self.coverage_stats.get('coverage_breadth', 0) / 100
        
        if expected_difficulty > 0:
            performance_ratio = actual_breadth / (1.0 - expected_difficulty * 0.5)
            if performance_ratio > 1.0:  # Better than expected
                difficulty_score = min(1.0, difficulty_score + (performance_ratio - 1.0) * 0.2)
        
        return max(0.0, min(1.0, difficulty_score))
    
    def _determine_category(self, overall_score: float) -> QualityCategory:
        """Determine quality category from overall score."""
        if overall_score >= self.thresholds['excellent']:
            return QualityCategory.EXCELLENT
        elif overall_score >= self.thresholds['good']:
            return QualityCategory.GOOD
        elif overall_score >= self.thresholds['fair']:
            return QualityCategory.FAIR
        else:
            return QualityCategory.POOR
    
    def _generate_quality_recommendations(self, component_scores: Dict[str, float]) -> List[str]:
        """Generate recommendations based on component scores."""
        recommendations = []
        
        # Coverage breadth recommendations
        if component_scores['coverage_breadth'] < 0.7:
            breadth = self.coverage_stats.get('coverage_breadth', 0)
            recommendations.append(
                f"Coverage breadth is {breadth:.1f}%. Consider adding more oligos or using "
                "'baitUtils fill' to target uncovered regions."
            )
        
        # Coverage depth recommendations
        if component_scores['coverage_depth'] < 0.7:
            mean_depth = self.coverage_stats.get('mean_depth', 0)
            cv = self.coverage_stats.get('coverage_cv', 0)
            
            if mean_depth < 10:
                recommendations.append(
                    f"Mean coverage depth is {mean_depth:.1f}x. Consider increasing oligo density "
                    "or optimizing oligo placement for better depth."
                )
            
            if cv > 1.0:
                recommendations.append(
                    f"Coverage depth is highly variable (CV={cv:.2f}). Consider redistributing "
                    "oligos for more uniform coverage."
                )
        
        # Mapping efficiency recommendations
        if component_scores['mapping_efficiency'] < 0.7:
            efficiency = self.coverage_stats.get('mapping_efficiency', 0)
            recommendations.append(
                f"Mapping efficiency is {efficiency:.1f}%. Review oligo design parameters "
                "or reference sequence quality."
            )
        
        # Gap characteristics recommendations
        if component_scores['gap_characteristics'] < 0.7:
            gap_pct = self.gap_analysis.get('gap_percentage', 0)
            largest_gap = self.gap_analysis.get('max_gap_size', 0)
            
            if gap_pct > 10:
                recommendations.append(
                    f"Gap percentage is {gap_pct:.1f}%. Focus on gap closure using iterative "
                    "oligo selection methods."
                )
            
            if largest_gap > 5000:
                recommendations.append(
                    f"Largest gap is {largest_gap:,} bp. Very large gaps may indicate "
                    "problematic reference regions requiring specialized approaches."
                )
        
        # Reference difficulty recommendations
        if component_scores['reference_difficulty'] < 0.6:
            recommendations.append(
                "Reference sequences contain challenging regions (extreme GC, repeats). "
                "Consider specialized oligo design strategies for difficult regions."
            )
        
        # Overall recommendations
        overall_score = sum(score * self.weights[component] 
                          for component, score in component_scores.items())
        
        if overall_score < 0.5:
            recommendations.append(
                "Overall coverage quality is poor. Consider comprehensive redesign of "
                "the oligo set with adjusted parameters."
            )
        elif overall_score < 0.7:
            recommendations.append(
                "Coverage quality can be improved. Focus on the lowest-scoring components "
                "for targeted improvements."
            )
        
        return recommendations
    
    def _calculate_benchmarks(self) -> Dict[str, float]:
        """Calculate how actual performance compares to benchmarks."""
        benchmarks = {}
        
        # Coverage breadth benchmark
        actual_breadth = self.coverage_stats.get('coverage_breadth', 0)
        benchmarks['coverage_breadth_ratio'] = actual_breadth / self.benchmarks['coverage_breadth']
        
        # Depth benchmark
        actual_depth = self.coverage_stats.get('mean_depth', 0)
        benchmarks['depth_ratio'] = actual_depth / self.benchmarks['mean_depth']
        
        # Mapping efficiency benchmark
        actual_efficiency = self.coverage_stats.get('mapping_efficiency', 0)
        benchmarks['efficiency_ratio'] = actual_efficiency / self.benchmarks['mapping_efficiency']
        
        # Gap benchmark (inverse - lower is better)
        actual_gap_pct = self.gap_analysis.get('gap_percentage', 100)
        benchmarks['gap_ratio'] = self.benchmarks['gap_percentage'] / max(actual_gap_pct, 0.1)
        
        return benchmarks
    
    def generate_quality_report(self, quality_score: QualityScore) -> str:
        """Generate a formatted quality report."""
        report_lines = []
        report_lines.append("=== COVERAGE QUALITY ASSESSMENT ===")
        report_lines.append("")
        
        # Overall score
        report_lines.append(f"Overall Quality Score: {quality_score.overall_score:.3f}")
        report_lines.append(f"Quality Category: {quality_score.category.value}")
        report_lines.append("")
        
        # Component breakdown
        report_lines.append("Component Scores:")
        for component, score in quality_score.component_scores.items():
            weighted_score = quality_score.weighted_scores[component]
            weight = self.weights[component]
            
            report_lines.append(
                f"  {component.replace('_', ' ').title()}: "
                f"{score:.3f} (weighted: {weighted_score:.3f}, weight: {weight:.2f})"
            )
        
        report_lines.append("")
        
        # Benchmarks
        report_lines.append("Performance vs Benchmarks:")
        for metric, ratio in quality_score.benchmarks.items():
            status = "✓" if ratio >= 0.8 else "✗" if ratio < 0.5 else "~"
            report_lines.append(f"  {status} {metric.replace('_', ' ').title()}: {ratio:.2f}x target")
        
        report_lines.append("")
        
        # Recommendations
        if quality_score.recommendations:
            report_lines.append("Recommendations:")
            for i, rec in enumerate(quality_score.recommendations, 1):
                report_lines.append(f"  {i}. {rec}")
        
        return "\n".join(report_lines)
    
    def compare_to_theoretical_optimal(self) -> Dict[str, Any]:
        """Compare current coverage to theoretical optimal."""
        # Calculate theoretical optimal based on oligo characteristics
        total_oligos = self.coverage_stats.get('total_oligos', 0)
        mapped_oligos = self.coverage_stats.get('mapped_oligos', 0)
        reference_length = self.coverage_stats.get('reference_length', 1)
        
        # Assume average oligo length of 120bp
        avg_oligo_length = 120
        total_oligo_bases = mapped_oligos * avg_oligo_length
        
        # Theoretical maximum coverage (if perfectly distributed)
        theoretical_max_coverage = min(100.0, (total_oligo_bases / reference_length) * 100)
        
        # Theoretical optimal metrics
        theoretical_optimal = {
            'max_possible_breadth': theoretical_max_coverage,
            'optimal_mean_depth': total_oligo_bases / reference_length,
            'optimal_uniformity': 0.2,  # CV of 0.2 is very good
            'optimal_gaps': max(1, int((100 - theoretical_max_coverage) / 100 * reference_length / 1000)),
            'efficiency_ceiling': 95.0  # Account for mapping challenges
        }
        
        # Current vs theoretical
        current_performance = {
            'breadth_efficiency': self.coverage_stats.get('coverage_breadth', 0) / theoretical_max_coverage,
            'depth_efficiency': min(1.0, self.coverage_stats.get('mean_depth', 0) / theoretical_optimal['optimal_mean_depth']),
            'uniformity_efficiency': max(0.0, 1.0 - self.coverage_stats.get('coverage_cv', 2.0) / 2.0),
            'gap_efficiency': max(0.0, 1.0 - self.gap_analysis.get('total_gaps', 1000) / 1000)
        }
        
        return {
            'theoretical_optimal': theoretical_optimal,
            'current_performance': current_performance,
            'overall_efficiency': np.mean(list(current_performance.values())),
            'improvement_potential': {
                metric: max(0.0, 1.0 - efficiency) 
                for metric, efficiency in current_performance.items()
            }
        }
    
    def export_quality_data(self) -> Dict[str, Any]:
        """Export comprehensive quality data for analysis."""
        quality_score = self.calculate_quality_score()
        theoretical_comparison = self.compare_to_theoretical_optimal()
        
        return {
            'quality_score': {
                'overall_score': quality_score.overall_score,
                'category': quality_score.category.value,
                'component_scores': quality_score.component_scores,
                'weighted_scores': quality_score.weighted_scores,
                'benchmarks': quality_score.benchmarks
            },
            'theoretical_comparison': theoretical_comparison,
            'recommendations': quality_score.recommendations,
            'scoring_weights': self.weights,
            'quality_thresholds': self.thresholds,
            'benchmark_targets': self.benchmarks
        }