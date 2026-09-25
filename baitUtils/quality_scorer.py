#!/usr/bin/env python3

"""
quality_scorer.py

Comprehensive coverage quality scoring system for oligo set evaluation.
Provides standardized metrics and benchmarking against theoretical optimal coverage.
"""

import logging
from typing import Dict, List, Any, Optional
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

    def to_dict(self) -> Dict[str, Any]:
        """Return a plain dictionary representation for reports and JSON output."""
        return {
            'overall_score': self.overall_score,
            'category': self.category.value,
            'component_scores': dict(self.component_scores),
            'weighted_scores': dict(self.weighted_scores),
            'benchmarks': dict(self.benchmarks),
            'recommendations': list(self.recommendations),
        }


class QualityScorer:
    """
    Quality scoring system for coverage evaluation.

    The overall score is a weighted sum of five component scores in 0-1.
    Weights, targets and category thresholds are stated here so that they can
    be cited and overridden. They are design conventions for 120 bp
    hybridization baits at moderate tiling density, not fitted parameters.

    Component            Weight  Rationale
    -------------------  ------  -----------------------------------------------
    coverage_breadth     0.30    Fraction of target covered is the primary goal
    coverage_depth       0.20    Adequate and even depth aids capture efficiency
    gap_characteristics  0.20    Few, short gaps are easier to close
    mapping_efficiency   0.15    Oligos that do not map are wasted synthesis
    reference_difficulty 0.15    Credit for hard references (repeats, extreme GC)

    Target               Value   Meaning
    -------------------  ------  -----------------------------------------------
    coverage_breadth     95 %    Breadth at which breadth scores 1.0
    mean_depth           20 x    Mean depth at which depth adequacy scores 1.0
    depth_uniformity     0.8     Target uniformity expressed as 1 - Gini
    mapping_efficiency   90 %    Mapped fraction at which efficiency scores 1.0
    gap_percentage       2 %     Gap fraction at or below which gaps score 1.0
    largest_gap          1000 bp Largest gap at or below which it scores 1.0

    Category thresholds on the overall score: Excellent >= 0.85, Good >= 0.70,
    Fair >= 0.50, otherwise Poor.
    """

    DEFAULT_WEIGHTS = {
        'coverage_breadth': 0.30,
        'coverage_depth': 0.20,
        'mapping_efficiency': 0.15,
        'gap_characteristics': 0.20,
        'reference_difficulty': 0.15,
    }

    DEFAULT_BENCHMARKS = {
        'coverage_breadth': 95.0,
        'mean_depth': 20.0,
        'depth_uniformity': 0.8,
        'mapping_efficiency': 90.0,
        'gap_percentage': 2.0,
        'largest_gap': 1000,
    }

    DEFAULT_THRESHOLDS = {
        'excellent': 0.85,
        'good': 0.70,
        'fair': 0.50,
        'poor': 0.0,
    }
    
    def __init__(
        self,
        coverage_stats: Dict[str, Any],
        gap_analysis: Dict[str, Any],
        reference_analysis: Dict[str, Any] = None,
        weights: Optional[Dict[str, float]] = None,
        benchmarks: Optional[Dict[str, float]] = None,
    ):
        """
        Initialize the quality scorer.
        
        Args:
            coverage_stats: Coverage statistics from CoverageAnalyzer
            gap_analysis: Gap analysis results from GapAnalyzer
            reference_analysis: Reference sequence analysis (optional)
            weights: Overrides for DEFAULT_WEIGHTS
            benchmarks: Overrides for DEFAULT_BENCHMARKS
        """
        self.coverage_stats = coverage_stats
        self.gap_analysis = gap_analysis
        self.reference_analysis = reference_analysis or {}
        
        self.weights = dict(self.DEFAULT_WEIGHTS)
        if weights:
            self.weights.update(weights)
        self.benchmarks = dict(self.DEFAULT_BENCHMARKS)
        if benchmarks:
            self.benchmarks.update(benchmarks)
        self.thresholds = dict(self.DEFAULT_THRESHOLDS)
    
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
    
