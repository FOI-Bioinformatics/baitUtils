#!/usr/bin/env python3

"""
benchmark.py

Benchmarking system for comparing actual coverage performance against theoretical optimal.
Provides insights into how well the current oligo set performs compared to ideal scenarios.

This module calculates theoretical optimal coverage metrics and compares them with
actual performance to identify improvement potential and design efficiency.
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Optional, NamedTuple
from pathlib import Path
from Bio import SeqIO
import math

from baitUtils.quality_scorer import QualityScore


class BenchmarkResult(NamedTuple):
    """Container for benchmark comparison results."""
    actual_score: float
    theoretical_optimal: float
    efficiency_ratio: float
    improvement_potential: float
    category: str
    recommendations: List[str]


class TheoreticalOptimal(NamedTuple):
    """Container for theoretical optimal metrics."""
    max_coverage_breadth: float
    optimal_depth_uniformity: float
    min_possible_gaps: int
    optimal_mapping_efficiency: float
    theoretical_quality_score: float


class BenchmarkAnalyzer:
    """
    Analyzes coverage performance against theoretical optimal benchmarks.
    
    Calculates what would be theoretically possible given the reference sequence
    characteristics and oligo design constraints, then compares actual performance.
    """
    
    def __init__(self, coverage_stats: Dict, gap_analysis: Dict, 
                 reference_analysis: Dict, quality_score: QualityScore,
                 oligo_length: int = 120, target_spacing: int = 500):
        """
        Initialize benchmark analyzer.
        
        Args:
            coverage_stats: Coverage analysis results
            gap_analysis: Gap analysis results  
            reference_analysis: Reference sequence analysis results
            quality_score: Calculated quality score
            oligo_length: Average oligo length for calculations
            target_spacing: Target spacing between oligos
        """
        self.coverage_stats = coverage_stats
        self.gap_analysis = gap_analysis
        self.reference_analysis = reference_analysis
        self.quality_score = quality_score
        self.oligo_length = oligo_length
        self.target_spacing = target_spacing
        
    def calculate_theoretical_optimal(self) -> TheoreticalOptimal:
        """Calculate theoretical optimal performance metrics."""
        
        # Get reference sequence characteristics
        ref_length = self.coverage_stats.get('reference_length', 0)
        challenging_regions = self.reference_analysis.get('challenging_regions', [])
        challenging_bp = sum(region['length'] for region in challenging_regions)
        
        # Calculate theoretical maximum coverage breadth
        # Assume we can cover all non-challenging regions perfectly
        max_coverage_breadth = ((ref_length - challenging_bp) / ref_length * 100) if ref_length > 0 else 0
        
        # Optimal depth uniformity (theoretical Gini coefficient approaching 0)
        optimal_depth_uniformity = 95.0  # Nearly perfect uniformity
        
        # Calculate minimum possible gaps
        # Based on challenging regions and unavoidable sequence gaps
        min_possible_gaps = len(challenging_regions)
        
        # Theoretical mapping efficiency
        # Assumes optimal oligo design with minimal off-target binding
        optimal_mapping_efficiency = 95.0
        
        # Calculate theoretical quality score
        # Based on optimal component scores
        theoretical_quality_score = self._calculate_theoretical_quality_score(
            max_coverage_breadth, optimal_depth_uniformity, 
            min_possible_gaps, optimal_mapping_efficiency
        )
        
        return TheoreticalOptimal(
            max_coverage_breadth=max_coverage_breadth,
            optimal_depth_uniformity=optimal_depth_uniformity,
            min_possible_gaps=min_possible_gaps,
            optimal_mapping_efficiency=optimal_mapping_efficiency,
            theoretical_quality_score=theoretical_quality_score
        )
    
    def _calculate_theoretical_quality_score(self, breadth: float, uniformity: float,
                                           gaps: int, mapping_eff: float) -> float:
        """Calculate theoretical optimal quality score."""
        
        # Use similar weighting as actual quality scorer
        breadth_score = min(breadth / 95.0 * 10, 10)  # Scale to /10
        depth_score = uniformity / 10  # Already scaled
        gap_score = max(10 - (gaps / 10), 1)  # Fewer gaps = higher score
        mapping_score = mapping_eff / 10  # Scale to /10
        
        # Weighted average (same weights as QualityScorer)
        weights = {'coverage': 0.3, 'depth': 0.25, 'gaps': 0.25, 'mapping': 0.2}
        
        theoretical_score = (
            breadth_score * weights['coverage'] +
            depth_score * weights['depth'] +
            gap_score * weights['gaps'] +
            mapping_score * weights['mapping']
        )
        
        return min(theoretical_score, 10.0)
    
    def benchmark_coverage_breadth(self, theoretical: TheoreticalOptimal) -> BenchmarkResult:
        """Benchmark coverage breadth against theoretical maximum."""
        
        actual_breadth = self.coverage_stats.get('coverage_breadth', 0)
        theoretical_breadth = theoretical.max_coverage_breadth
        
        if theoretical_breadth > 0:
            efficiency_ratio = actual_breadth / theoretical_breadth
            improvement_potential = theoretical_breadth - actual_breadth
        else:
            efficiency_ratio = 1.0
            improvement_potential = 0
        
        # Categorize performance
        if efficiency_ratio >= 0.9:
            category = "Excellent"
        elif efficiency_ratio >= 0.8:
            category = "Good" 
        elif efficiency_ratio >= 0.7:
            category = "Fair"
        else:
            category = "Poor"
        
        # Generate recommendations
        recommendations = []
        if improvement_potential > 10:
            recommendations.append(
                f"Coverage breadth could be improved by {improvement_potential:.1f}% "
                f"by targeting challenging regions more effectively"
            )
        if efficiency_ratio < 0.8:
            recommendations.append(
                "Consider adding more oligos to increase coverage breadth, "
                "focusing on the largest uncovered regions"
            )
        
        return BenchmarkResult(
            actual_score=actual_breadth,
            theoretical_optimal=theoretical_breadth,
            efficiency_ratio=efficiency_ratio,
            improvement_potential=improvement_potential,
            category=category,
            recommendations=recommendations
        )
    
    def benchmark_depth_uniformity(self, theoretical: TheoreticalOptimal) -> BenchmarkResult:
        """Benchmark depth uniformity against theoretical optimal."""
        
        # Convert Gini coefficient to uniformity score
        gini = self.coverage_stats.get('gini_coefficient', 0.5)
        actual_uniformity = (1 - gini) * 100  # Higher = more uniform
        theoretical_uniformity = theoretical.optimal_depth_uniformity
        
        if theoretical_uniformity > 0:
            efficiency_ratio = actual_uniformity / theoretical_uniformity
            improvement_potential = theoretical_uniformity - actual_uniformity
        else:
            efficiency_ratio = 1.0
            improvement_potential = 0
        
        # Categorize performance
        if efficiency_ratio >= 0.85:
            category = "Excellent"
        elif efficiency_ratio >= 0.7:
            category = "Good"
        elif efficiency_ratio >= 0.55:
            category = "Fair"
        else:
            category = "Poor"
        
        # Generate recommendations
        recommendations = []
        if improvement_potential > 15:
            recommendations.append(
                f"Coverage uniformity could be improved by {improvement_potential:.1f}% "
                f"by redistributing oligo density more evenly"
            )
        if efficiency_ratio < 0.7:
            recommendations.append(
                "High coverage variability detected. Consider rebalancing oligo "
                "distribution to achieve more uniform coverage"
            )
        
        return BenchmarkResult(
            actual_score=actual_uniformity,
            theoretical_optimal=theoretical_uniformity,
            efficiency_ratio=efficiency_ratio,
            improvement_potential=improvement_potential,
            category=category,
            recommendations=recommendations
        )
    
    def benchmark_gap_reduction(self, theoretical: TheoreticalOptimal) -> BenchmarkResult:
        """Benchmark gap count against theoretical minimum."""
        
        actual_gaps = self.gap_analysis.get('total_gaps', 0)
        theoretical_gaps = theoretical.min_possible_gaps
        
        # Calculate efficiency (lower gaps = higher efficiency)
        if actual_gaps > 0:
            efficiency_ratio = max(theoretical_gaps / actual_gaps, 0)
            improvement_potential = actual_gaps - theoretical_gaps
        else:
            efficiency_ratio = 1.0
            improvement_potential = 0
        
        # Categorize performance
        if efficiency_ratio >= 0.8:
            category = "Excellent"
        elif efficiency_ratio >= 0.6:
            category = "Good"
        elif efficiency_ratio >= 0.4:
            category = "Fair"
        else:
            category = "Poor"
        
        # Generate recommendations
        recommendations = []
        if improvement_potential > 50:
            recommendations.append(
                f"Gap count could potentially be reduced by {improvement_potential} gaps "
                f"through strategic oligo placement"
            )
        if efficiency_ratio < 0.6:
            recommendations.append(
                "High number of coverage gaps detected. Use 'baitUtils fill' "
                "to systematically close gaps with additional oligos"
            )
        
        return BenchmarkResult(
            actual_score=float(actual_gaps),
            theoretical_optimal=float(theoretical_gaps),
            efficiency_ratio=efficiency_ratio,
            improvement_potential=improvement_potential,
            category=category,
            recommendations=recommendations
        )
    
    def benchmark_overall_quality(self, theoretical: TheoreticalOptimal) -> BenchmarkResult:
        """Benchmark overall quality score against theoretical optimal."""
        
        actual_quality = self.quality_score.overall_score
        theoretical_quality = theoretical.theoretical_quality_score
        
        if theoretical_quality > 0:
            efficiency_ratio = actual_quality / theoretical_quality
            improvement_potential = theoretical_quality - actual_quality
        else:
            efficiency_ratio = 1.0
            improvement_potential = 0
        
        # Categorize performance
        if efficiency_ratio >= 0.9:
            category = "Excellent"
        elif efficiency_ratio >= 0.8:
            category = "Good"
        elif efficiency_ratio >= 0.7:
            category = "Fair"
        else:
            category = "Poor"
        
        # Generate recommendations based on quality score components
        recommendations = []
        if improvement_potential > 1.0:
            recommendations.append(
                f"Overall quality could be improved by {improvement_potential:.1f} points "
                f"({improvement_potential/10*100:.1f}%) through targeted optimization"
            )
        
        # Add specific recommendations based on component scores
        component_recs = self._generate_component_recommendations()
        recommendations.extend(component_recs)
        
        return BenchmarkResult(
            actual_score=actual_quality,
            theoretical_optimal=theoretical_quality,
            efficiency_ratio=efficiency_ratio,
            improvement_potential=improvement_potential,
            category=category,
            recommendations=recommendations
        )
    
    def _generate_component_recommendations(self) -> List[str]:
        """Generate recommendations based on quality score components."""
        recommendations = []
        
        # Check which components are underperforming
        if hasattr(self.quality_score, 'component_scores'):
            components = self.quality_score.component_scores
            
            if components.get('coverage_score', 10) < 7:
                recommendations.append(
                    "Coverage breadth is limiting overall quality - add oligos to uncovered regions"
                )
            
            if components.get('depth_uniformity_score', 10) < 7:
                recommendations.append(
                    "Coverage uniformity is limiting quality - redistribute oligo density"
                )
            
            if components.get('gap_score', 10) < 7:
                recommendations.append(
                    "Gap count is limiting quality - use gap-filling strategies"
                )
            
            if components.get('mapping_efficiency_score', 10) < 7:
                recommendations.append(
                    "Mapping efficiency is limiting quality - review oligo design parameters"
                )
        
        return recommendations
    
    def run_full_benchmark(self) -> Dict[str, BenchmarkResult]:
        """Run complete benchmark analysis."""
        
        logging.info("Calculating theoretical optimal performance...")
        theoretical = self.calculate_theoretical_optimal()
        
        logging.info("Benchmarking coverage performance...")
        benchmarks = {
            'coverage_breadth': self.benchmark_coverage_breadth(theoretical),
            'depth_uniformity': self.benchmark_depth_uniformity(theoretical),
            'gap_reduction': self.benchmark_gap_reduction(theoretical),
            'overall_quality': self.benchmark_overall_quality(theoretical)
        }
        
        return benchmarks
    
    def generate_benchmark_report(self, benchmarks: Dict[str, BenchmarkResult]) -> str:
        """Generate a formatted benchmark report."""
        
        report = []
        report.append("=" * 80)
        report.append("PERFORMANCE BENCHMARK ANALYSIS")
        report.append("=" * 80)
        report.append("")
        
        report.append("Performance vs Theoretical Optimal:")
        report.append("-" * 40)
        
        for metric, result in benchmarks.items():
            metric_name = metric.replace('_', ' ').title()
            report.append(f"{metric_name}:")
            report.append(f"  Actual:      {result.actual_score:8.1f}")
            report.append(f"  Theoretical: {result.theoretical_optimal:8.1f}")
            report.append(f"  Efficiency:  {result.efficiency_ratio*100:8.1f}%")
            report.append(f"  Category:    {result.category:>8s}")
            report.append("")
        
        report.append("Improvement Recommendations:")
        report.append("-" * 40)
        
        all_recommendations = []
        for result in benchmarks.values():
            all_recommendations.extend(result.recommendations)
        
        if all_recommendations:
            for i, rec in enumerate(all_recommendations, 1):
                report.append(f"{i}. {rec}")
        else:
            report.append("Performance is near optimal - no major improvements identified.")
        
        report.append("")
        report.append("=" * 80)
        
        return "\n".join(report)