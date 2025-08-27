#!/usr/bin/env python3

"""
differential_analysis.py

Statistical differential analysis tools for comparing oligo sets.
Provides statistical testing, effect size calculations, and significance assessment
for differences between oligo set performance metrics.

This module enables rigorous statistical comparison of oligo sets to identify
statistically significant differences in coverage patterns and quality metrics.
"""

import logging
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, NamedTuple
from dataclasses import dataclass
import scipy.stats as stats
from pathlib import Path

from baitUtils.comparative_analyzer import OligoSetResult, ComparisonMetrics


@dataclass
class StatisticalTest:
    """Container for statistical test results."""
    test_name: str
    statistic: float
    p_value: float
    effect_size: float
    significance_level: str  # 'ns', '*', '**', '***'
    interpretation: str


@dataclass
class CoverageDistributionComparison:
    """Container for coverage distribution comparison results."""
    set1_name: str
    set2_name: str
    ks_test: StatisticalTest
    mann_whitney_test: StatisticalTest
    levene_test: StatisticalTest  # Test for equal variances
    summary_stats: Dict[str, Dict[str, float]]


class DifferentialAnalyzer:
    """
    Statistical analysis toolkit for comparing oligo set performance.
    
    Provides comprehensive statistical testing including:
    - Coverage distribution comparisons
    - Gap pattern analysis
    - Quality metric significance testing
    - Effect size calculations
    - Multiple comparison corrections
    """
    
    def __init__(self, significance_level: float = 0.05):
        """
        Initialize differential analyzer.
        
        Args:
            significance_level: Statistical significance threshold
        """
        self.significance_level = significance_level
        self.alpha_levels = {
            0.001: '***',
            0.01: '**', 
            0.05: '*',
            1.0: 'ns'
        }
    
    def compare_coverage_distributions(self, set1: OligoSetResult, 
                                     set2: OligoSetResult) -> CoverageDistributionComparison:
        """
        Compare coverage depth distributions between two oligo sets.
        
        Args:
            set1: First oligo set results
            set2: Second oligo set results
            
        Returns:
            Statistical comparison results
        """
        # Extract coverage arrays (mock implementation - would need actual coverage arrays)
        # In real implementation, this would extract per-position coverage data
        coverage1 = self._extract_coverage_array(set1)
        coverage2 = self._extract_coverage_array(set2)
        
        # Kolmogorov-Smirnov test (distribution shape)
        ks_stat, ks_p = stats.kstest(coverage1, coverage2)
        ks_effect_size = self._calculate_ks_effect_size(ks_stat, len(coverage1), len(coverage2))
        
        ks_test = StatisticalTest(
            test_name="Kolmogorov-Smirnov",
            statistic=ks_stat,
            p_value=ks_p,
            effect_size=ks_effect_size,
            significance_level=self._get_significance_level(ks_p),
            interpretation=self._interpret_ks_test(ks_stat, ks_p)
        )
        
        # Mann-Whitney U test (median differences)
        mw_stat, mw_p = stats.mannwhitneyu(coverage1, coverage2, alternative='two-sided')
        mw_effect_size = self._calculate_mannwhitney_effect_size(mw_stat, len(coverage1), len(coverage2))
        
        mw_test = StatisticalTest(
            test_name="Mann-Whitney U",
            statistic=mw_stat,
            p_value=mw_p,
            effect_size=mw_effect_size,
            significance_level=self._get_significance_level(mw_p),
            interpretation=self._interpret_mannwhitney_test(mw_stat, mw_p)
        )
        
        # Levene's test for equal variances
        levene_stat, levene_p = stats.levene(coverage1, coverage2)
        levene_effect_size = self._calculate_levene_effect_size(levene_stat, len(coverage1), len(coverage2))
        
        levene_test = StatisticalTest(
            test_name="Levene's Test",
            statistic=levene_stat,
            p_value=levene_p,
            effect_size=levene_effect_size,
            significance_level=self._get_significance_level(levene_p),
            interpretation=self._interpret_levene_test(levene_stat, levene_p)
        )
        
        # Summary statistics
        summary_stats = {
            set1.name: {
                'mean': np.mean(coverage1),
                'median': np.median(coverage1),
                'std': np.std(coverage1),
                'q25': np.percentile(coverage1, 25),
                'q75': np.percentile(coverage1, 75),
                'skewness': stats.skew(coverage1),
                'kurtosis': stats.kurtosis(coverage1)
            },
            set2.name: {
                'mean': np.mean(coverage2),
                'median': np.median(coverage2), 
                'std': np.std(coverage2),
                'q25': np.percentile(coverage2, 25),
                'q75': np.percentile(coverage2, 75),
                'skewness': stats.skew(coverage2),
                'kurtosis': stats.kurtosis(coverage2)
            }
        }
        
        return CoverageDistributionComparison(
            set1_name=set1.name,
            set2_name=set2.name,
            ks_test=ks_test,
            mann_whitney_test=mw_test,
            levene_test=levene_test,
            summary_stats=summary_stats
        )
    
    def _extract_coverage_array(self, oligo_set: OligoSetResult) -> np.ndarray:
        """Extract coverage depth array from oligo set results."""
        # Mock implementation - in reality would extract per-position coverage
        # For now, simulate based on coverage statistics
        
        mean_depth = oligo_set.coverage_stats.get('mean_depth', 5.0)
        coverage_breadth = oligo_set.coverage_stats.get('coverage_breadth', 80.0)
        gini_coeff = oligo_set.coverage_stats.get('gini_coefficient', 0.3)
        
        # Simulate coverage array based on statistics
        ref_length = oligo_set.coverage_stats.get('reference_length', 10000)
        covered_positions = int(ref_length * coverage_breadth / 100)
        
        # Generate coverage values with appropriate distribution
        # Use gamma distribution to simulate realistic coverage
        if gini_coeff < 0.2:  # Uniform coverage
            coverage = np.random.gamma(shape=mean_depth**2, scale=1/mean_depth, size=covered_positions)
        elif gini_coeff > 0.5:  # Variable coverage
            coverage = np.random.gamma(shape=1, scale=mean_depth, size=covered_positions)
        else:  # Moderate variability
            coverage = np.random.gamma(shape=mean_depth, scale=1, size=covered_positions)
        
        # Add zeros for uncovered positions
        zeros = np.zeros(ref_length - covered_positions)
        full_coverage = np.concatenate([coverage, zeros])
        np.random.shuffle(full_coverage)
        
        return full_coverage
    
    def compare_quality_metrics(self, oligo_sets: List[OligoSetResult]) -> Dict[str, StatisticalTest]:
        """
        Compare quality metrics across multiple oligo sets.
        
        Args:
            oligo_sets: List of oligo set results to compare
            
        Returns:
            Dictionary of statistical tests for each metric
        """
        if len(oligo_sets) < 2:
            raise ValueError("Need at least 2 oligo sets for comparison")
        
        results = {}
        
        # Extract metrics for comparison
        metrics = {
            'quality_score': [s.quality_score.overall_score for s in oligo_sets],
            'coverage_breadth': [s.coverage_stats.get('coverage_breadth', 0) for s in oligo_sets],
            'mean_depth': [s.coverage_stats.get('mean_depth', 0) for s in oligo_sets],
            'mapping_efficiency': [s.coverage_stats.get('mapping_efficiency', 0) for s in oligo_sets],
            'gap_count': [s.gap_analysis.get('total_gaps', 0) for s in oligo_sets],
            'gini_coefficient': [s.coverage_stats.get('gini_coefficient', 0) for s in oligo_sets]
        }
        
        for metric_name, values in metrics.items():
            if len(oligo_sets) == 2:
                # Use t-test for two groups
                stat, p_value = stats.ttest_ind(values[:1], values[1:])
                effect_size = self._calculate_cohens_d(values[:1], values[1:])
                test_name = "Independent t-test"
            else:
                # Use ANOVA for multiple groups
                groups = [[v] for v in values]  # Each oligo set as a group
                stat, p_value = stats.f_oneway(*groups)
                effect_size = self._calculate_eta_squared(groups)
                test_name = "One-way ANOVA"
            
            results[metric_name] = StatisticalTest(
                test_name=test_name,
                statistic=stat,
                p_value=p_value,
                effect_size=effect_size,
                significance_level=self._get_significance_level(p_value),
                interpretation=self._interpret_metric_test(metric_name, stat, p_value)
            )
        
        return results
    
    def analyze_gap_patterns(self, set1: OligoSetResult, set2: OligoSetResult) -> Dict[str, StatisticalTest]:
        """
        Analyze differences in gap patterns between two oligo sets.
        
        Args:
            set1: First oligo set results
            set2: Second oligo set results
            
        Returns:
            Statistical tests for gap pattern differences
        """
        results = {}
        
        # Extract gap size distributions
        gaps1 = set1.gap_analysis.get('gaps', [])
        gaps2 = set2.gap_analysis.get('gaps', [])
        
        if not gaps1 or not gaps2:
            logging.warning("Insufficient gap data for statistical comparison")
            return results
        
        gap_sizes1 = [gap['length'] for gap in gaps1]
        gap_sizes2 = [gap['length'] for gap in gaps2]
        
        # Mann-Whitney test for gap size distributions
        if len(gap_sizes1) > 0 and len(gap_sizes2) > 0:
            mw_stat, mw_p = stats.mannwhitneyu(gap_sizes1, gap_sizes2, alternative='two-sided')
            mw_effect_size = self._calculate_mannwhitney_effect_size(mw_stat, len(gap_sizes1), len(gap_sizes2))
            
            results['gap_size_distribution'] = StatisticalTest(
                test_name="Mann-Whitney U (Gap Sizes)",
                statistic=mw_stat,
                p_value=mw_p,
                effect_size=mw_effect_size,
                significance_level=self._get_significance_level(mw_p),
                interpretation=self._interpret_gap_size_test(mw_stat, mw_p)
            )
        
        # Chi-square test for gap count differences
        gap_counts = [len(gaps1), len(gaps2)]
        total_gaps = sum(gap_counts)
        expected = [total_gaps / 2, total_gaps / 2]
        
        if total_gaps > 0:
            chi2_stat, chi2_p = stats.chisquare(gap_counts, expected)
            chi2_effect_size = np.sqrt(chi2_stat / total_gaps)
            
            results['gap_count_difference'] = StatisticalTest(
                test_name="Chi-square (Gap Count)",
                statistic=chi2_stat,
                p_value=chi2_p,
                effect_size=chi2_effect_size,
                significance_level=self._get_significance_level(chi2_p),
                interpretation=self._interpret_gap_count_test(chi2_stat, chi2_p)
            )
        
        return results
    
    def multiple_comparison_correction(self, p_values: List[float], 
                                     method: str = 'bonferroni') -> List[float]:
        """
        Apply multiple comparison correction to p-values.
        
        Args:
            p_values: List of uncorrected p-values
            method: Correction method ('bonferroni', 'holm', 'fdr')
            
        Returns:
            Corrected p-values
        """
        p_array = np.array(p_values)
        
        if method == 'bonferroni':
            corrected = p_array * len(p_values)
            corrected = np.minimum(corrected, 1.0)
        elif method == 'holm':
            corrected = self._holm_correction(p_array)
        elif method == 'fdr':
            corrected = self._fdr_correction(p_array)
        else:
            raise ValueError(f"Unknown correction method: {method}")
        
        return corrected.tolist()
    
    def _holm_correction(self, p_values: np.ndarray) -> np.ndarray:
        """Apply Holm-Bonferroni correction."""
        sorted_indices = np.argsort(p_values)
        sorted_p = p_values[sorted_indices]
        n = len(p_values)
        
        corrected = np.zeros_like(sorted_p)
        for i, p in enumerate(sorted_p):
            corrected[i] = min(1.0, p * (n - i))
        
        # Ensure monotonicity
        for i in range(1, len(corrected)):
            corrected[i] = max(corrected[i], corrected[i-1])
        
        # Restore original order
        result = np.zeros_like(corrected)
        result[sorted_indices] = corrected
        return result
    
    def _fdr_correction(self, p_values: np.ndarray) -> np.ndarray:
        """Apply False Discovery Rate (Benjamini-Hochberg) correction."""
        sorted_indices = np.argsort(p_values)
        sorted_p = p_values[sorted_indices]
        n = len(p_values)
        
        corrected = np.zeros_like(sorted_p)
        for i in range(n-1, -1, -1):
            if i == n-1:
                corrected[i] = sorted_p[i]
            else:
                corrected[i] = min(corrected[i+1], sorted_p[i] * n / (i+1))
        
        # Restore original order
        result = np.zeros_like(corrected)
        result[sorted_indices] = corrected
        return result
    
    # Effect size calculation methods
    def _calculate_cohens_d(self, group1: List[float], group2: List[float]) -> float:
        """Calculate Cohen's d effect size."""
        mean1, mean2 = np.mean(group1), np.mean(group2)
        pooled_std = np.sqrt(((len(group1) - 1) * np.var(group1, ddof=1) + 
                             (len(group2) - 1) * np.var(group2, ddof=1)) / 
                            (len(group1) + len(group2) - 2))
        return (mean1 - mean2) / pooled_std if pooled_std > 0 else 0
    
    def _calculate_eta_squared(self, groups: List[List[float]]) -> float:
        """Calculate eta-squared effect size for ANOVA."""
        # Simplified implementation
        all_values = [val for group in groups for val in group]
        grand_mean = np.mean(all_values)
        
        ss_between = sum(len(group) * (np.mean(group) - grand_mean)**2 for group in groups)
        ss_total = sum((val - grand_mean)**2 for val in all_values)
        
        return ss_between / ss_total if ss_total > 0 else 0
    
    def _calculate_ks_effect_size(self, ks_stat: float, n1: int, n2: int) -> float:
        """Calculate effect size for KS test."""
        return ks_stat  # KS statistic itself is a measure of effect size
    
    def _calculate_mannwhitney_effect_size(self, mw_stat: float, n1: int, n2: int) -> float:
        """Calculate effect size for Mann-Whitney test."""
        return (mw_stat / (n1 * n2)) - 0.5  # r = (U / (n1 * n2)) - 0.5
    
    def _calculate_levene_effect_size(self, levene_stat: float, n1: int, n2: int) -> float:
        """Calculate effect size for Levene's test."""
        # Use eta-squared approximation
        df_between = 1
        df_total = n1 + n2 - 1
        return (levene_stat * df_between) / (levene_stat * df_between + df_total)
    
    # Significance level determination
    def _get_significance_level(self, p_value: float) -> str:
        """Determine significance level based on p-value."""
        for threshold, level in self.alpha_levels.items():
            if p_value <= threshold:
                return level
        return 'ns'
    
    # Interpretation methods
    def _interpret_ks_test(self, statistic: float, p_value: float) -> str:
        """Interpret KS test results."""
        if p_value <= self.significance_level:
            return f"Significant difference in coverage distributions (D={statistic:.3f})"
        else:
            return f"No significant difference in coverage distributions (D={statistic:.3f})"
    
    def _interpret_mannwhitney_test(self, statistic: float, p_value: float) -> str:
        """Interpret Mann-Whitney test results."""
        if p_value <= self.significance_level:
            return f"Significant difference in coverage medians (U={statistic:.0f})"
        else:
            return f"No significant difference in coverage medians (U={statistic:.0f})"
    
    def _interpret_levene_test(self, statistic: float, p_value: float) -> str:
        """Interpret Levene's test results."""
        if p_value <= self.significance_level:
            return f"Significant difference in coverage variability (W={statistic:.3f})"
        else:
            return f"No significant difference in coverage variability (W={statistic:.3f})"
    
    def _interpret_metric_test(self, metric: str, statistic: float, p_value: float) -> str:
        """Interpret metric comparison test results."""
        metric_names = {
            'quality_score': 'quality scores',
            'coverage_breadth': 'coverage breadth',
            'mean_depth': 'mean coverage depth',
            'mapping_efficiency': 'mapping efficiency',
            'gap_count': 'gap counts',
            'gini_coefficient': 'coverage uniformity'
        }
        
        metric_name = metric_names.get(metric, metric)
        
        if p_value <= self.significance_level:
            return f"Significant difference in {metric_name} (stat={statistic:.3f})"
        else:
            return f"No significant difference in {metric_name} (stat={statistic:.3f})"
    
    def _interpret_gap_size_test(self, statistic: float, p_value: float) -> str:
        """Interpret gap size distribution test results."""
        if p_value <= self.significance_level:
            return f"Significant difference in gap size distributions (U={statistic:.0f})"
        else:
            return f"No significant difference in gap size distributions (U={statistic:.0f})"
    
    def _interpret_gap_count_test(self, statistic: float, p_value: float) -> str:
        """Interpret gap count test results."""
        if p_value <= self.significance_level:
            return f"Significant difference in gap counts (χ²={statistic:.3f})"
        else:
            return f"No significant difference in gap counts (χ²={statistic:.3f})"