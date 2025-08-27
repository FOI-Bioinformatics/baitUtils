#!/usr/bin/env python3

"""
comparative_visualizations.py

Comprehensive visualization suite for comparing multiple oligo sets.
Provides side-by-side comparisons, statistical plots, and interactive dashboards
for analyzing differences between oligo set performance.

This module creates publication-quality comparative plots and interactive
dashboards for Phase 3 comparative analysis functionality.
"""

import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.figure_factory as ff

from baitUtils.comparative_analyzer import OligoSetResult, ComparativeAnalyzer
from baitUtils.differential_analysis import DifferentialAnalyzer, CoverageDistributionComparison


class ComparativeVisualizer:
    """
    Comprehensive visualization suite for comparative oligo set analysis.
    
    Provides both static (matplotlib/seaborn) and interactive (plotly) visualizations
    for comparing multiple oligo sets across various performance metrics.
    """
    
    def __init__(self, output_dir: Path, plot_format: str = 'png', 
                 dpi: int = 300, style: str = 'whitegrid'):
        """
        Initialize comparative visualizer.
        
        Args:
            output_dir: Directory for saving plots
            plot_format: Plot format ('png', 'pdf', 'svg')
            dpi: Plot resolution
            style: Seaborn style theme
        """
        self.output_dir = Path(output_dir)
        self.plot_format = plot_format
        self.dpi = dpi
        
        # Create plots directory
        self.plots_dir = self.output_dir / "comparative_plots"
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        
        self.interactive_dir = self.output_dir / "interactive_comparative_plots"
        self.interactive_dir.mkdir(parents=True, exist_ok=True)
        
        # Set style
        sns.set_style(style)
        plt.rcParams['figure.dpi'] = dpi
        plt.rcParams['savefig.dpi'] = dpi
    
    def create_comparison_dashboard(self, analyzer: ComparativeAnalyzer) -> str:
        """Create comprehensive comparison dashboard."""
        
        if not analyzer.oligo_sets:
            raise ValueError("No oligo sets found in analyzer")
        
        # Generate comparison matrix
        comparison_matrix = analyzer.generate_comparison_matrix()
        
        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=(
                'Quality Score Comparison', 'Coverage Breadth Comparison',
                'Gap Count Comparison', 'Mapping Efficiency Comparison',
                'Coverage Depth Comparison', 'Performance Ranking'
            ),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}]]
        )
        
        set_names = comparison_matrix['Name'].tolist()
        colors = px.colors.qualitative.Set3[:len(set_names)]
        
        # Quality Score Comparison
        fig.add_trace(
            go.Bar(
                name='Quality Score',
                x=set_names,
                y=comparison_matrix['Quality_Score'].tolist(),
                marker_color=colors,
                text=comparison_matrix['Quality_Grade'].tolist(),
                textposition='auto'
            ),
            row=1, col=1
        )
        
        # Coverage Breadth Comparison
        fig.add_trace(
            go.Bar(
                name='Coverage Breadth',
                x=set_names,
                y=comparison_matrix['Coverage_Breadth_%'].tolist(),
                marker_color=colors,
                text=[f"{x:.1f}%" for x in comparison_matrix['Coverage_Breadth_%']],
                textposition='auto'
            ),
            row=1, col=2
        )
        
        # Gap Count Comparison
        fig.add_trace(
            go.Bar(
                name='Total Gaps',
                x=set_names,
                y=comparison_matrix['Total_Gaps'].tolist(),
                marker_color=colors,
                text=comparison_matrix['Total_Gaps'].tolist(),
                textposition='auto'
            ),
            row=2, col=1
        )
        
        # Mapping Efficiency Comparison
        fig.add_trace(
            go.Bar(
                name='Mapping Efficiency',
                x=set_names,
                y=comparison_matrix['Mapping_Efficiency_%'].tolist(),
                marker_color=colors,
                text=[f"{x:.1f}%" for x in comparison_matrix['Mapping_Efficiency_%']],
                textposition='auto'
            ),
            row=2, col=2
        )
        
        # Coverage Depth Comparison
        fig.add_trace(
            go.Bar(
                name='Mean Depth',
                x=set_names,
                y=comparison_matrix['Mean_Depth_x'].tolist(),
                marker_color=colors,
                text=[f"{x:.1f}x" for x in comparison_matrix['Mean_Depth_x']],
                textposition='auto'
            ),
            row=3, col=1
        )
        
        # Performance Ranking (radar-like using polar bar chart simulation)
        ranking = analyzer.generate_ranking()
        rank_names, rank_scores = zip(*ranking) if ranking else ([], [])
        
        fig.add_trace(
            go.Bar(
                name='Composite Score',
                x=list(rank_names),
                y=list(rank_scores),
                marker_color=colors[:len(rank_names)],
                text=[f"{x:.2f}" for x in rank_scores],
                textposition='auto'
            ),
            row=3, col=2
        )
        
        fig.update_layout(
            height=1200,
            title_text="Oligo Set Comparative Analysis Dashboard",
            showlegend=False
        )
        
        # Save interactive plot
        dashboard_file = self.interactive_dir / "comparison_dashboard.html"
        fig.write_html(str(dashboard_file))
        
        return str(dashboard_file)
    
    def create_coverage_distribution_plots(self, oligo_sets: List[OligoSetResult]) -> str:
        """Create coverage distribution comparison plots."""
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Coverage Distribution Comparisons', fontsize=16, fontweight='bold')
        
        # Extract coverage statistics for each set
        data = []
        for result in oligo_sets:
            data.append({
                'Name': result.name,
                'Coverage_Breadth': result.coverage_stats.get('coverage_breadth', 0),
                'Mean_Depth': result.coverage_stats.get('mean_depth', 0),
                'Median_Depth': result.coverage_stats.get('median_depth', 0),
                'Gini_Coefficient': result.coverage_stats.get('gini_coefficient', 0)
            })
        
        df = pd.DataFrame(data)
        
        # Coverage Breadth comparison
        sns.barplot(data=df, x='Name', y='Coverage_Breadth', ax=axes[0,0])
        axes[0,0].set_title('Coverage Breadth Comparison')
        axes[0,0].set_ylabel('Coverage Breadth (%)')
        axes[0,0].tick_params(axis='x', rotation=45)
        
        # Add values on bars
        for i, v in enumerate(df['Coverage_Breadth']):
            axes[0,0].text(i, v + 1, f'{v:.1f}%', ha='center', va='bottom')
        
        # Mean vs Median Depth
        x = np.arange(len(df))
        width = 0.35
        
        axes[0,1].bar(x - width/2, df['Mean_Depth'], width, label='Mean Depth', alpha=0.8)
        axes[0,1].bar(x + width/2, df['Median_Depth'], width, label='Median Depth', alpha=0.8)
        axes[0,1].set_title('Coverage Depth Distribution')
        axes[0,1].set_ylabel('Coverage Depth (x)')
        axes[0,1].set_xticks(x)
        axes[0,1].set_xticklabels(df['Name'], rotation=45)
        axes[0,1].legend()
        
        # Coverage Uniformity (Gini Coefficient)
        sns.barplot(data=df, x='Name', y='Gini_Coefficient', ax=axes[1,0])
        axes[1,0].set_title('Coverage Uniformity (Gini Coefficient)')
        axes[1,0].set_ylabel('Gini Coefficient (lower = more uniform)')
        axes[1,0].tick_params(axis='x', rotation=45)
        
        # Add values on bars
        for i, v in enumerate(df['Gini_Coefficient']):
            axes[1,0].text(i, v + 0.01, f'{v:.3f}', ha='center', va='bottom')
        
        # Quality Score vs Coverage Breadth scatter
        quality_scores = [result.quality_score.overall_score for result in oligo_sets]
        coverage_breadths = [result.coverage_stats.get('coverage_breadth', 0) for result in oligo_sets]
        
        scatter = axes[1,1].scatter(coverage_breadths, quality_scores, 
                                  s=100, alpha=0.7, c=range(len(oligo_sets)), 
                                  cmap='viridis')
        axes[1,1].set_xlabel('Coverage Breadth (%)')
        axes[1,1].set_ylabel('Quality Score')
        axes[1,1].set_title('Quality Score vs Coverage Breadth')
        
        # Add labels to points
        for i, (x, y) in enumerate(zip(coverage_breadths, quality_scores)):
            axes[1,1].annotate(oligo_sets[i].name, (x, y), 
                             xytext=(5, 5), textcoords='offset points', fontsize=8)
        
        plt.tight_layout()
        
        # Save plot
        coverage_file = self.plots_dir / f"coverage_distributions.{self.plot_format}"
        plt.savefig(coverage_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return str(coverage_file)
    
    def create_gap_analysis_comparison(self, oligo_sets: List[OligoSetResult]) -> str:
        """Create gap analysis comparison plots."""
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Gap Analysis Comparison', fontsize=16, fontweight='bold')
        
        # Extract gap data
        gap_data = []
        for result in oligo_sets:
            gap_data.append({
                'Name': result.name,
                'Total_Gaps': result.gap_analysis.get('total_gaps', 0),
                'Mean_Gap_Size': result.gap_analysis.get('mean_gap_size', 0),
                'Max_Gap_Size': result.gap_analysis.get('max_gap_size', 0),
                'Gap_Percentage': result.gap_analysis.get('gap_percentage', 0)
            })
        
        gap_df = pd.DataFrame(gap_data)
        
        # Total Gaps comparison
        bars = axes[0,0].bar(gap_df['Name'], gap_df['Total_Gaps'])
        axes[0,0].set_title('Total Gap Count')
        axes[0,0].set_ylabel('Number of Gaps')
        axes[0,0].tick_params(axis='x', rotation=45)
        
        # Color bars based on values (red = more gaps)
        max_gaps = gap_df['Total_Gaps'].max()
        for bar, gaps in zip(bars, gap_df['Total_Gaps']):
            normalized_color = gaps / max_gaps if max_gaps > 0 else 0
            bar.set_color(plt.cm.Reds(0.3 + 0.7 * normalized_color))
        
        # Add values on bars
        for i, v in enumerate(gap_df['Total_Gaps']):
            axes[0,0].text(i, v + max_gaps * 0.01, str(int(v)), ha='center', va='bottom')
        
        # Gap size distribution
        axes[0,1].bar(gap_df['Name'], gap_df['Mean_Gap_Size'], alpha=0.7, label='Mean Gap Size')
        axes[0,1].set_title('Mean Gap Size')
        axes[0,1].set_ylabel('Gap Size (bp)')
        axes[0,1].tick_params(axis='x', rotation=45)
        
        # Maximum gap size comparison
        bars = axes[1,0].bar(gap_df['Name'], gap_df['Max_Gap_Size'])
        axes[1,0].set_title('Largest Gap Size')
        axes[1,0].set_ylabel('Gap Size (bp)')
        axes[1,0].tick_params(axis='x', rotation=45)
        
        # Color bars based on values (red = larger gaps)
        max_gap_size = gap_df['Max_Gap_Size'].max()
        for bar, gap_size in zip(bars, gap_df['Max_Gap_Size']):
            normalized_color = gap_size / max_gap_size if max_gap_size > 0 else 0
            bar.set_color(plt.cm.Oranges(0.3 + 0.7 * normalized_color))
        
        # Gap percentage (proportion of reference uncovered)
        sns.barplot(data=gap_df, x='Name', y='Gap_Percentage', ax=axes[1,1])
        axes[1,1].set_title('Percentage of Reference Uncovered')
        axes[1,1].set_ylabel('Gap Percentage (%)')
        axes[1,1].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        # Save plot
        gap_file = self.plots_dir / f"gap_analysis_comparison.{self.plot_format}"
        plt.savefig(gap_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return str(gap_file)
    
    def create_performance_heatmap(self, analyzer: ComparativeAnalyzer) -> str:
        """Create performance heatmap comparing all metrics."""
        
        comparison_matrix = analyzer.generate_comparison_matrix()
        
        # Select metrics for heatmap
        metrics = ['Coverage_Breadth_%', 'Mean_Depth_x', 'Mapping_Efficiency_%', 
                  'Quality_Score', 'Total_Gaps', 'Gini_Coefficient']
        
        heatmap_data = comparison_matrix[['Name'] + metrics].copy()
        
        # Normalize metrics (0-1 scale, higher is better)
        for metric in metrics:
            if metric in ['Total_Gaps', 'Gini_Coefficient']:  # Lower is better
                heatmap_data[metric] = 1 - (heatmap_data[metric] / heatmap_data[metric].max())
            else:  # Higher is better
                heatmap_data[metric] = heatmap_data[metric] / heatmap_data[metric].max()
        
        # Set names as index
        heatmap_data = heatmap_data.set_index('Name')
        
        # Create heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(heatmap_data.T, annot=True, cmap='RdYlGn', fmt='.3f',
                   cbar_kws={'label': 'Normalized Performance (0-1)'})
        plt.title('Oligo Set Performance Heatmap\n(Normalized Metrics: Green=Better, Red=Worse)')
        plt.ylabel('Performance Metrics')
        plt.xlabel('Oligo Sets')
        plt.tight_layout()
        
        # Save plot
        heatmap_file = self.plots_dir / f"performance_heatmap.{self.plot_format}"
        plt.savefig(heatmap_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return str(heatmap_file)
    
    def create_statistical_comparison_plot(self, differential_analyzer: DifferentialAnalyzer,
                                         oligo_sets: List[OligoSetResult]) -> str:
        """Create statistical comparison visualization."""
        
        if len(oligo_sets) < 2:
            raise ValueError("Need at least 2 oligo sets for statistical comparison")
        
        # Perform statistical tests
        quality_tests = differential_analyzer.compare_quality_metrics(oligo_sets)
        
        # Create significance plot
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Statistical Comparison Results', fontsize=16, fontweight='bold')
        
        # P-values plot
        metrics = list(quality_tests.keys())
        p_values = [test.p_value for test in quality_tests.values()]
        effect_sizes = [abs(test.effect_size) for test in quality_tests.values()]
        
        bars = axes[0,0].bar(metrics, [-np.log10(p) for p in p_values])
        axes[0,0].axhline(y=-np.log10(0.05), color='red', linestyle='--', 
                         label='α = 0.05')
        axes[0,0].axhline(y=-np.log10(0.01), color='orange', linestyle='--', 
                         label='α = 0.01')
        axes[0,0].set_ylabel('-log₁₀(p-value)')
        axes[0,0].set_title('Statistical Significance')
        axes[0,0].tick_params(axis='x', rotation=45)
        axes[0,0].legend()
        
        # Color bars by significance
        for bar, p_val in zip(bars, p_values):
            if p_val <= 0.001:
                bar.set_color('darkgreen')
            elif p_val <= 0.01:
                bar.set_color('green')
            elif p_val <= 0.05:
                bar.set_color('orange')
            else:
                bar.set_color('gray')
        
        # Effect sizes plot
        bars = axes[0,1].bar(metrics, effect_sizes)
        axes[0,1].set_ylabel('Effect Size (|d|)')
        axes[0,1].set_title('Effect Sizes')
        axes[0,1].tick_params(axis='x', rotation=45)
        
        # Color bars by effect size magnitude
        for bar, effect in zip(bars, effect_sizes):
            if effect >= 0.8:  # Large effect
                bar.set_color('darkblue')
            elif effect >= 0.5:  # Medium effect
                bar.set_color('blue')
            elif effect >= 0.2:  # Small effect
                bar.set_color('lightblue')
            else:  # Negligible effect
                bar.set_color('gray')
        
        # Pairwise comparison (if 2 sets)
        if len(oligo_sets) == 2:
            # Coverage distribution comparison
            dist_comparison = differential_analyzer.compare_coverage_distributions(
                oligo_sets[0], oligo_sets[1]
            )
            
            # Plot summary statistics
            stats_data = pd.DataFrame(dist_comparison.summary_stats).T
            stats_data[['mean', 'median']].plot(kind='bar', ax=axes[1,0])
            axes[1,0].set_title('Coverage Distribution Statistics')
            axes[1,0].set_ylabel('Coverage Depth')
            axes[1,0].tick_params(axis='x', rotation=45)
            axes[1,0].legend()
            
            # Test results summary
            test_results = [
                ('KS Test', dist_comparison.ks_test.p_value, dist_comparison.ks_test.effect_size),
                ('Mann-Whitney', dist_comparison.mann_whitney_test.p_value, 
                 dist_comparison.mann_whitney_test.effect_size),
                ('Levene', dist_comparison.levene_test.p_value, 
                 dist_comparison.levene_test.effect_size)
            ]
            
            test_names, test_p_vals, test_effects = zip(*test_results)
            
            x_pos = np.arange(len(test_names))
            axes[1,1].bar(x_pos, [-np.log10(p) for p in test_p_vals], alpha=0.7)
            axes[1,1].axhline(y=-np.log10(0.05), color='red', linestyle='--')
            axes[1,1].set_xticks(x_pos)
            axes[1,1].set_xticklabels(test_names)
            axes[1,1].set_ylabel('-log₁₀(p-value)')
            axes[1,1].set_title('Distribution Tests')
        
        else:
            # Multiple sets - show ANOVA results
            axes[1,0].text(0.5, 0.5, 'ANOVA Results\n(Multiple Groups)', 
                          ha='center', va='center', transform=axes[1,0].transAxes,
                          fontsize=14)
            axes[1,0].set_title('Multiple Group Analysis')
            
            axes[1,1].text(0.5, 0.5, 'See Quality Tests\nfor Detailed Results', 
                          ha='center', va='center', transform=axes[1,1].transAxes,
                          fontsize=14)
            axes[1,1].set_title('Statistical Summary')
        
        plt.tight_layout()
        
        # Save plot
        stats_file = self.plots_dir / f"statistical_comparison.{self.plot_format}"
        plt.savefig(stats_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return str(stats_file)
    
    def create_interactive_ranking_plot(self, analyzer: ComparativeAnalyzer) -> str:
        """Create interactive ranking visualization."""
        
        ranking = analyzer.generate_ranking()
        if not ranking:
            raise ValueError("No ranking data available")
        
        names, scores = zip(*ranking)
        
        # Create radar chart for top performers
        comparison_matrix = analyzer.generate_comparison_matrix()
        
        fig = go.Figure()
        
        # Add traces for each oligo set
        for _, row in comparison_matrix.iterrows():
            fig.add_trace(go.Scatterpolar(
                r=[
                    row['Coverage_Breadth_%'] / 100,  # Normalize to 0-1
                    row['Mean_Depth_x'] / 20,  # Normalize assuming max ~20x
                    row['Mapping_Efficiency_%'] / 100,
                    row['Quality_Score'] / 10,
                    1 - (row['Total_Gaps'] / comparison_matrix['Total_Gaps'].max()),  # Invert
                    1 - row['Gini_Coefficient']  # Invert (higher = more uniform)
                ],
                theta=[
                    'Coverage<br>Breadth',
                    'Mean<br>Depth', 
                    'Mapping<br>Efficiency',
                    'Quality<br>Score',
                    'Gap<br>Reduction',
                    'Coverage<br>Uniformity'
                ],
                fill='toself',
                name=row['Name']
            ))
        
        fig.update_layout(
            polar=dict(
                radialaxis=dict(
                    visible=True,
                    range=[0, 1]
                )),
            title="Oligo Set Performance Comparison (Radar Chart)",
            showlegend=True,
            height=600
        )
        
        # Save interactive plot
        radar_file = self.interactive_dir / "performance_radar.html"
        fig.write_html(str(radar_file))
        
        return str(radar_file)
    
    def generate_all_comparative_plots(self, analyzer: ComparativeAnalyzer,
                                     differential_analyzer: Optional[DifferentialAnalyzer] = None) -> Dict[str, str]:
        """Generate all comparative visualization plots."""
        
        plots = {}
        
        # Main comparison dashboard
        try:
            plots['dashboard'] = self.create_comparison_dashboard(analyzer)
        except Exception as e:
            logging.error(f"Failed to create comparison dashboard: {e}")
        
        # Coverage distribution plots
        try:
            plots['coverage_distributions'] = self.create_coverage_distribution_plots(analyzer.oligo_sets)
        except Exception as e:
            logging.error(f"Failed to create coverage distribution plots: {e}")
        
        # Gap analysis comparison
        try:
            plots['gap_analysis'] = self.create_gap_analysis_comparison(analyzer.oligo_sets)
        except Exception as e:
            logging.error(f"Failed to create gap analysis plots: {e}")
        
        # Performance heatmap
        try:
            plots['performance_heatmap'] = self.create_performance_heatmap(analyzer)
        except Exception as e:
            logging.error(f"Failed to create performance heatmap: {e}")
        
        # Interactive ranking plot
        try:
            plots['ranking_radar'] = self.create_interactive_ranking_plot(analyzer)
        except Exception as e:
            logging.error(f"Failed to create ranking plot: {e}")
        
        # Statistical comparison (if differential analyzer provided)
        if differential_analyzer and len(analyzer.oligo_sets) >= 2:
            try:
                plots['statistical_comparison'] = self.create_statistical_comparison_plot(
                    differential_analyzer, analyzer.oligo_sets
                )
            except Exception as e:
                logging.error(f"Failed to create statistical comparison plot: {e}")
        
        logging.info(f"Generated {len(plots)} comparative visualization plots")
        return plots