#!/usr/bin/env python3

"""
coverage_viz.py

Enhanced visualization functions for oligo set coverage evaluation.
Generates comprehensive plots including heatmaps, distributions,
and gap analysis visualizations.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import warnings

# Suppress matplotlib warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')


class CoverageVisualizer:
    """Enhanced visualization generator for coverage analysis."""
    
    def __init__(
        self,
        coverage_stats: Dict[str, Any],
        gap_analysis: Dict[str, Any],
        output_dir: Path,
        plot_format: str = 'png',
        dpi: int = 300
    ):
        """
        Initialize the coverage visualizer.
        
        Args:
            coverage_stats: Coverage statistics from CoverageAnalyzer
            gap_analysis: Gap analysis results from GapAnalyzer
            output_dir: Directory to save plots
            plot_format: Output format for plots ('png', 'pdf', 'svg')
            dpi: Plot resolution in DPI
        """
        self.coverage_stats = coverage_stats
        self.gap_analysis = gap_analysis
        self.output_dir = Path(output_dir)
        self.plot_format = plot_format
        self.dpi = dpi
        
        # Create plots directory
        self.plots_dir = self.output_dir / "plots"
        self.plots_dir.mkdir(exist_ok=True)
        
        # Set up matplotlib style
        plt.style.use('seaborn-v0_8')
        sns.set_palette("husl")
    
    def generate_all_plots(self) -> None:
        """Generate all coverage visualization plots."""
        logging.info("Generating coverage overview plot...")
        self.plot_coverage_overview()
        
        logging.info("Generating depth distribution plots...")
        self.plot_depth_distribution()
        
        logging.info("Generating gap analysis plots...")
        self.plot_gap_analysis()
        
        logging.info("Generating coverage heatmap...")
        self.plot_coverage_heatmap()
        
        logging.info("Generating per-reference plots...")
        self.plot_per_reference_coverage()
        
        logging.info("All visualizations generated successfully")
    
    def plot_coverage_overview(self) -> None:
        """Generate overview plot with key coverage metrics."""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Coverage breadth pie chart
        breadth = self.coverage_stats.get('coverage_breadth', 0)
        uncovered = 100 - breadth
        
        ax1.pie([breadth, uncovered], 
                labels=[f'Covered ({breadth:.1f}%)', f'Uncovered ({uncovered:.1f}%)'],
                autopct='%1.1f%%',
                colors=['#2E8B57', '#DC143C'],
                startangle=90)
        ax1.set_title('Coverage Breadth', fontsize=14, fontweight='bold')
        
        # Depth distribution bar plot
        depth_dist = self.coverage_stats.get('depth_distribution', {})
        if depth_dist:
            thresholds = list(depth_dist.keys())
            percentages = list(depth_dist.values())
            
            bars = ax2.bar(range(len(thresholds)), percentages, 
                          color=plt.cm.viridis(np.linspace(0, 1, len(thresholds))))
            ax2.set_xlabel('Minimum Coverage Depth')
            ax2.set_ylabel('Percentage of Bases (%)')
            ax2.set_title('Coverage Depth Distribution', fontsize=14, fontweight='bold')
            ax2.set_xticks(range(len(thresholds)))
            ax2.set_xticklabels([f'≥{t}x' for t in thresholds])
            
            # Add value labels on bars
            for bar, pct in zip(bars, percentages):
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2., height + 1,
                        f'{pct:.1f}%', ha='center', va='bottom', fontsize=10)
        
        # Gap statistics
        gap_stats = [
            ('Total Gaps', self.gap_analysis.get('total_gaps', 0)),
            ('Mean Gap Size', f"{self.gap_analysis.get('mean_gap_size', 0):.0f} bp"),
            ('Largest Gap', f"{self.gap_analysis.get('max_gap_size', 0):,} bp"),
            ('Gap Percentage', f"{self.gap_analysis.get('gap_percentage', 0):.1f}%")
        ]
        
        ax3.axis('off')
        y_pos = 0.8
        for label, value in gap_stats:
            ax3.text(0.1, y_pos, f'{label}:', fontweight='bold', fontsize=12)
            ax3.text(0.6, y_pos, str(value), fontsize=12)
            y_pos -= 0.15
        
        ax3.set_title('Gap Analysis Summary', fontsize=14, fontweight='bold')
        ax3.set_xlim(0, 1)
        ax3.set_ylim(0, 1)
        
        # Coverage uniformity metrics
        uniformity_stats = [
            ('Mean Depth', f"{self.coverage_stats.get('mean_depth', 0):.1f}x"),
            ('Std Deviation', f"{self.coverage_stats.get('depth_std', 0):.1f}x"),
            ('CV', f"{self.coverage_stats.get('coverage_cv', 0):.2f}"),
            ('Uniformity Score', f"{self.coverage_stats.get('uniformity_score', 0):.2f}")
        ]
        
        ax4.axis('off')
        y_pos = 0.8
        for label, value in uniformity_stats:
            ax4.text(0.1, y_pos, f'{label}:', fontweight='bold', fontsize=12)
            ax4.text(0.6, y_pos, str(value), fontsize=12)
            y_pos -= 0.15
        
        ax4.set_title('Coverage Uniformity', fontsize=14, fontweight='bold')
        ax4.set_xlim(0, 1)
        ax4.set_ylim(0, 1)
        
        plt.tight_layout()
        self._save_plot(fig, 'coverage_overview')
    
    def plot_depth_distribution(self) -> None:
        """Generate detailed coverage depth distribution plots."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Histogram of coverage depths
        hist_data = self.coverage_stats.get('depth_histogram', {})
        if hist_data:
            counts = hist_data.get('counts', [])
            bin_edges = hist_data.get('bin_edges', [])
            
            if counts and bin_edges:
                bin_centers = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in range(len(counts))]
                ax1.bar(bin_centers, counts, width=np.diff(bin_edges), 
                       alpha=0.7, color='skyblue', edgecolor='black')
        
        ax1.set_xlabel('Coverage Depth')
        ax1.set_ylabel('Number of Positions')
        ax1.set_title('Coverage Depth Histogram', fontweight='bold')
        ax1.grid(True, alpha=0.3)
        
        # Cumulative distribution
        depth_dist = self.coverage_stats.get('depth_distribution', {})
        if depth_dist:
            thresholds = sorted(depth_dist.keys())
            cumulative = [depth_dist[t] for t in thresholds]
            
            ax2.plot(thresholds, cumulative, 'o-', linewidth=2, markersize=8)
            ax2.fill_between(thresholds, cumulative, alpha=0.3)
            ax2.set_xlabel('Minimum Coverage Depth')
            ax2.set_ylabel('Percentage of Bases ≥ Depth (%)')
            ax2.set_title('Cumulative Coverage Distribution', fontweight='bold')
            ax2.grid(True, alpha=0.3)
            
            # Add target coverage line
            target = self.coverage_stats['parameters'].get('target_coverage', 10)
            if target in depth_dist:
                ax2.axvline(target, color='red', linestyle='--', 
                           label=f'Target Coverage ({target}x)')
                ax2.legend()
        
        plt.tight_layout()
        self._save_plot(fig, 'depth_distribution')
    
    def plot_gap_analysis(self) -> None:
        """Generate gap analysis visualization plots."""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Gap size distribution
        size_dist = self.gap_analysis.get('size_distribution', {})
        if size_dist:
            ranges = list(size_dist.keys())
            counts = list(size_dist.values())
            
            bars = ax1.bar(range(len(ranges)), counts, color='coral')
            ax1.set_xlabel('Gap Size Range (bp)')
            ax1.set_ylabel('Number of Gaps')
            ax1.set_title('Gap Size Distribution', fontweight='bold')
            ax1.set_xticks(range(len(ranges)))
            ax1.set_xticklabels(ranges, rotation=45)
            
            # Add count labels on bars
            for bar, count in zip(bars, counts):
                height = bar.get_height()
                if height > 0:
                    ax1.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                            str(count), ha='center', va='bottom')
        
        # Largest gaps bar chart
        largest_gaps = self.gap_analysis.get('largest_gaps', [])[:10]
        if largest_gaps:
            gap_names = [f"{g['chromosome']}:{g['start']}-{g['end']}" 
                        for g in largest_gaps]
            gap_sizes = [g['length'] for g in largest_gaps]
            
            bars = ax2.barh(range(len(gap_names)), gap_sizes, color='lightcoral')
            ax2.set_xlabel('Gap Size (bp)')
            ax2.set_ylabel('Gap Location')
            ax2.set_title('Top 10 Largest Gaps', fontweight='bold')
            ax2.set_yticks(range(len(gap_names)))
            ax2.set_yticklabels([name.split(':')[0][:20] + '...' if len(name) > 20 else name 
                               for name in gap_names])
            
            # Add size labels
            for i, (bar, size) in enumerate(zip(bars, gap_sizes)):
                ax2.text(bar.get_width() + max(gap_sizes) * 0.01, bar.get_y() + bar.get_height()/2,
                        f'{size:,}', va='center', fontsize=10)
        
        # Gap density by chromosome
        per_ref_stats = self.coverage_stats.get('per_reference', {})
        if per_ref_stats:
            chromosomes = list(per_ref_stats.keys())[:15]  # Show top 15
            gap_counts = [per_ref_stats[chrom]['gaps'] for chrom in chromosomes]
            
            bars = ax3.bar(range(len(chromosomes)), gap_counts, color='lightblue')
            ax3.set_xlabel('Reference Sequence')
            ax3.set_ylabel('Number of Gaps')
            ax3.set_title('Gap Count by Reference', fontweight='bold')
            ax3.set_xticks(range(len(chromosomes)))
            ax3.set_xticklabels([chrom[:10] + '...' if len(chrom) > 10 else chrom 
                               for chrom in chromosomes], rotation=45)
        
        # Gap percentage vs reference length
        if per_ref_stats:
            ref_lengths = []
            gap_percentages = []
            
            for chrom, stats in per_ref_stats.items():
                ref_lengths.append(stats['length'])
                coverage_breadth = stats['coverage_breadth']
                gap_percentages.append(100 - coverage_breadth)
            
            ax4.scatter(ref_lengths, gap_percentages, alpha=0.6, s=50, color='orange')
            ax4.set_xlabel('Reference Length (bp)')
            ax4.set_ylabel('Gap Percentage (%)')
            ax4.set_title('Gap Percentage vs Reference Length', fontweight='bold')
            ax4.grid(True, alpha=0.3)
            
            if ref_lengths:
                ax4.set_xscale('log')
        
        plt.tight_layout()
        self._save_plot(fig, 'gap_analysis')
    
    def plot_coverage_heatmap(self) -> None:
        """Generate coverage heatmap for reference sequences."""
        per_ref_stats = self.coverage_stats.get('per_reference', {})
        
        if not per_ref_stats:
            logging.warning("No per-reference statistics available for heatmap")
            return
        
        # Prepare data for heatmap
        refs = list(per_ref_stats.keys())[:20]  # Show top 20 references
        
        metrics = ['coverage_breadth', 'mean_depth', 'max_depth', 'gaps']
        metric_labels = ['Coverage Breadth (%)', 'Mean Depth', 'Max Depth', 'Gap Count']
        
        heatmap_data = []
        for metric in metrics:
            row = [per_ref_stats[ref][metric] for ref in refs]
            
            # Normalize each metric to 0-1 scale for better visualization
            if metric == 'coverage_breadth':
                row = [x / 100.0 for x in row]  # Already percentage
            elif metric in ['mean_depth', 'max_depth']:
                max_val = max(row) if row else 1
                row = [x / max_val for x in row]
            else:  # gaps
                max_val = max(row) if row else 1
                row = [1 - (x / max_val) for x in row]  # Invert so fewer gaps = better
            
            heatmap_data.append(row)
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(max(12, len(refs) * 0.8), 8))
        
        heatmap_array = np.array(heatmap_data)
        
        # Custom colormap
        colors = ['#d32f2f', '#ff9800', '#ffc107', '#8bc34a', '#4caf50']
        n_bins = 100
        cmap = LinearSegmentedColormap.from_list('coverage', colors, N=n_bins)
        
        im = ax.imshow(heatmap_array, cmap=cmap, aspect='auto', vmin=0, vmax=1)
        
        # Set labels
        ax.set_xticks(range(len(refs)))
        ax.set_xticklabels([ref[:15] + '...' if len(ref) > 15 else ref for ref in refs], 
                          rotation=45, ha='right')
        ax.set_yticks(range(len(metric_labels)))
        ax.set_yticklabels(metric_labels)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Normalized Score (0=Poor, 1=Excellent)', rotation=270, labelpad=20)
        
        # Add value annotations
        for i in range(len(metric_labels)):
            for j in range(len(refs)):
                value = heatmap_array[i, j]
                color = 'white' if value < 0.5 else 'black'
                ax.text(j, i, f'{value:.2f}', ha='center', va='center', 
                       color=color, fontsize=8)
        
        ax.set_title('Coverage Quality Heatmap by Reference Sequence', 
                    fontweight='bold', pad=20)
        
        plt.tight_layout()
        self._save_plot(fig, 'coverage_heatmap')
    
    def plot_per_reference_coverage(self) -> None:
        """Generate per-reference coverage statistics plots."""
        per_ref_stats = self.coverage_stats.get('per_reference', {})
        
        if not per_ref_stats:
            logging.warning("No per-reference statistics available")
            return
        
        # Select top references by length
        refs_sorted = sorted(per_ref_stats.items(), 
                           key=lambda x: x[1]['length'], reverse=True)[:15]
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Coverage breadth by reference
        refs = [ref for ref, _ in refs_sorted]
        breadths = [stats['coverage_breadth'] for _, stats in refs_sorted]
        
        bars = ax1.barh(range(len(refs)), breadths, color='lightgreen')
        ax1.set_xlabel('Coverage Breadth (%)')
        ax1.set_ylabel('Reference Sequence')
        ax1.set_title('Coverage Breadth by Reference', fontweight='bold')
        ax1.set_yticks(range(len(refs)))
        ax1.set_yticklabels([ref[:20] + '...' if len(ref) > 20 else ref for ref in refs])
        
        # Add percentage labels
        for i, (bar, breadth) in enumerate(zip(bars, breadths)):
            ax1.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
                    f'{breadth:.1f}%', va='center', fontsize=10)
        
        # Mean depth by reference
        mean_depths = [stats['mean_depth'] for _, stats in refs_sorted]
        
        bars = ax2.barh(range(len(refs)), mean_depths, color='lightblue')
        ax2.set_xlabel('Mean Coverage Depth')
        ax2.set_ylabel('Reference Sequence')
        ax2.set_title('Mean Depth by Reference', fontweight='bold')
        ax2.set_yticks(range(len(refs)))
        ax2.set_yticklabels([ref[:20] + '...' if len(ref) > 20 else ref for ref in refs])
        
        # Reference length vs coverage breadth scatter
        lengths = [stats['length'] for _, stats in per_ref_stats.items()]
        all_breadths = [stats['coverage_breadth'] for _, stats in per_ref_stats.items()]
        
        ax3.scatter(lengths, all_breadths, alpha=0.6, s=50, color='purple')
        ax3.set_xlabel('Reference Length (bp)')
        ax3.set_ylabel('Coverage Breadth (%)')
        ax3.set_title('Coverage vs Reference Length', fontweight='bold')
        ax3.grid(True, alpha=0.3)
        
        if lengths:
            ax3.set_xscale('log')
        
        # Coverage efficiency (breadth vs mean depth)
        all_mean_depths = [stats['mean_depth'] for _, stats in per_ref_stats.items()]
        
        scatter = ax4.scatter(all_mean_depths, all_breadths, alpha=0.6, s=50, 
                             c=lengths, cmap='viridis')
        ax4.set_xlabel('Mean Coverage Depth')
        ax4.set_ylabel('Coverage Breadth (%)')
        ax4.set_title('Coverage Efficiency (colored by length)', fontweight='bold')
        ax4.grid(True, alpha=0.3)
        
        # Add colorbar for length
        cbar = plt.colorbar(scatter, ax=ax4)
        cbar.set_label('Reference Length (bp)')
        
        plt.tight_layout()
        self._save_plot(fig, 'per_reference_coverage')
    
    def _save_plot(self, fig, filename: str) -> None:
        """Save plot to file with proper formatting."""
        filepath = self.plots_dir / f"{filename}.{self.plot_format}"
        
        try:
            fig.savefig(filepath, dpi=self.dpi, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
            plt.close(fig)
            logging.debug(f"Plot saved: {filepath}")
            
        except Exception as e:
            logging.error(f"Error saving plot {filename}: {e}")
            plt.close(fig)