#!/usr/bin/env python3

"""
interactive_plots.py

Interactive plotting utilities using Plotly for enhanced coverage visualization.
Creates dynamic, zoomable, and interactive plots for comprehensive analysis.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
import numpy as np
import pandas as pd

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.offline as pyo


class InteractivePlotter:
    """Generate interactive plots for coverage analysis."""
    
    def __init__(
        self,
        coverage_stats: Dict[str, Any],
        gap_analysis: Dict[str, Any],
        quality_scores: Dict[str, Any],
        reference_analysis: Dict[str, Any],
        output_dir: Path
    ):
        """
        Initialize the interactive plotter.
        
        Args:
            coverage_stats: Coverage statistics
            gap_analysis: Gap analysis results
            quality_scores: Quality assessment results
            reference_analysis: Reference sequence analysis
            output_dir: Directory to save plots
        """
        self.coverage_stats = coverage_stats
        self.gap_analysis = gap_analysis
        self.quality_scores = quality_scores
        self.reference_analysis = reference_analysis
        self.output_dir = Path(output_dir)
        
        # Create plots directory
        self.plots_dir = self.output_dir / "interactive_plots"
        self.plots_dir.mkdir(exist_ok=True)
        
        # Plot configuration
        self.default_config = {
            'displayModeBar': True,
            'modeBarButtonsToRemove': ['pan2d', 'lasso2d'],
            'displaylogo': False,
            'toImageButtonOptions': {
                'format': 'png',
                'filename': 'baitUtils_plot',
                'height': 600,
                'width': 800,
                'scale': 2
            }
        }
    
    def create_coverage_dashboard(self) -> str:
        """Create comprehensive coverage analysis dashboard."""
        logging.info("Creating interactive coverage dashboard...")
        
        # Create subplot figure with multiple panels
        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=[
                'Coverage Breadth Distribution', 'Coverage Depth Histogram',
                'Gap Size Distribution', 'Quality Component Scores',
                'Reference Sequence Features', 'Coverage vs Reference Length'
            ],
            specs=[
                [{"type": "pie"}, {"type": "histogram"}],
                [{"type": "bar"}, {"type": "scatterpolar"}],
                [{"type": "heatmap"}, {"type": "scatter"}]
            ],
            vertical_spacing=0.08,
            horizontal_spacing=0.1
        )
        
        # 1. Coverage Breadth Pie Chart
        breadth = self.coverage_stats.get('coverage_breadth', 0)
        uncovered = 100 - breadth
        
        fig.add_trace(go.Pie(
            labels=['Covered', 'Uncovered'],
            values=[breadth, uncovered],
            hole=0.4,
            marker_colors=['#2E8B57', '#DC143C'],
            textinfo='label+percent+value',
            textfont_size=12,
            showlegend=False
        ), row=1, col=1)
        
        # 2. Coverage Depth Histogram
        depth_dist = self.coverage_stats.get('depth_distribution', {})
        if depth_dist:
            thresholds = list(depth_dist.keys())
            percentages = list(depth_dist.values())
            
            fig.add_trace(go.Bar(
                x=[f"≥{t}x" for t in thresholds],
                y=percentages,
                marker_color='rgba(55, 128, 191, 0.7)',
                text=[f'{p:.1f}%' for p in percentages],
                textposition='auto',
                showlegend=False
            ), row=1, col=2)
        
        # 3. Gap Size Distribution
        size_dist = self.gap_analysis.get('size_distribution', {})
        if size_dist:
            ranges = list(size_dist.keys())
            counts = list(size_dist.values())
            
            fig.add_trace(go.Bar(
                x=ranges,
                y=counts,
                marker_color='rgba(255, 99, 132, 0.7)',
                text=counts,
                textposition='auto',
                showlegend=False
            ), row=2, col=1)
        
        # 4. Quality Component Scores (Radar Chart)
        component_scores = self.quality_scores.get('component_scores', {})
        if component_scores:
            categories = [cat.replace('_', ' ').title() for cat in component_scores.keys()]
            values = list(component_scores.values())
            
            fig.add_trace(go.Scatterpolar(
                r=values,
                theta=categories,
                fill='toself',
                marker_color='rgb(32, 201, 151)',
                line_color='rgb(32, 201, 151)',
                showlegend=False
            ), row=2, col=2)
        
        # 5. Reference Sequence Features Heatmap
        ref_features = self.reference_analysis.get('sequence_features', {})
        if ref_features:
            # Prepare heatmap data
            feature_matrix, ref_names, feature_names = self._prepare_reference_heatmap_data(ref_features)
            
            fig.add_trace(go.Heatmap(
                z=feature_matrix,
                x=ref_names,
                y=feature_names,
                colorscale='Viridis',
                showscale=False,
                hoverongaps=False
            ), row=3, col=1)
        
        # 6. Coverage vs Reference Length
        per_ref_stats = self.coverage_stats.get('per_reference', {})
        if per_ref_stats:
            lengths = []
            breadths = []
            ref_ids = []
            
            for ref_id, stats in per_ref_stats.items():
                lengths.append(stats.get('length', 0))
                breadths.append(stats.get('coverage_breadth', 0))
                ref_ids.append(ref_id[:20] + '...' if len(ref_id) > 20 else ref_id)
            
            fig.add_trace(go.Scatter(
                x=lengths,
                y=breadths,
                mode='markers',
                marker=dict(
                    size=8,
                    color=breadths,
                    colorscale='RdYlGn',
                    showscale=False
                ),
                text=ref_ids,
                hovertemplate='<b>%{text}</b><br>' +
                             'Length: %{x:,} bp<br>' +
                             'Coverage: %{y:.1f}%<extra></extra>',
                showlegend=False
            ), row=3, col=2)
        
        # Update layout
        fig.update_layout(
            title_text="baitUtils Coverage Analysis Dashboard",
            title_font_size=20,
            height=1000,
            showlegend=False,
            template="plotly_white"
        )
        
        # Update axes labels
        fig.update_xaxes(title_text="Coverage Depth", row=1, col=2)
        fig.update_yaxes(title_text="Percentage of Bases (%)", row=1, col=2)
        
        fig.update_xaxes(title_text="Gap Size Range (bp)", row=2, col=1)
        fig.update_yaxes(title_text="Number of Gaps", row=2, col=1)
        
        fig.update_xaxes(title_text="Reference Sequences", row=3, col=1)
        fig.update_yaxes(title_text="Features", row=3, col=1)
        
        fig.update_xaxes(title_text="Reference Length (bp)", type="log", row=3, col=2)
        fig.update_yaxes(title_text="Coverage Breadth (%)", row=3, col=2)
        
        # Save dashboard
        dashboard_file = self.plots_dir / "coverage_dashboard.html"
        fig.write_html(str(dashboard_file), config=self.default_config)
        
        logging.info(f"Interactive dashboard saved: {dashboard_file}")
        return str(dashboard_file)
    
    def create_detailed_gap_analysis(self) -> str:
        """Create detailed interactive gap analysis plot."""
        logging.info("Creating detailed gap analysis plot...")
        
        # Create figure with secondary y-axis
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=[
                'Largest Coverage Gaps', 'Gap Features Correlation',
                'Gap Distribution by Chromosome', 'Gap Size vs Sequence Features'
            ],
            specs=[[{"type": "bar"}, {"type": "scatter"}],
                   [{"type": "bar"}, {"type": "scatter"}]]
        )
        
        # 1. Largest Coverage Gaps
        largest_gaps = self.gap_analysis.get('largest_gaps', [])[:15]
        if largest_gaps:
            gap_ids = [f"{g.get('chromosome', 'chr')}:{g.get('start', 0)//1000}k" 
                      for g in largest_gaps]
            gap_sizes = [g.get('length', 0) for g in largest_gaps]
            
            fig.add_trace(go.Bar(
                y=gap_ids,
                x=gap_sizes,
                orientation='h',
                marker_color='lightcoral',
                text=[f'{size:,} bp' for size in gap_sizes],
                textposition='auto',
                showlegend=False,
                hovertemplate='<b>%{y}</b><br>Size: %{x:,} bp<extra></extra>'
            ), row=1, col=1)
        
        # 2. Gap Features Correlation (if feature analysis available)
        gap_features = self.gap_analysis.get('gap_features', [])
        if gap_features:
            gc_contents = [f.get('gc_content', 50) for f in gap_features]
            gap_lengths = [f.get('length', 0) for f in gap_features]
            complexities = [f.get('complexity', 1.5) for f in gap_features]
            
            fig.add_trace(go.Scatter(
                x=gc_contents,
                y=gap_lengths,
                mode='markers',
                marker=dict(
                    size=8,
                    color=complexities,
                    colorscale='Plasma',
                    colorbar=dict(title="Complexity", x=0.48, len=0.4),
                    showscale=True
                ),
                showlegend=False,
                hovertemplate='GC: %{x:.1f}%<br>Length: %{y:,} bp<br>' +
                             'Complexity: %{marker.color:.2f}<extra></extra>'
            ), row=1, col=2)
        
        # 3. Gap Distribution by Chromosome
        per_ref_stats = self.coverage_stats.get('per_reference', {})
        if per_ref_stats:
            chromosomes = []
            gap_counts = []
            coverage_breadths = []
            
            for ref_id, stats in list(per_ref_stats.items())[:20]:  # Top 20
                chromosomes.append(ref_id[:15] + '...' if len(ref_id) > 15 else ref_id)
                gap_counts.append(stats.get('gaps', 0))
                coverage_breadths.append(stats.get('coverage_breadth', 0))
            
            fig.add_trace(go.Bar(
                x=chromosomes,
                y=gap_counts,
                marker_color=coverage_breadths,
                marker_colorscale='RdYlGn',
                showlegend=False,
                hovertemplate='<b>%{x}</b><br>Gaps: %{y}<br>' +
                             'Coverage: %{marker.color:.1f}%<extra></extra>'
            ), row=2, col=1)
        
        # 4. Gap Size vs Sequence Features (scatter matrix style)
        if gap_features:
            repeat_contents = [f.get('repeat_content', 0) for f in gap_features]
            n_contents = [f.get('n_content', 0) for f in gap_features]
            
            fig.add_trace(go.Scatter(
                x=repeat_contents,
                y=n_contents,
                mode='markers',
                marker=dict(
                    size=[min(50, max(5, length/100)) for length in gap_lengths],
                    color='orange',
                    opacity=0.6,
                    line=dict(width=1, color='DarkSlateGrey')
                ),
                showlegend=False,
                hovertemplate='Repeat Content: %{x:.1f}%<br>' +
                             'N Content: %{y:.1f}%<br>' +
                             'Gap Size: %{marker.size}<extra></extra>'
            ), row=2, col=2)
        
        # Update layout
        fig.update_layout(
            title_text="Detailed Gap Analysis",
            title_font_size=18,
            height=800,
            template="plotly_white"
        )
        
        # Update axes
        fig.update_xaxes(title_text="Gap Size (bp)", row=1, col=1)
        fig.update_yaxes(title_text="Gaps (Ranked)", row=1, col=1)
        
        fig.update_xaxes(title_text="GC Content (%)", row=1, col=2)
        fig.update_yaxes(title_text="Gap Length (bp)", type="log", row=1, col=2)
        
        fig.update_xaxes(title_text="Reference Sequence", row=2, col=1)
        fig.update_yaxes(title_text="Number of Gaps", row=2, col=1)
        
        fig.update_xaxes(title_text="Repeat Content (%)", row=2, col=2)
        fig.update_yaxes(title_text="N Content (%)", row=2, col=2)
        
        # Save plot
        gap_file = self.plots_dir / "detailed_gap_analysis.html"
        fig.write_html(str(gap_file), config=self.default_config)
        
        logging.info(f"Detailed gap analysis saved: {gap_file}")
        return str(gap_file)
    
    def create_quality_assessment_plot(self) -> str:
        """Create interactive quality assessment visualization."""
        logging.info("Creating quality assessment plot...")
        
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=[
                'Quality Score Breakdown', 'Performance vs Benchmarks',
                'Theoretical vs Actual Performance', 'Quality Trends'
            ],
            specs=[
                [{"type": "bar"}, {"type": "scatterpolar"}],
                [{"type": "bar"}, {"type": "scatter"}]
            ]
        )
        
        # 1. Quality Score Breakdown
        component_scores = self.quality_scores.get('component_scores', {})
        weighted_scores = self.quality_scores.get('weighted_scores', {})
        
        if component_scores and weighted_scores:
            components = [comp.replace('_', ' ').title() for comp in component_scores.keys()]
            raw_scores = list(component_scores.values())
            weighted_vals = list(weighted_scores.values())
            
            fig.add_trace(go.Bar(
                x=components,
                y=raw_scores,
                name='Raw Score',
                marker_color='lightblue',
                text=[f'{score:.2f}' for score in raw_scores],
                textposition='auto'
            ), row=1, col=1)
            
            fig.add_trace(go.Bar(
                x=components,
                y=weighted_vals,
                name='Weighted Score',
                marker_color='navy',
                text=[f'{score:.3f}' for score in weighted_vals],
                textposition='auto'
            ), row=1, col=1)
        
        # 2. Performance vs Benchmarks (Radar)
        benchmarks = self.quality_scores.get('benchmarks', {})
        if benchmarks:
            bench_labels = [label.replace('_', ' ').title() for label in benchmarks.keys()]
            bench_values = list(benchmarks.values())
            
            fig.add_trace(go.Scatterpolar(
                r=bench_values,
                theta=bench_labels,
                fill='toself',
                name='Actual Performance',
                marker_color='rgb(255, 165, 0)',
                line_color='rgb(255, 140, 0)'
            ), row=1, col=2)
            
            # Add benchmark line (1.0 = target)
            fig.add_trace(go.Scatterpolar(
                r=[1.0] * len(bench_labels),
                theta=bench_labels,
                mode='lines',
                name='Target Benchmark',
                line=dict(color='red', dash='dash')
            ), row=1, col=2)
        
        # 3. Theoretical vs Actual Performance
        theoretical = self.quality_scores.get('theoretical_comparison', {})
        if theoretical:
            current_perf = theoretical.get('current_performance', {})
            metrics = list(current_perf.keys())
            efficiencies = list(current_perf.values())
            
            fig.add_trace(go.Bar(
                x=metrics,
                y=efficiencies,
                marker_color=['green' if e >= 0.8 else 'orange' if e >= 0.5 else 'red' 
                             for e in efficiencies],
                text=[f'{e:.1%}' for e in efficiencies],
                textposition='auto',
                showlegend=False
            ), row=2, col=1)
        
        # 4. Quality Trends (simulated time series)
        # This would be actual historical data in a real implementation
        time_points = ['Initial', 'Current', 'Potential']
        quality_trend = [
            0.4,  # Hypothetical initial
            self.quality_scores.get('overall_score', 0.5),  # Current
            min(1.0, self.quality_scores.get('overall_score', 0.5) + 0.2)  # Potential improvement
        ]
        
        fig.add_trace(go.Scatter(
            x=time_points,
            y=quality_trend,
            mode='lines+markers',
            marker=dict(size=10, color='purple'),
            line=dict(width=3, color='purple'),
            showlegend=False,
            hovertemplate='%{x}: %{y:.3f}<extra></extra>'
        ), row=2, col=2)
        
        # Add quality thresholds
        thresholds = self.quality_scores.get('quality_thresholds', {})
        if thresholds:
            for category, threshold in thresholds.items():
                if category != 'poor':  # Skip poor threshold (0.0)
                    fig.add_hline(
                        y=threshold, 
                        line_dash="dash", 
                        line_color="gray",
                        row=2, col=2,
                        annotation_text=category.title()
                    )
        
        # Update layout
        fig.update_layout(
            title_text="Coverage Quality Assessment",
            title_font_size=18,
            height=800,
            template="plotly_white",
            showlegend=True
        )
        
        # Update axes
        fig.update_xaxes(title_text="Quality Components", row=1, col=1)
        fig.update_yaxes(title_text="Score", row=1, col=1)
        
        fig.update_xaxes(title_text="Performance Metrics", row=2, col=1)
        fig.update_yaxes(title_text="Efficiency", row=2, col=1)
        
        fig.update_xaxes(title_text="Analysis Stage", row=2, col=2)
        fig.update_yaxes(title_text="Quality Score", row=2, col=2)
        
        # Save plot
        quality_file = self.plots_dir / "quality_assessment.html"
        fig.write_html(str(quality_file), config=self.default_config)
        
        logging.info(f"Quality assessment plot saved: {quality_file}")
        return str(quality_file)
    
    def create_reference_analysis_plot(self) -> str:
        """Create interactive reference sequence analysis plot."""
        logging.info("Creating reference analysis plot...")
        
        ref_features = self.reference_analysis.get('sequence_features', {})
        if not ref_features:
            logging.warning("No reference analysis data available")
            return ""
        
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=[
                'GC Content Distribution', 'Complexity vs Repeat Content',
                'Challenging Regions', 'Feature Correlations'
            ]
        )
        
        # Extract feature data
        gc_contents = [features.get('gc_content', 50) for features in ref_features.values()]
        complexities = [features.get('shannon_entropy', 1.5) for features in ref_features.values()]
        repeat_contents = [features.get('repeat_content', 0) for features in ref_features.values()]
        ref_ids = list(ref_features.keys())
        
        # 1. GC Content Distribution
        fig.add_trace(go.Histogram(
            x=gc_contents,
            nbinsx=20,
            marker_color='lightgreen',
            opacity=0.7,
            showlegend=False,
            hovertemplate='GC Content: %{x:.1f}%<br>Count: %{y}<extra></extra>'
        ), row=1, col=1)
        
        # Add optimal range
        fig.add_vline(x=40, line_dash="dash", line_color="red", row=1, col=1, 
                      annotation_text="Optimal Range")
        fig.add_vline(x=60, line_dash="dash", line_color="red", row=1, col=1)
        
        # 2. Complexity vs Repeat Content
        fig.add_trace(go.Scatter(
            x=complexities,
            y=repeat_contents,
            mode='markers',
            marker=dict(
                size=8,
                color=gc_contents,
                colorscale='RdYlBu',
                colorbar=dict(title="GC Content (%)", x=0.48),
                showscale=True
            ),
            text=ref_ids,
            showlegend=False,
            hovertemplate='<b>%{text}</b><br>Complexity: %{x:.2f}<br>' +
                         'Repeat Content: %{y:.1f}%<br>GC: %{marker.color:.1f}%<extra></extra>'
        ), row=1, col=2)
        
        # 3. Challenging Regions
        challenging_regions = self.reference_analysis.get('challenging_regions', {})
        if challenging_regions:
            ref_names = []
            difficulty_scores = []
            difficulty_levels = []
            
            for ref_id, challenge_data in list(challenging_regions.items())[:20]:
                ref_names.append(ref_id[:15] + '...' if len(ref_id) > 15 else ref_id)
                difficulty_scores.append(challenge_data.get('challenging_score', 0))
                difficulty_levels.append(challenge_data.get('difficulty_level', 'Easy'))
            
            color_map = {'Easy': 'green', 'Moderate': 'yellow', 'Difficult': 'orange', 'Very Difficult': 'red'}
            colors = [color_map.get(level, 'gray') for level in difficulty_levels]
            
            fig.add_trace(go.Bar(
                x=ref_names,
                y=difficulty_scores,
                marker_color=colors,
                text=difficulty_levels,
                textposition='auto',
                showlegend=False,
                hovertemplate='<b>%{x}</b><br>Score: %{y}<br>Level: %{text}<extra></extra>'
            ), row=2, col=1)
        
        # 4. Feature Correlations Heatmap
        if len(ref_features) > 1:
            feature_names = ['GC Content', 'Complexity', 'Repeat Content', 'Max Homopolymer']
            correlation_matrix = np.corrcoef([
                gc_contents,
                complexities,
                repeat_contents,
                [features.get('max_homopolymer', 0) for features in ref_features.values()]
            ])
            
            fig.add_trace(go.Heatmap(
                z=correlation_matrix,
                x=feature_names,
                y=feature_names,
                colorscale='RdBu',
                zmid=0,
                showscale=False,
                hovertemplate='%{x} vs %{y}<br>Correlation: %{z:.2f}<extra></extra>'
            ), row=2, col=2)
        
        # Update layout
        fig.update_layout(
            title_text="Reference Sequence Analysis",
            title_font_size=18,
            height=800,
            template="plotly_white"
        )
        
        # Update axes
        fig.update_xaxes(title_text="GC Content (%)", row=1, col=1)
        fig.update_yaxes(title_text="Frequency", row=1, col=1)
        
        fig.update_xaxes(title_text="Shannon Entropy", row=1, col=2)
        fig.update_yaxes(title_text="Repeat Content (%)", row=1, col=2)
        
        fig.update_xaxes(title_text="Reference Sequences", row=2, col=1)
        fig.update_yaxes(title_text="Difficulty Score", row=2, col=1)
        
        # Save plot
        ref_file = self.plots_dir / "reference_analysis.html"
        fig.write_html(str(ref_file), config=self.default_config)
        
        logging.info(f"Reference analysis plot saved: {ref_file}")
        return str(ref_file)
    
    def _prepare_reference_heatmap_data(self, ref_features: Dict[str, Dict]) -> Tuple[np.ndarray, List[str], List[str]]:
        """Prepare data for reference features heatmap."""
        if not ref_features:
            return np.array([]), [], []
        
        # Select top references by length or take first 15
        selected_refs = list(ref_features.keys())[:15]
        
        # Feature names to include
        feature_names = ['GC Content', 'Complexity', 'Repeat Content', 'Homopolymers']
        feature_keys = ['gc_content', 'shannon_entropy', 'repeat_content', 'max_homopolymer']
        
        # Build matrix
        matrix = []
        for feature_key in feature_keys:
            row = []
            for ref_id in selected_refs:
                value = ref_features[ref_id].get(feature_key, 0)
                # Normalize values for better visualization
                if feature_key == 'gc_content':
                    value = value / 100  # 0-1 scale
                elif feature_key == 'shannon_entropy':
                    value = value / 2.0  # 0-1 scale (max entropy ~2 for DNA)
                elif feature_key == 'repeat_content':
                    value = min(value / 50, 1.0)  # Cap at 50%
                elif feature_key == 'max_homopolymer':
                    value = min(value / 20, 1.0)  # Cap at 20bp
                row.append(value)
            matrix.append(row)
        
        ref_names = [ref[:10] + '...' if len(ref) > 10 else ref for ref in selected_refs]
        
        return np.array(matrix), ref_names, feature_names
    
    def create_all_interactive_plots(self) -> Dict[str, str]:
        """Create all interactive plots and return file paths."""
        plot_files = {}
        
        try:
            plot_files['dashboard'] = self.create_coverage_dashboard()
        except Exception as e:
            logging.error(f"Error creating dashboard: {e}")
        
        try:
            plot_files['gap_analysis'] = self.create_detailed_gap_analysis()
        except Exception as e:
            logging.error(f"Error creating gap analysis: {e}")
        
        try:
            plot_files['quality_assessment'] = self.create_quality_assessment_plot()
        except Exception as e:
            logging.error(f"Error creating quality assessment: {e}")
        
        try:
            plot_files['reference_analysis'] = self.create_reference_analysis_plot()
        except Exception as e:
            logging.error(f"Error creating reference analysis: {e}")
        
        return plot_files