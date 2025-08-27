#!/usr/bin/env python3

"""
report_generator.py

Interactive HTML report generator for oligo set coverage evaluation.
Creates comprehensive, self-contained HTML reports with embedded interactive plots
and collapsible sections for detailed analysis.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
from datetime import datetime
import json
import base64
import io

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.offline as pyo


class InteractiveReportGenerator:
    """Generate comprehensive interactive HTML reports for coverage analysis."""
    
    def __init__(
        self,
        coverage_stats: Dict[str, Any],
        gap_analysis: Dict[str, Any],
        quality_scores: Dict[str, Any],
        reference_analysis: Dict[str, Any],
        output_dir: Path
    ):
        """
        Initialize the interactive report generator.
        
        Args:
            coverage_stats: Coverage statistics from CoverageAnalyzer
            gap_analysis: Gap analysis results from GapAnalyzer
            quality_scores: Quality scores from QualityScorer
            reference_analysis: Reference sequence analysis results
            output_dir: Directory to save the report
        """
        self.coverage_stats = coverage_stats
        self.gap_analysis = gap_analysis
        self.quality_scores = quality_scores
        self.reference_analysis = reference_analysis
        self.output_dir = Path(output_dir)
        
        # Report metadata
        self.report_title = "baitUtils Coverage Evaluation Report"
        self.generation_time = datetime.now().isoformat()
        
        # HTML template components
        self.html_components = []
    
    def generate_report(self) -> Path:
        """
        Generate the complete interactive HTML report.
        
        Returns:
            Path to the generated HTML report file
        """
        logging.info("Generating interactive HTML report...")
        
        # Generate all report sections
        self._create_header()
        self._create_executive_summary()
        self._create_coverage_overview()
        self._create_detailed_statistics()
        self._create_gap_analysis_section()
        self._create_quality_assessment()
        self._create_reference_analysis()
        self._create_recommendations()
        self._create_data_exports()
        self._create_methodology()
        
        # Assemble final HTML
        html_content = self._assemble_html()
        
        # Save report
        report_file = self.output_dir / "coverage_evaluation_report.html"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        logging.info(f"Interactive report generated: {report_file}")
        return report_file
    
    def _create_header(self) -> None:
        """Create report header with metadata."""
        header_html = f"""
        <div class="report-header">
            <h1>{self.report_title}</h1>
            <div class="metadata">
                <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p><strong>Analysis Parameters:</strong></p>
                <ul class="parameters">
                    <li>Minimum Identity: {self.coverage_stats.get('parameters', {}).get('min_identity', 'N/A')}%</li>
                    <li>Minimum Coverage: {self.coverage_stats.get('parameters', {}).get('min_coverage', 'N/A')}x</li>
                    <li>Target Coverage: {self.coverage_stats.get('parameters', {}).get('target_coverage', 'N/A')}x</li>
                </ul>
            </div>
        </div>
        """
        self.html_components.append(header_html)
    
    def _create_executive_summary(self) -> None:
        """Create executive summary with key metrics."""
        # Calculate summary metrics
        coverage_breadth = self.coverage_stats.get('coverage_breadth', 0)
        mean_depth = self.coverage_stats.get('mean_depth', 0)
        total_gaps = self.gap_analysis.get('total_gaps', 0)
        mapping_efficiency = self.coverage_stats.get('mapping_efficiency', 0)
        quality_score = self.quality_scores.get('overall_score', 0)
        
        # Determine quality category
        if quality_score >= 0.8:
            quality_category = "Excellent"
            quality_color = "#4CAF50"
        elif quality_score >= 0.6:
            quality_category = "Good"
            quality_color = "#FF9800"
        elif quality_score >= 0.4:
            quality_category = "Fair"
            quality_color = "#FF5722"
        else:
            quality_category = "Poor"
            quality_color = "#F44336"
        
        summary_html = f"""
        <div class="executive-summary">
            <h2>Executive Summary</h2>
            <div class="summary-grid">
                <div class="metric-card">
                    <div class="metric-value">{coverage_breadth:.1f}%</div>
                    <div class="metric-label">Coverage Breadth</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{mean_depth:.1f}x</div>
                    <div class="metric-label">Mean Depth</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{total_gaps:,}</div>
                    <div class="metric-label">Coverage Gaps</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{mapping_efficiency:.1f}%</div>
                    <div class="metric-label">Mapping Efficiency</div>
                </div>
                <div class="metric-card quality-card" style="border-color: {quality_color}">
                    <div class="metric-value" style="color: {quality_color}">{quality_category}</div>
                    <div class="metric-label">Overall Quality</div>
                    <div class="quality-score">Score: {quality_score:.2f}/1.00</div>
                </div>
            </div>
        </div>
        """
        self.html_components.append(summary_html)
    
    def _create_coverage_overview(self) -> None:
        """Create interactive coverage overview section."""
        # Generate interactive coverage overview plot
        coverage_plot = self._create_interactive_coverage_plot()
        depth_plot = self._create_interactive_depth_plot()
        
        overview_html = f"""
        <div class="section">
            <h2>Coverage Overview</h2>
            <div class="plot-container">
                {coverage_plot}
            </div>
            <div class="plot-container">
                {depth_plot}
            </div>
        </div>
        """
        self.html_components.append(overview_html)
    
    def _create_interactive_coverage_plot(self) -> str:
        """Create interactive coverage breadth visualization."""
        # Create pie chart for coverage breadth
        breadth = self.coverage_stats.get('coverage_breadth', 0)
        uncovered = 100 - breadth
        
        fig = go.Figure(data=[go.Pie(
            labels=['Covered', 'Uncovered'],
            values=[breadth, uncovered],
            hole=0.3,
            marker_colors=['#4CAF50', '#F44336'],
            textinfo='label+percent',
            textfont_size=14
        )])
        
        fig.update_layout(
            title="Reference Coverage Breadth",
            title_font_size=16,
            showlegend=True,
            height=400,
            margin=dict(t=50, b=50, l=50, r=50)
        )
        
        return fig.to_html(full_html=False, include_plotlyjs=False, div_id="coverage-pie")
    
    def _create_interactive_depth_plot(self) -> str:
        """Create interactive depth distribution plot."""
        depth_dist = self.coverage_stats.get('depth_distribution', {})
        
        if not depth_dist:
            return "<p>No depth distribution data available</p>"
        
        thresholds = list(depth_dist.keys())
        percentages = list(depth_dist.values())
        
        fig = go.Figure()
        
        # Bar plot
        fig.add_trace(go.Bar(
            x=[f"≥{t}x" for t in thresholds],
            y=percentages,
            name="Coverage Distribution",
            marker_color='rgba(55, 128, 191, 0.7)',
            text=[f'{p:.1f}%' for p in percentages],
            textposition='auto'
        ))
        
        fig.update_layout(
            title="Coverage Depth Distribution",
            xaxis_title="Minimum Coverage Depth",
            yaxis_title="Percentage of Bases (%)",
            title_font_size=16,
            height=400,
            showlegend=False
        )
        
        return fig.to_html(full_html=False, include_plotlyjs=False, div_id="depth-dist")
    
    def _create_detailed_statistics(self) -> None:
        """Create collapsible detailed statistics section."""
        # Create comprehensive statistics table
        stats_table = self._create_statistics_table()
        
        detailed_html = f"""
        <div class="section">
            <h2>Detailed Statistics</h2>
            <div class="collapsible">
                <button class="collapsible-header">View Detailed Statistics</button>
                <div class="collapsible-content">
                    {stats_table}
                </div>
            </div>
        </div>
        """
        self.html_components.append(detailed_html)
    
    def _create_statistics_table(self) -> str:
        """Create comprehensive statistics table."""
        stats_data = []
        
        # Coverage statistics
        stats_data.extend([
            ("Coverage Breadth", f"{self.coverage_stats.get('coverage_breadth', 0):.1f}%"),
            ("Mean Coverage Depth", f"{self.coverage_stats.get('mean_depth', 0):.1f}x"),
            ("Median Coverage Depth", f"{self.coverage_stats.get('median_depth', 0):.1f}x"),
            ("Coverage Std Dev", f"{self.coverage_stats.get('depth_std', 0):.1f}x"),
            ("Maximum Depth", f"{self.coverage_stats.get('max_depth', 0)}x"),
            ("Coverage CV", f"{self.coverage_stats.get('coverage_cv', 0):.2f}"),
            ("Coverage Gini", f"{self.coverage_stats.get('coverage_gini', 0):.2f}"),
            ("Uniformity Score", f"{self.coverage_stats.get('uniformity_score', 0):.2f}")
        ])
        
        # Mapping statistics
        stats_data.extend([
            ("Total Oligos", f"{self.coverage_stats.get('total_oligos', 0):,}"),
            ("Mapped Oligos", f"{self.coverage_stats.get('mapped_oligos', 0):,}"),
            ("Mapping Efficiency", f"{self.coverage_stats.get('mapping_efficiency', 0):.1f}%"),
            ("Total Mappings", f"{self.coverage_stats.get('total_mappings', 0):,}")
        ])
        
        # Reference statistics
        stats_data.extend([
            ("Reference Sequences", f"{self.coverage_stats.get('reference_sequences', 0):,}"),
            ("Total Reference Length", f"{self.coverage_stats.get('reference_length', 0):,} bp"),
            ("Covered Bases", f"{self.coverage_stats.get('covered_bases', 0):,} bp"),
            ("Uncovered Bases", f"{self.coverage_stats.get('uncovered_bases', 0):,} bp")
        ])
        
        # Gap statistics
        stats_data.extend([
            ("Total Gaps", f"{self.gap_analysis.get('total_gaps', 0):,}"),
            ("Gap Percentage", f"{self.gap_analysis.get('gap_percentage', 0):.1f}%"),
            ("Mean Gap Size", f"{self.gap_analysis.get('mean_gap_size', 0):.0f} bp"),
            ("Median Gap Size", f"{self.gap_analysis.get('median_gap_size', 0):.0f} bp"),
            ("Largest Gap", f"{self.gap_analysis.get('max_gap_size', 0):,} bp")
        ])
        
        # Create HTML table
        table_rows = []
        for metric, value in stats_data:
            table_rows.append(f"<tr><td>{metric}</td><td>{value}</td></tr>")
        
        table_html = f"""
        <table class="stats-table">
            <thead>
                <tr><th>Metric</th><th>Value</th></tr>
            </thead>
            <tbody>
                {''.join(table_rows)}
            </tbody>
        </table>
        """
        
        return table_html
    
    def _create_gap_analysis_section(self) -> None:
        """Create interactive gap analysis section."""
        gap_plot = self._create_interactive_gap_plot()
        largest_gaps_table = self._create_largest_gaps_table()
        
        gap_html = f"""
        <div class="section">
            <h2>Gap Analysis</h2>
            <div class="plot-container">
                {gap_plot}
            </div>
            <div class="collapsible">
                <button class="collapsible-header">Largest Coverage Gaps</button>
                <div class="collapsible-content">
                    {largest_gaps_table}
                </div>
            </div>
        </div>
        """
        self.html_components.append(gap_html)
    
    def _create_interactive_gap_plot(self) -> str:
        """Create interactive gap size distribution plot."""
        size_dist = self.gap_analysis.get('size_distribution', {})
        
        if not size_dist:
            return "<p>No gap distribution data available</p>"
        
        ranges = list(size_dist.keys())
        counts = list(size_dist.values())
        
        fig = go.Figure()
        
        fig.add_trace(go.Bar(
            x=ranges,
            y=counts,
            marker_color='rgba(255, 99, 132, 0.7)',
            text=counts,
            textposition='auto'
        ))
        
        fig.update_layout(
            title="Gap Size Distribution",
            xaxis_title="Gap Size Range (bp)",
            yaxis_title="Number of Gaps",
            title_font_size=16,
            height=400,
            showlegend=False
        )
        
        return fig.to_html(full_html=False, include_plotlyjs=False, div_id="gap-dist")
    
    def _create_largest_gaps_table(self) -> str:
        """Create table of largest coverage gaps."""
        largest_gaps = self.gap_analysis.get('largest_gaps', [])[:10]
        
        if not largest_gaps:
            return "<p>No gap data available</p>"
        
        table_rows = []
        for i, gap in enumerate(largest_gaps, 1):
            table_rows.append(
                f"<tr>"
                f"<td>{i}</td>"
                f"<td>{gap.get('chromosome', 'Unknown')}</td>"
                f"<td>{gap.get('start', 0):,}</td>"
                f"<td>{gap.get('end', 0):,}</td>"
                f"<td>{gap.get('length', 0):,}</td>"
                f"</tr>"
            )
        
        table_html = f"""
        <table class="gaps-table">
            <thead>
                <tr>
                    <th>Rank</th>
                    <th>Chromosome</th>
                    <th>Start</th>
                    <th>End</th>
                    <th>Length (bp)</th>
                </tr>
            </thead>
            <tbody>
                {''.join(table_rows)}
            </tbody>
        </table>
        """
        
        return table_html
    
    def _create_quality_assessment(self) -> None:
        """Create quality assessment section with scoring breakdown."""
        quality_breakdown = self._create_quality_breakdown()
        
        quality_html = f"""
        <div class="section">
            <h2>Quality Assessment</h2>
            {quality_breakdown}
        </div>
        """
        self.html_components.append(quality_html)
    
    def _create_quality_breakdown(self) -> str:
        """Create quality score breakdown visualization."""
        scores = self.quality_scores.get('component_scores', {})
        
        if not scores:
            return "<p>No quality scores available</p>"
        
        # Create radar chart for quality components
        categories = list(scores.keys())
        values = list(scores.values())
        
        fig = go.Figure()
        
        fig.add_trace(go.Scatterpolar(
            r=values,
            theta=categories,
            fill='toself',
            name='Quality Scores',
            line_color='rgb(32, 201, 151)'
        ))
        
        fig.update_layout(
            polar=dict(
                radialaxis=dict(
                    visible=True,
                    range=[0, 1]
                )),
            title="Quality Score Breakdown",
            title_font_size=16,
            height=500,
            showlegend=False
        )
        
        return fig.to_html(full_html=False, include_plotlyjs=False, div_id="quality-radar")
    
    def _create_reference_analysis(self) -> None:
        """Create reference sequence analysis section."""
        ref_analysis = self._create_reference_features_plot()
        
        ref_html = f"""
        <div class="section">
            <h2>Reference Sequence Analysis</h2>
            <div class="collapsible">
                <button class="collapsible-header">Reference Sequence Features</button>
                <div class="collapsible-content">
                    {ref_analysis}
                </div>
            </div>
        </div>
        """
        self.html_components.append(ref_html)
    
    def _create_reference_features_plot(self) -> str:
        """Create reference sequence features visualization."""
        features = self.reference_analysis.get('sequence_features', {})
        
        if not features:
            return "<p>No reference analysis data available</p>"
        
        # Create heatmap of reference features
        feature_data = []
        for ref_id, ref_features in features.items():
            feature_data.append({
                'Reference': ref_id[:20] + '...' if len(ref_id) > 20 else ref_id,
                'GC Content': ref_features.get('gc_content', 0),
                'Complexity': ref_features.get('complexity', 0),
                'Repeat Content': ref_features.get('repeat_content', 0),
                'Max Homopolymer': min(ref_features.get('max_homopolymer', 0), 20)  # Cap for viz
            })
        
        if not feature_data:
            return "<p>No feature data to display</p>"
        
        df = pd.DataFrame(feature_data)
        
        fig = go.Figure(data=go.Heatmap(
            z=df.iloc[:, 1:].values.T,
            x=df['Reference'],
            y=df.columns[1:],
            colorscale='Viridis',
            showscale=True
        ))
        
        fig.update_layout(
            title="Reference Sequence Features",
            title_font_size=16,
            height=400,
            xaxis_title="Reference Sequences",
            yaxis_title="Features"
        )
        
        return fig.to_html(full_html=False, include_plotlyjs=False, div_id="ref-features")
    
    def _create_recommendations(self) -> None:
        """Create recommendations section."""
        recommendations = self.gap_analysis.get('suggestions', [])
        quality_recommendations = self.quality_scores.get('recommendations', [])
        
        all_recommendations = recommendations + quality_recommendations
        
        if not all_recommendations:
            all_recommendations = ["Coverage analysis complete. No specific recommendations at this time."]
        
        rec_items = []
        for i, rec in enumerate(all_recommendations, 1):
            rec_items.append(f"<li>{rec}</li>")
        
        rec_html = f"""
        <div class="section">
            <h2>Recommendations</h2>
            <div class="recommendations">
                <ol>
                    {''.join(rec_items)}
                </ol>
            </div>
        </div>
        """
        self.html_components.append(rec_html)
    
    def _create_data_exports(self) -> None:
        """Create data exports section."""
        exports_html = """
        <div class="section">
            <h2>Data Exports</h2>
            <div class="exports-grid">
                <div class="export-item">
                    <h4>Coverage Data</h4>
                    <p>Per-position coverage depths</p>
                    <a href="data/coverage_per_position.csv" class="download-btn">Download CSV</a>
                </div>
                <div class="export-item">
                    <h4>Gap Regions</h4>
                    <p>Coverage gaps in BED format</p>
                    <a href="data/gap_regions.bed" class="download-btn">Download BED</a>
                </div>
                <div class="export-item">
                    <h4>Statistics Summary</h4>
                    <p>Complete statistics in JSON format</p>
                    <a href="data/statistics_summary.json" class="download-btn">Download JSON</a>
                </div>
            </div>
        </div>
        """
        self.html_components.append(exports_html)
    
    def _create_methodology(self) -> None:
        """Create methodology section."""
        method_html = """
        <div class="section">
            <h2>Methodology</h2>
            <div class="collapsible">
                <button class="collapsible-header">Analysis Methodology</button>
                <div class="collapsible-content">
                    <h4>Mapping</h4>
                    <p>Oligos were mapped to reference sequences using pblat with specified identity and length thresholds.</p>
                    
                    <h4>Coverage Calculation</h4>
                    <p>Coverage depth was calculated at each reference position. Coverage breadth represents the percentage of reference positions with coverage ≥ minimum threshold.</p>
                    
                    <h4>Gap Analysis</h4>
                    <p>Coverage gaps were identified as contiguous regions with coverage below the minimum threshold. Gap features were analyzed including GC content, sequence complexity, and repetitive elements.</p>
                    
                    <h4>Quality Scoring</h4>
                    <p>Overall quality score combines coverage breadth, depth uniformity, mapping efficiency, and gap characteristics using weighted averaging.</p>
                </div>
            </div>
        </div>
        """
        self.html_components.append(method_html)
    
    def _assemble_html(self) -> str:
        """Assemble final HTML document."""
        css_styles = self._get_css_styles()
        js_scripts = self._get_js_scripts()
        
        html_content = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>{self.report_title}</title>
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            <style>{css_styles}</style>
        </head>
        <body>
            <div class="container">
                {''.join(self.html_components)}
            </div>
            <script>{js_scripts}</script>
        </body>
        </html>
        """
        
        return html_content
    
    def _get_css_styles(self) -> str:
        """Get CSS styles for the report."""
        return """
        * { box-sizing: border-box; margin: 0; padding: 0; }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background-color: #f5f5f5;
        }
        
        .container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: white;
            box-shadow: 0 0 20px rgba(0,0,0,0.1);
        }
        
        .report-header {
            text-align: center;
            margin-bottom: 30px;
            padding: 30px 0;
            border-bottom: 3px solid #4CAF50;
        }
        
        .report-header h1 {
            color: #2c3e50;
            font-size: 2.5rem;
            margin-bottom: 15px;
        }
        
        .metadata {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            text-align: left;
            display: inline-block;
        }
        
        .parameters {
            list-style: none;
            margin-top: 10px;
        }
        
        .parameters li {
            padding: 2px 0;
        }
        
        .executive-summary {
            margin-bottom: 30px;
        }
        
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }
        
        .metric-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 25px;
            border-radius: 10px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        }
        
        .metric-card.quality-card {
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
            border: 3px solid;
        }
        
        .metric-value {
            font-size: 2.2rem;
            font-weight: bold;
            margin-bottom: 8px;
        }
        
        .metric-label {
            font-size: 1rem;
            opacity: 0.9;
        }
        
        .quality-score {
            font-size: 0.9rem;
            margin-top: 5px;
            opacity: 0.8;
        }
        
        .section {
            margin-bottom: 40px;
            padding: 20px;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.05);
        }
        
        .section h2 {
            color: #2c3e50;
            border-bottom: 2px solid #4CAF50;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }
        
        .plot-container {
            margin: 20px 0;
            padding: 15px;
            background: #fafafa;
            border-radius: 5px;
        }
        
        .collapsible {
            margin: 20px 0;
        }
        
        .collapsible-header {
            background-color: #4CAF50;
            color: white;
            cursor: pointer;
            padding: 15px;
            width: 100%;
            border: none;
            text-align: left;
            outline: none;
            font-size: 16px;
            border-radius: 5px;
            transition: background-color 0.3s;
        }
        
        .collapsible-header:hover {
            background-color: #45a049;
        }
        
        .collapsible-header:after {
            content: '+';
            float: right;
            font-weight: bold;
        }
        
        .collapsible-header.active:after {
            content: '-';
        }
        
        .collapsible-content {
            padding: 0;
            max-height: 0;
            overflow: hidden;
            transition: max-height 0.3s ease-out;
            background-color: #f1f1f1;
        }
        
        .collapsible-content.active {
            padding: 15px;
            max-height: 1000px;
        }
        
        .stats-table, .gaps-table {
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }
        
        .stats-table th, .stats-table td,
        .gaps-table th, .gaps-table td {
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }
        
        .stats-table th, .gaps-table th {
            background-color: #f2f2f2;
            font-weight: bold;
        }
        
        .stats-table tr:hover, .gaps-table tr:hover {
            background-color: #f5f5f5;
        }
        
        .recommendations ol {
            padding-left: 20px;
        }
        
        .recommendations li {
            margin: 10px 0;
            padding: 10px;
            background: #f8f9fa;
            border-left: 4px solid #4CAF50;
        }
        
        .exports-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }
        
        .export-item {
            padding: 20px;
            background: #f8f9fa;
            border-radius: 8px;
            text-align: center;
        }
        
        .download-btn {
            display: inline-block;
            padding: 10px 20px;
            background: #4CAF50;
            color: white;
            text-decoration: none;
            border-radius: 5px;
            margin-top: 10px;
            transition: background-color 0.3s;
        }
        
        .download-btn:hover {
            background: #45a049;
        }
        
        @media (max-width: 768px) {
            .summary-grid {
                grid-template-columns: 1fr;
            }
            
            .container {
                padding: 10px;
            }
            
            .report-header h1 {
                font-size: 2rem;
            }
        }
        """
    
    def _get_js_scripts(self) -> str:
        """Get JavaScript for interactive functionality."""
        return """
        // Collapsible sections
        document.querySelectorAll('.collapsible-header').forEach(function(header) {
            header.addEventListener('click', function() {
                this.classList.toggle('active');
                var content = this.nextElementSibling;
                content.classList.toggle('active');
                
                if (content.classList.contains('active')) {
                    content.style.maxHeight = content.scrollHeight + 'px';
                } else {
                    content.style.maxHeight = '0';
                }
            });
        });
        
        // Smooth scrolling for internal links
        document.querySelectorAll('a[href^="#"]').forEach(anchor => {
            anchor.addEventListener('click', function (e) {
                e.preventDefault();
                document.querySelector(this.getAttribute('href')).scrollIntoView({
                    behavior: 'smooth'
                });
            });
        });
        """
    
    def export_data_files(self) -> None:
        """Export data files referenced in the report."""
        data_dir = self.output_dir / "data"
        data_dir.mkdir(exist_ok=True)
        
        # Export statistics summary as JSON
        stats_summary = {
            'coverage_stats': self.coverage_stats,
            'gap_analysis': self.gap_analysis,
            'quality_scores': self.quality_scores,
            'reference_analysis': self.reference_analysis,
            'generation_time': self.generation_time
        }
        
        stats_file = data_dir / "statistics_summary.json"
        with open(stats_file, 'w') as f:
            json.dump(stats_summary, f, indent=2, default=str)
        
        logging.info(f"Data files exported to {data_dir}")