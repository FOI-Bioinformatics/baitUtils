#!/usr/bin/env python3

"""
comparative_report_generator.py

Comprehensive reporting system for comparative oligo set analysis.
Generates detailed HTML reports comparing multiple oligo sets with
statistical analysis, visualizations, and actionable recommendations.

This module provides the final reporting layer for Phase 3 comparative
analysis functionality.
"""

import logging
from datetime import datetime
from typing import Dict, Optional
from pathlib import Path
import base64

from baitUtils.comparative_analyzer import ComparativeAnalyzer
from baitUtils.differential_analysis import DifferentialAnalyzer
from baitUtils.templates import load_template
from baitUtils.comparative_visualizations import ComparativeVisualizer


class ComparativeReportGenerator:
    """
    Generate comprehensive HTML reports for comparative oligo set analysis.
    
    Creates professional, interactive reports with embedded visualizations,
    statistical analysis results, and detailed recommendations for oligo set
    selection and optimization.
    """
    
    def __init__(self, analyzer: ComparativeAnalyzer, output_dir: Path,
                 differential_analyzer: Optional[DifferentialAnalyzer] = None):
        """
        Initialize comparative report generator.
        
        Args:
            analyzer: Comparative analyzer with oligo set results
            output_dir: Directory for output files
            differential_analyzer: Optional statistical analysis
        """
        self.analyzer = analyzer
        self.differential_analyzer = differential_analyzer
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.report_file = self.output_dir / "comparative_analysis_report.html"
        self.visualizer = ComparativeVisualizer(output_dir)
        
    def generate_report(self) -> str:
        """Generate comprehensive comparative analysis report."""
        
        logging.info("Generating comparative analysis report...")
        
        # Generate all visualizations
        plots = self.visualizer.generate_all_comparative_plots(
            self.analyzer, self.differential_analyzer
        )
        
        # Generate HTML content
        html_content = self._generate_html_content(plots)
        
        # Write report
        with open(self.report_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        logging.info(f"Comparative analysis report saved to {self.report_file}")
        return str(self.report_file)
    
    def _generate_html_content(self, plots: Dict[str, str]) -> str:
        """Generate complete HTML report content."""
        
        # Get comparison data
        comparison_matrix = self.analyzer.generate_comparison_matrix()
        ranking = self.analyzer.generate_ranking()
        gap_analysis = self.analyzer.analyze_gap_overlap()
        
        html_parts = [
            self._get_html_header(),
            self._generate_executive_summary(),
            self._generate_comparison_overview(comparison_matrix),
            self._generate_performance_ranking(ranking),
            self._generate_detailed_analysis(),
            self._generate_statistical_analysis() if self.differential_analyzer else "",
            self._generate_gap_overlap_analysis(gap_analysis),
            self._generate_recommendations(),
            self._generate_visualizations_section(plots),
            self._generate_methodology_section(),
            self._get_html_footer()
        ]
        
        return '\n'.join(filter(None, html_parts))
    
    def _get_html_header(self) -> str:
        """Generate HTML header with CSS styling."""
        return self._html_header_template().replace("__CSS__", load_template("comparative_report.css"))
    
    def _html_header_template(self) -> str:
        return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Comparative Oligo Set Analysis Report</title>
    <style>__CSS__</style>
</head>
<body>
<div class="container">
    <div class="header">
        <h1>Comparative Oligo Set Analysis Report</h1>
        <div class="subtitle">Comprehensive comparison of multiple oligo set designs</div>
        <div class="subtitle">Generated on """ + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + """</div>
    </div>
"""
    
    def _generate_executive_summary(self) -> str:
        """Generate executive summary section."""
        
        # Identify best performer
        best_performer = self.analyzer.identify_best_performer('quality_score')
        ranking = self.analyzer.generate_ranking()
        
        n_sets = len(self.analyzer.oligo_sets)
        
        summary = f"""
    <div class="section">
        <h2>📊 Executive Summary</h2>
        
        <div class="metric-card">
            <h4>Oligo Sets Compared</h4>
            <div class="value">{n_sets}</div>
        </div>
        
        <div class="metric-card">
            <h4>Best Performer</h4>
            <div class="value">{best_performer.name}</div>
        </div>
        
        <div class="metric-card">
            <h4>Top Quality Score</h4>
            <div class="value">{best_performer.quality_score.overall_score:.1f}/10</div>
        </div>
        
        <div class="metric-card">
            <h4>Top Quality Grade</h4>
            <div class="value grade-{best_performer.quality_score.category.value.lower()}">{best_performer.quality_score.category.value}</div>
        </div>
        
        <h3>🎯 Key Findings</h3>
        <ul>
            <li><strong>{best_performer.name}</strong> achieved the highest overall quality score ({best_performer.quality_score.overall_score:.1f}/10, Grade {best_performer.quality_score.category.value})</li>
            <li>Coverage breadth ranges from {min(result.coverage_stats.get('coverage_breadth', 0) for result in self.analyzer.oligo_sets):.1f}% to {max(result.coverage_stats.get('coverage_breadth', 0) for result in self.analyzer.oligo_sets):.1f}%</li>
            <li>Gap counts vary from {min(result.gap_analysis.get('total_gaps', 0) for result in self.analyzer.oligo_sets):,} to {max(result.gap_analysis.get('total_gaps', 0) for result in self.analyzer.oligo_sets):,} gaps</li>
            <li>Mapping efficiency ranges from {min(result.coverage_stats.get('mapping_efficiency', 0) for result in self.analyzer.oligo_sets):.1f}% to {max(result.coverage_stats.get('mapping_efficiency', 0) for result in self.analyzer.oligo_sets):.1f}%</li>
        </ul>
        
        <h3>📈 Performance Ranking</h3>
        <ol>
        """
        
        for i, (name, score) in enumerate(ranking[:5], 1):  # Top 5
            oligo_set = next(result for result in self.analyzer.oligo_sets if result.name == name)
            grade_class = f"grade-{oligo_set.quality_score.category.value.lower()}"
            summary += f'<li><strong>{name}</strong> (Score: {score:.2f}, Grade: <span class="{grade_class}">{oligo_set.quality_score.category.value}</span>)</li>'
        
        if len(ranking) > 5:
            summary += f'<li><em>...and {len(ranking) - 5} more</em></li>'
        
        summary += """
        </ol>
    </div>
        """
        
        return summary
    
    def _generate_comparison_overview(self, comparison_matrix) -> str:
        """Generate comparison overview table."""
        
        if comparison_matrix is None:
            return ""
        
        # Generate HTML table
        html = """
    <div class="section">
        <h2>📋 Comparison Overview</h2>
        <p>Detailed metrics comparison across all oligo sets:</p>
        
        <table class="comparison-table">
            <thead>
                <tr>
                    <th>Oligo Set</th>
                    <th>Quality Score</th>
                    <th>Grade</th>
                    <th>Coverage Breadth (%)</th>
                    <th>Mean Depth (x)</th>
                    <th>Total Gaps</th>
                    <th>Mapping Efficiency (%)</th>
                    <th>Gini Coefficient</th>
                </tr>
            </thead>
            <tbody>
        """
        
        # Find best performers for highlighting
        best_quality = comparison_matrix['Quality_Score'].max()
        
        for _, row in comparison_matrix.iterrows():
            # Highlight best performers
            row_class = ""
            if row['Quality_Score'] == best_quality:
                row_class = "best-performer"
            
            grade_class = f"grade-{row['Quality_Grade'].lower()}"
            
            html += f"""
                <tr class="{row_class}">
                    <td><strong>{row['Name']}</strong></td>
                    <td>{row['Quality_Score']:.1f}/10</td>
                    <td><span class="{grade_class}">{row['Quality_Grade']}</span></td>
                    <td>{row['Coverage_Breadth_%']:.1f}%</td>
                    <td>{row['Mean_Depth_x']:.1f}x</td>
                    <td>{row['Total_Gaps']:,}</td>
                    <td>{row['Mapping_Efficiency_%']:.1f}%</td>
                    <td>{row['Gini_Coefficient']:.3f}</td>
                </tr>
            """
        
        html += """
            </tbody>
        </table>
    </div>
        """
        
        return html
    
    def _generate_performance_ranking(self, ranking) -> str:
        """Generate performance ranking section."""
        
        html = """
    <div class="section">
        <h2>🏆 Performance Ranking</h2>
        <p>Oligo sets ranked by composite performance score:</p>
        
        <div style="display: flex; flex-wrap: wrap; justify-content: space-around;">
        """
        
        for i, (name, score) in enumerate(ranking, 1):
            # Find corresponding oligo set for quality grade
            oligo_set = next(result for result in self.analyzer.oligo_sets if result.name == name)
            grade_class = f"grade-{oligo_set.quality_score.category.value.lower()}"
            
            # Color based on ranking
            if i == 1:
                bg_color = "linear-gradient(135deg, #ffd700, #ffed4a)"  # Gold
            elif i == 2:
                bg_color = "linear-gradient(135deg, #c0c0c0, #e2e8f0)"  # Silver
            elif i == 3:
                bg_color = "linear-gradient(135deg, #cd7f32, #d69e2e)"  # Bronze
            else:
                bg_color = "linear-gradient(135deg, #6c757d, #adb5bd)"  # Gray
            
            html += f"""
            <div class="metric-card" style="background: {bg_color}; color: #333;">
                <h4>#{i} {name}</h4>
                <div class="value">{score:.2f}</div>
                <div>Grade: <span class="{grade_class}">{oligo_set.quality_score.category.value}</span></div>
            </div>
            """
        
        html += """
        </div>
    </div>
        """
        
        return html
    
    def _generate_detailed_analysis(self) -> str:
        """Generate detailed analysis section."""
        
        pairwise_comparisons = self.analyzer.calculate_pairwise_comparisons()
        
        html = """
    <div class="section">
        <h2>🔍 Detailed Analysis</h2>
        
        <button class="collapsible">Pairwise Comparisons</button>
        <div class="content">
            <h3>Pairwise Comparison Metrics</h3>
            <p>Detailed comparison between all pairs of oligo sets:</p>
            
            <table class="comparison-table">
                <thead>
                    <tr>
                        <th>Comparison</th>
                        <th>Coverage Breadth Diff (%)</th>
                        <th>Mean Depth Diff (x)</th>
                        <th>Gap Count Diff</th>
                        <th>Quality Score Diff</th>
                        <th>Mapping Efficiency Diff (%)</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for (set1, set2), metrics in pairwise_comparisons.items():
            html += f"""
                <tr>
                    <td><strong>{set1}</strong> vs <strong>{set2}</strong></td>
                    <td>{metrics.coverage_breadth_diff:+.1f}%</td>
                    <td>{metrics.mean_depth_diff:+.1f}x</td>
                    <td>{metrics.gap_count_diff:+d}</td>
                    <td>{metrics.quality_score_diff:+.1f}</td>
                    <td>{metrics.mapping_efficiency_diff:+.1f}%</td>
                </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
    </div>
        """
        
        return html
    
    def _generate_statistical_analysis(self) -> str:
        """Generate statistical analysis section."""
        
        if not self.differential_analyzer:
            return ""
        
        alpha = self.differential_analyzer.significance_level
        method = self.differential_analyzer.correction_method
        quality_tests = self.differential_analyzer.compare_quality_metrics(self.analyzer.oligo_sets)
        
        html = f"""
    <div class="section">
        <h2>Statistical Analysis</h2>
        <p>Paired tests across reference sequences (each reference is one observation).
        P-values are adjusted for multiple comparisons ({method}); alpha = {alpha}.</p>
        
        """
        
        html += "".join(self._render_test(metric.replace('_', ' ').title(), test, alpha)
                        for metric, test in quality_tests.items())
        
        # Add coverage distribution analysis for 2-set comparisons
        if len(self.analyzer.oligo_sets) == 2:
            dist_comparison = self.differential_analyzer.compare_coverage_distributions(
                self.analyzer.oligo_sets[0], self.analyzer.oligo_sets[1]
            )
            
            html += f"""
            <h3>Coverage Distribution Analysis</h3>
            <p>Mean depth per {dist_comparison.window_size} bp window.</p>
            """
            
            for test_name, test in [
                ("Kolmogorov-Smirnov Test", dist_comparison.ks_test),
                ("Mann-Whitney U Test", dist_comparison.mann_whitney_test),
                ("Levene's Test", dist_comparison.levene_test)
            ]:
                html += self._render_test(test_name, test, alpha)
            
            identity_test = self.differential_analyzer.compare_oligo_identity(
                self.analyzer.oligo_sets[0], self.analyzer.oligo_sets[1]
            )
            html += "<h3>Bait Identity</h3>"
            html += self._render_test("Best-hit identity per bait", identity_test, alpha)
        
        html += """
    </div>
        """
        
        return html
    
    @staticmethod
    def _render_test(title: str, test, alpha: float) -> str:
        """Render one StatisticalTest as an HTML block."""
        if not test.applicable:
            return f"""
            <div class="statistical-result non-significant">
                <h4>{title}</h4>
                <p><strong>Test:</strong> {test.test_name}</p>
                <p><strong>Result:</strong> {test.interpretation}</p>
            </div>
            """
        significance_class = "significant" if test.p_reported <= alpha else "non-significant"
        adjusted = f" (adjusted {test.p_adjusted:.4f})" if test.p_adjusted is not None else ""
        return f"""
            <div class="statistical-result {significance_class}">
                <h4>{title}</h4>
                <p><strong>Test:</strong> {test.test_name} (n = {test.n})</p>
                <p><strong>P-value:</strong> {test.p_value:.4f}{adjusted} {test.significance_level}</p>
                <p><strong>Effect Size:</strong> {test.effect_size:.3f}</p>
                <p><strong>Interpretation:</strong> {test.interpretation}</p>
            </div>
            """
    
    def _generate_gap_overlap_analysis(self, gap_analysis) -> str:
        """Generate gap overlap analysis section."""
        
        if not gap_analysis:
            return ""
        
        html = """
    <div class="section">
        <h2>🕳️ Gap Overlap Analysis</h2>
        <p>Analysis of overlapping and unique gaps between oligo sets:</p>
        
        """
        
        html += f"""
        <div class="metric-card">
            <h4>Total Unique Gap Regions</h4>
            <div class="value">{gap_analysis.get('total_unique_regions', 0):,}</div>
        </div>
        """
        
        # Unique gaps by set
        html += """
        <h3>Unique Gaps by Set</h3>
        <table class="comparison-table">
            <thead>
                <tr>
                    <th>Oligo Set</th>
                    <th>Unique Gaps</th>
                    <th>Example Regions</th>
                </tr>
            </thead>
            <tbody>
        """
        
        for set_name, gaps in gap_analysis.get('unique_gaps', {}).items():
            examples = "; ".join([
                f"{gap['chromosome']}:{gap['start']}-{gap['end']}"
                for gap in gaps[:3]  # Show first 3
            ])
            if len(gaps) > 3:
                examples += f" (+{len(gaps) - 3} more)"
            
            html += f"""
                <tr>
                    <td><strong>{set_name}</strong></td>
                    <td>{len(gaps):,}</td>
                    <td>{examples}</td>
                </tr>
            """
        
        html += """
            </tbody>
        </table>
    </div>
        """
        
        return html
    
    def _generate_recommendations(self) -> str:
        """Generate recommendations section."""
        
        best_performer = self.analyzer.identify_best_performer('quality_score')
        ranking = self.analyzer.generate_ranking()
        
        html = """
    <div class="section">
        <h2>💡 Recommendations</h2>
        
        """
        
        # Primary recommendation
        html += f"""
        <div class="recommendation high-priority">
            <h3>🥇 Primary Recommendation</h3>
            <p><strong>Select {best_performer.name}</strong> as your primary oligo set based on its superior overall performance:</p>
            <ul>
                <li>Highest quality score: {best_performer.quality_score.overall_score:.1f}/10 (Grade {best_performer.quality_score.category.value})</li>
                <li>Coverage breadth: {best_performer.coverage_stats.get('coverage_breadth', 0):.1f}%</li>
                <li>Gap count: {best_performer.gap_analysis.get('total_gaps', 0):,}</li>
                <li>Mapping efficiency: {best_performer.coverage_stats.get('mapping_efficiency', 0):.1f}%</li>
            </ul>
        </div>
        """
        
        # Secondary recommendations based on ranking
        if len(ranking) > 1:
            second_best = ranking[1][0]
            second_best_result = next(result for result in self.analyzer.oligo_sets if result.name == second_best)
            
            html += f"""
            <div class="recommendation medium-priority">
                <h3>🥈 Alternative Option</h3>
                <p><strong>{second_best}</strong> is the second-best performer and could be considered if:</p>
                <ul>
                    <li>Cost or availability constraints affect the primary choice</li>
                    <li>Specific design requirements favor this approach</li>
                    <li>Quality score: {second_best_result.quality_score.overall_score:.1f}/10 (Grade {second_best_result.quality_score.category.value})</li>
                </ul>
            </div>
            """
        
        # Improvement recommendations
        html += """
        <div class="recommendation low-priority">
            <h3>🔧 General Improvement Strategies</h3>
            <ul>
        """
        
        # Analysis-based recommendations
        coverage_ranges = [result.coverage_stats.get('coverage_breadth', 0) for result in self.analyzer.oligo_sets]
        if max(coverage_ranges) - min(coverage_ranges) > 20:
            html += "<li>Large variation in coverage breadth suggests some designs could benefit from additional oligos in uncovered regions</li>"
        
        gap_counts = [result.gap_analysis.get('total_gaps', 0) for result in self.analyzer.oligo_sets]
        if max(gap_counts) > 100:
            html += "<li>High gap counts in some sets suggest need for gap-filling optimization using 'baitUtils fill' command</li>"
        
        mapping_effs = [result.coverage_stats.get('mapping_efficiency', 0) for result in self.analyzer.oligo_sets]
        if min(mapping_effs) < 80:
            html += "<li>Low mapping efficiency in some sets suggests reviewing oligo design parameters or reference quality</li>"
        
        html += """
                <li>Consider hybrid approaches combining strengths of top-performing sets</li>
                <li>Use gap overlap analysis to identify consistently problematic regions</li>
                <li>Validate top performers with experimental data before final selection</li>
            </ul>
        </div>
    </div>
        """
        
        return html
    
    def _generate_visualizations_section(self, plots: Dict[str, str]) -> str:
        """Generate visualizations section with embedded plots."""
        
        html = """
    <div class="section">
        <h2>📈 Visualizations</h2>
        
        """
        
        # Embed static plots as base64 images
        plot_titles = {
            'coverage_distributions': 'Coverage Distribution Comparisons',
            'gap_analysis': 'Gap Analysis Comparison',
            'performance_heatmap': 'Performance Heatmap',
            'statistical_comparison': 'Statistical Comparison Results'
        }
        
        for plot_key, plot_file in plots.items():
            if plot_key in plot_titles and Path(plot_file).exists():
                title = plot_titles[plot_key]
                
                # Read and encode image
                try:
                    with open(plot_file, 'rb') as f:
                        img_data = base64.b64encode(f.read()).decode('utf-8')
                    
                    html += f"""
                    <div class="plot-container">
                        <h3>{title}</h3>
                        <img src="data:image/png;base64,{img_data}" alt="{title}">
                    </div>
                    """
                except Exception as e:
                    logging.warning(f"Failed to embed plot {plot_file}: {e}")
        
        # Links to interactive plots
        interactive_plots = {k: v for k, v in plots.items() if 'html' in str(v)}
        if interactive_plots:
            html += """
            <h3>Interactive Visualizations</h3>
            <p>The following interactive plots are available as separate files:</p>
            <ul>
            """
            
            for plot_key, plot_file in interactive_plots.items():
                filename = Path(plot_file).name
                html += f'<li><a href="{filename}" target="_blank">{plot_key.replace("_", " ").title()}</a></li>'
            
            html += "</ul>"
        
        html += "</div>"
        return html
    
    def _generate_methodology_section(self) -> str:
        """Generate methodology section."""
        
        return """
    <div class="section">
        <h2>🔬 Methodology</h2>
        
        <button class="collapsible">Analysis Pipeline</button>
        <div class="content">
            <h3>Comparative Analysis Workflow</h3>
            <ol>
                <li><strong>Mapping:</strong> Oligos mapped to reference using pblat with specified identity thresholds</li>
                <li><strong>Coverage Analysis:</strong> Per-position coverage calculated and statistics computed</li>
                <li><strong>Gap Detection:</strong> Uncovered regions identified and characterized</li>
                <li><strong>Quality Scoring:</strong> Multi-component quality assessment with weighted scores</li>
                <li><strong>Statistical Testing:</strong> Significance testing of differences between sets</li>
                <li><strong>Benchmarking:</strong> Comparison against theoretical optimal performance</li>
                <li><strong>Ranking:</strong> Composite scoring for overall performance ranking</li>
            </ol>
        </div>
        
        <button class="collapsible">Quality Scoring</button>
        <div class="content">
            <h3>Quality Score Components</h3>
            <ul>
                <li><strong>Coverage Breadth (30%):</strong> Percentage of reference sequence covered</li>
                <li><strong>Depth Uniformity (25%):</strong> Evenness of coverage distribution (1 - Gini coefficient)</li>
                <li><strong>Gap Analysis (25%):</strong> Number and size of uncovered regions</li>
                <li><strong>Mapping Efficiency (20%):</strong> Percentage of oligos successfully mapped</li>
            </ul>
            <p><strong>Grading Scale:</strong> A (9.0-10.0), B (8.0-8.9), C (7.0-7.9), D (6.0-6.9), F (<6.0)</p>
        </div>
        
        <button class="collapsible">Statistical Methods</button>
        <div class="content">
            <h3>Statistical Tests Used</h3>
            <ul>
                <li><strong>Kolmogorov-Smirnov Test:</strong> Compare coverage distribution shapes</li>
                <li><strong>Mann-Whitney U Test:</strong> Compare median coverage depths</li>
                <li><strong>Levene's Test:</strong> Test for equal variances in coverage</li>
                <li><strong>t-test/ANOVA:</strong> Compare quality metric means</li>
                <li><strong>Effect Sizes:</strong> Cohen's d, eta-squared for practical significance</li>
            </ul>
            <p><strong>Significance Levels:</strong> *** p≤0.001, ** p≤0.01, * p≤0.05, ns p>0.05</p>
        </div>
    </div>
        """
    
    def _get_html_footer(self) -> str:
        return self._get_html_footer_template().replace("__JS__", load_template("comparative_report.js"))
    
    def _get_html_footer_template(self) -> str:
        return """
    <div class="footer">
        <p>🤖 Generated with <a href="https://claude.ai/code" target="_blank">Claude Code</a></p>
        <p>Co-Authored-By: Claude &lt;noreply@anthropic.com&gt;</p>
        <p>Report generated on """ + datetime.now().strftime("%Y-%m-%d at %H:%M:%S") + """</p>
    </div>

</div>

<script>__JS__</script>

</body>
</html>
        """