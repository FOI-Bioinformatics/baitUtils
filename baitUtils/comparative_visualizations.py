#!/usr/bin/env python3

"""
comparative_visualizations.py

Figures of the compare command. Every set keeps one colour (its position in
--sets maps to a categorical slot) across every figure, so a reader who
learns that a set is blue can follow it through the report.

Static (matplotlib):
- coverage_distributions: breadth, mean depth, uniformity and mapping
  efficiency per set
- gap_analysis_comparison: gap count, share of the reference in gaps, mean
  and largest gap per set
- performance_heatmap: the comparison matrix as a table with a sequential
  colour scale per metric
- statistical_comparison: adjusted p-values, effect sizes, window depth
  distributions and the paired per-reference breadth

Interactive (plotly):
- comparison_dashboard.html: coverage profiles of all sets, reference by
  reference
- metric_comparison.html: every metric relative to the best set
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from baitUtils import design
from baitUtils.design import CHROME, SEQUENTIAL, PLOTLY_CONFIG, style_axes, empty_panel
from baitUtils.comparative_analyzer import OligoSetResult, ComparativeAnalyzer
from baitUtils.differential_analysis import DifferentialAnalyzer

ACCENT = CHROME["accent"]

#: Metrics of the comparison matrix: column, label, higher-is-better, formatter
MATRIX_METRICS = [
    ("Coverage_Breadth_%", "Coverage breadth (%)", True, "{:.1f}"),
    ("Mean_Depth_x", "Mean depth (x)", True, "{:.2f}"),
    ("Mapping_Efficiency_%", "Mapping efficiency (%)", True, "{:.1f}"),
    ("Quality_Score", "Quality score (0-1)", True, "{:.2f}"),
    ("Total_Gaps", "Gaps", False, "{:,.0f}"),
    ("Gini_Coefficient", "Gini coefficient", False, "{:.3f}"),
]


def set_colors(names: List[str]) -> Dict[str, str]:
    """One categorical colour per set, by position."""
    return {name: design.series_color(i) for i, name in enumerate(names)}


class ComparativeVisualizer:
    """Static and interactive figures for the comparison of oligo sets."""

    def __init__(self, output_dir: Path, plot_format: str = 'png', dpi: int = 300, style: str = 'whitegrid'):
        """
        Args:
            output_dir: Directory that receives comparative_plots/ and interactive_comparative_plots/
            plot_format: png, pdf or svg
            dpi: Resolution of raster output
            style: Kept for API compatibility; the design system sets the style
        """
        self.output_dir = Path(output_dir)
        self.plot_format = plot_format
        self.dpi = dpi
        self.plots_dir = self.output_dir / "comparative_plots"
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.interactive_dir = self.output_dir / "interactive_comparative_plots"
        self.interactive_dir.mkdir(parents=True, exist_ok=True)
        design.apply_matplotlib_style(dpi=dpi)

    # ----------------------------------------------------------- helpers
    def _save(self, fig, name: str) -> str:
        path = self.plots_dir / f"{name}.{self.plot_format}"
        plt.savefig(path, dpi=self.dpi, bbox_inches='tight', facecolor=CHROME["surface"])
        plt.close(fig)
        return str(path)

    @staticmethod
    def _bars(ax, names: List[str], values: List[float], colors: Dict[str, str], title: str, xlabel: str,
              fmt=lambda v: f"{v:.1f}", xmax: Optional[float] = None) -> None:
        """Horizontal bars, one per set, in the set's colour, value at the tip."""
        y = np.arange(len(names))[::-1]
        ax.barh(y, values, height=0.58, color=[colors[n] for n in names])
        ax.set_yticks(y)
        ax.set_yticklabels([design.truncate(n, 22) for n in names])
        top = max([v for v in values if np.isfinite(v)] + [0])
        limit = xmax if xmax is not None else top * 1.25 if top > 0 else 1
        ax.set_xlim(0, limit)
        for yi, v in zip(y, values):
            if np.isfinite(v):
                ax.text(v + limit * 0.012, yi, fmt(v), va="center", fontsize=8.5, color=CHROME["ink_secondary"])
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        style_axes(ax, grid="x")
        ax.spines["left"].set_visible(False)

    # ------------------------------------------------------------ static
    def create_coverage_distribution_plots(self, oligo_sets: List[OligoSetResult]) -> str:
        """Breadth, mean depth, uniformity and mapping efficiency per set."""
        names = [s.name for s in oligo_sets]
        colors = set_colors(names)
        fig, axes = plt.subplots(2, 2, figsize=(12, 3.0 + 0.5 * len(names) + 3.0))
        stats = [s.coverage_stats for s in oligo_sets]
        self._bars(axes[0, 0], names, [c.get('coverage_breadth', 0) for c in stats], colors,
                   "Coverage breadth", "Bases covered (%)", lambda v: f"{v:.1f}%", xmax=115)
        self._bars(axes[0, 1], names, [c.get('mean_depth', 0) for c in stats], colors,
                   "Mean depth", "Aligned baits per base", lambda v: f"{v:.2f}x")
        self._bars(axes[1, 0], names, [1 - c.get('coverage_gini', 1) for c in stats], colors,
                   "Uniformity", "1 - Gini (1 = perfectly even)", lambda v: f"{v:.2f}", xmax=1.15)
        self._bars(axes[1, 1], names, [c.get('mapping_efficiency', 0) for c in stats], colors,
                   "Mapping efficiency", "Baits with a hit (%)", lambda v: f"{v:.1f}%", xmax=115)
        fig.suptitle("Coverage metrics by set", x=0.01, ha="left")
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        return self._save(fig, "coverage_distributions")

    def create_gap_analysis_comparison(self, oligo_sets: List[OligoSetResult]) -> str:
        """Gap count, share of the reference in gaps, mean and largest gap per set."""
        names = [s.name for s in oligo_sets]
        colors = set_colors(names)
        gaps = [s.gap_analysis for s in oligo_sets]
        fig, axes = plt.subplots(2, 2, figsize=(12, 3.0 + 0.5 * len(names) + 3.0))
        self._bars(axes[0, 0], names, [g.get('total_gaps', 0) for g in gaps], colors,
                   "Gaps", "Gaps at or above the minimum gap size", lambda v: f"{v:,.0f}")
        self._bars(axes[0, 1], names, [g.get('gap_percentage', 0) for g in gaps], colors,
                   "Share of the reference in gaps", "Reference in gaps (%)", lambda v: f"{v:.1f}%")
        self._bars(axes[1, 0], names, [g.get('mean_gap_size', 0) for g in gaps], colors,
                   "Mean gap length", "bp", design.fmt_bp)
        self._bars(axes[1, 1], names, [g.get('max_gap_size', 0) for g in gaps], colors,
                   "Largest gap", "bp", design.fmt_bp)
        fig.suptitle("Gap metrics by set", x=0.01, ha="left")
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        return self._save(fig, "gap_analysis_comparison")

    def create_performance_heatmap(self, analyzer: ComparativeAnalyzer) -> str:
        """The comparison matrix as a table: raw values, shaded per metric from worst to best."""
        matrix = analyzer.generate_comparison_matrix()
        names = matrix['Name'].tolist()
        columns = [m for m in MATRIX_METRICS if m[0] in matrix.columns]
        shade = pd.DataFrame(index=names)
        labels = pd.DataFrame(index=names)
        for column, label, higher_better, fmt in columns:
            values = matrix[column].astype(float).to_numpy()
            lo, hi = np.nanmin(values), np.nanmax(values)
            scaled = (values - lo) / (hi - lo) if hi > lo else np.full(values.shape, 0.5)
            shade[label] = scaled if higher_better else 1 - scaled
            labels[label] = [fmt.format(v) for v in values]
        fig, ax = plt.subplots(figsize=(1.6 * len(columns) + 2.5, 0.6 * len(names) + 1.8))
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("baitutils_blue", ["#ffffff"] + SEQUENTIAL[:9])
        sns.heatmap(shade, annot=labels.to_numpy(), fmt="", cmap=cmap, vmin=0, vmax=1, cbar=False,
                    linewidths=2, linecolor=CHROME["surface"], ax=ax,
                    annot_kws={"size": 9, "color": CHROME["ink"]})
        # seaborn writes the labels row by row; dark cells get light text
        for text, value in zip(ax.texts, shade.to_numpy().ravel()):
            text.set_color(design.label_in_fill(matplotlib.colors.to_hex(cmap(float(value)))))
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.xaxis.tick_top()
        ax.tick_params(length=0, labelsize=9)
        plt.setp(ax.get_xticklabels(), rotation=0, ha="center")
        ax.set_title("Comparison matrix (darker is better within each column)", pad=28)
        fig.tight_layout()
        return self._save(fig, "performance_heatmap")

    def create_statistical_comparison_plot(self, differential_analyzer: DifferentialAnalyzer,
                                           oligo_sets: List[OligoSetResult]) -> str:
        """Adjusted p-values, effect sizes, window depth distributions and paired per-reference breadth."""
        if len(oligo_sets) < 2:
            raise ValueError("Need at least 2 oligo sets for statistical comparison")
        names = [s.name for s in oligo_sets]
        colors = set_colors(names)
        alpha = differential_analyzer.significance_level
        tests = differential_analyzer.compare_quality_metrics(oligo_sets)
        rows: List[tuple] = [(m.replace('_', ' '), t) for m, t in tests.items()]
        if len(oligo_sets) == 2:
            dist = differential_analyzer.compare_coverage_distributions(oligo_sets[0], oligo_sets[1])
            rows += [("window depth: KS", dist.ks_test), ("window depth: Mann-Whitney", dist.mann_whitney_test),
                     ("window depth: Levene", dist.levene_test)]
            identity = differential_analyzer.compare_oligo_identity(oligo_sets[0], oligo_sets[1])
            rows.append(("bait identity", identity))
            gap_test = differential_analyzer.analyze_gap_patterns(oligo_sets[0], oligo_sets[1])
            rows.append(("gap lengths", gap_test['gap_size_distribution']))
        live = [(label, t) for label, t in rows if t.applicable]

        fig, axes = plt.subplots(2, 2, figsize=(12, 9))
        ax = axes[0, 0]
        if live:
            y = np.arange(len(live))[::-1]
            p = np.array([max(t.p_reported, 1e-300) for _, t in live])
            ax.hlines(y, 0, -np.log10(p), color=CHROME["grid"], linewidth=2)
            ax.plot(-np.log10(p), y, "o", color=ACCENT, markersize=7, markeredgecolor=CHROME["surface"], markeredgewidth=1.5)
            ax.axvline(-np.log10(alpha), color=CHROME["ink_muted"], linewidth=0.8, linestyle=(0, (3, 3)))
            ax.text(-np.log10(alpha), len(live) - 0.4, f" alpha = {alpha:g}", fontsize=8, color=CHROME["ink_muted"])
            ax.set_yticks(y)
            ax.set_yticklabels([f"{label} (n={t.n})" for label, t in live], fontsize=8.5)
            ax.set_xlabel("-log10(adjusted p)")
            ax.set_xlim(0, max(float(np.max(-np.log10(p))), -np.log10(alpha)) * 1.3 + 0.1)
            style_axes(ax, grid="x")
            ax.spines["left"].set_visible(False)
        else:
            empty_panel(ax, "No applicable test")
        ax.set_title("Tests (adjusted p-values)")

        ax = axes[0, 1]
        effects = [(label, t) for label, t in live if np.isfinite(t.effect_size)]
        if effects:
            y = np.arange(len(effects))[::-1]
            values = [t.effect_size for _, t in effects]
            ax.barh(y, values, height=0.55, color=[ACCENT if v >= 0 else CHROME["gap_fill"] for v in values])
            ax.axvline(0, color=CHROME["axis"], linewidth=0.8)
            ax.set_yticks(y)
            ax.set_yticklabels([label for label, _ in effects], fontsize=8.5)
            ax.set_xlabel("Effect size (rank-biserial r, paired d, or Kendall's W)")
            style_axes(ax, grid="x")
            ax.spines["left"].set_visible(False)
        else:
            empty_panel(ax, "No effect size to show")
        ax.set_title("Effect sizes (positive: first set larger)")

        ax = axes[1, 0]
        window = differential_analyzer.window_size
        data = []
        for s in oligo_sets:
            try:
                data.append(differential_analyzer._window_means(s))
            except ValueError:
                data.append(np.array([]))
        if any(d.size for d in data):
            parts = ax.boxplot([d if d.size else [np.nan] for d in data], vert=True, widths=0.5, patch_artist=True,
                               showfliers=len(oligo_sets) <= 4, medianprops=dict(color=CHROME["ink"], linewidth=1.2),
                               whiskerprops=dict(color=CHROME["ink_muted"]), capprops=dict(color=CHROME["ink_muted"]),
                               flierprops=dict(marker="o", markersize=3, markerfacecolor=CHROME["ink_muted"],
                                               markeredgecolor="none", alpha=0.6))
            for box, name in zip(parts['boxes'], names):
                box.set(facecolor=design.rgba(colors[name], 0.35) if False else colors[name], alpha=0.45, linewidth=0)
            ax.set_xticks(range(1, len(names) + 1))
            ax.set_xticklabels([design.truncate(n, 16) for n in names])
            ax.set_ylabel(f"Mean depth per {window:,} bp window")
            style_axes(ax, grid="y")
        else:
            empty_panel(ax, "No per-base coverage for window depths")
        ax.set_title("Window depth distribution")

        ax = axes[1, 1]
        refs = sorted(set.intersection(*(set(s.coverage_stats.get('per_reference', {})) for s in oligo_sets))) \
            if oligo_sets else []
        if refs:
            x = np.arange(len(names))
            for ref in refs[:60]:
                y = [s.coverage_stats['per_reference'][ref].get('coverage_breadth', np.nan) for s in oligo_sets]
                ax.plot(x, y, color=CHROME["emphasis_gray"], linewidth=1, alpha=0.7, zorder=1)
            for i, name in enumerate(names):
                values = [s.coverage_stats['per_reference'][r].get('coverage_breadth', np.nan)
                          for s in oligo_sets if s.name == name for r in refs[:60]]
                ax.plot([i] * len(values), values, "o", color=colors[name], markersize=5,
                        markeredgecolor=CHROME["surface"], markeredgewidth=1, zorder=2)
            ax.set_xticks(x)
            ax.set_xticklabels([design.truncate(n, 16) for n in names])
            ax.set_ylabel("Breadth per reference (%)")
            ax.set_ylim(-3, 105)
            style_axes(ax, grid="y")
            note = f"{len(refs)} references" + (", first 60 drawn" if len(refs) > 60 else "")
            ax.set_title(f"Paired breadth per reference ({note})")
        else:
            empty_panel(ax, "No shared references")
            ax.set_title("Paired breadth per reference")
        fig.suptitle("Statistical comparison", x=0.01, ha="left")
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        return self._save(fig, "statistical_comparison")

    # ------------------------------------------------------- interactive
    def coverage_profiles_figure(self, oligo_sets: List[OligoSetResult], max_references: int = 6):
        """Coverage profiles of every set, one panel per reference; None without arrays."""
        with_arrays = [s for s in oligo_sets if s.coverage_arrays]
        if not with_arrays:
            return None
        lengths: Dict[str, int] = {}
        for s in with_arrays:
            for ref, array in s.coverage_arrays.items():
                lengths[ref] = max(lengths.get(ref, 0), len(array))
        refs = design.select_references({r: {'length': n} for r, n in lengths.items()}, max_references)
        names = [s.name for s in oligo_sets]
        colors = set_colors(names)
        fig = make_subplots(rows=len(refs), cols=1, subplot_titles=[design.truncate(r, 40) for r in refs],
                            vertical_spacing=min(0.12, 0.6 / max(len(refs), 1)))
        for i, ref in enumerate(refs, 1):
            for s in oligo_sets:
                array = (s.coverage_arrays or {}).get(ref)
                if array is None:
                    continue
                x, y = design.bin_coverage(array, 1000)
                fig.add_trace(go.Scatter(x=x, y=y, mode="lines", line=dict(color=colors[s.name], width=1.5),
                                         name=s.name, legendgroup=s.name, showlegend=(i == 1),
                                         hovertemplate=s.name + "<br>position %{x:,.0f}<br>depth %{y:.2f}<extra></extra>"),
                              row=i, col=1)
            fig.update_xaxes(range=[0, lengths[ref]], row=i, col=1)
            fig.update_yaxes(rangemode="tozero", title_text="Depth" if i == 1 else None, row=i, col=1)
        fig.update_xaxes(title_text="Position along reference (bp)", row=len(refs), col=1)
        for annotation in fig.layout.annotations:
            annotation.update(font=dict(size=12, color=CHROME["ink"]), x=0, xanchor="left")
        title = "Coverage profiles by set"
        if len(lengths) > len(refs):
            title += f" ({len(refs)} longest of {len(lengths)} references)"
        design.plotly_layout(fig, title, 170 * len(refs) + 140, legend_above_subplots=True)
        return fig

    def create_comparison_dashboard(self, analyzer: ComparativeAnalyzer) -> str:
        """Coverage profiles of all sets, reference by reference."""
        if not analyzer.oligo_sets:
            raise ValueError("No oligo sets found in analyzer")
        fig = self.coverage_profiles_figure(analyzer.oligo_sets)
        if fig is None:
            fig = go.Figure()
            design.plotly_layout(fig, "Coverage profiles by set", 360)
            fig.add_annotation(text="No per-base coverage available", showarrow=False, x=0.5, y=0.5,
                               xref="paper", yref="paper", font=dict(color=CHROME["ink_muted"]))
        path = self.interactive_dir / "comparison_dashboard.html"
        fig.write_html(str(path), config=PLOTLY_CONFIG, include_plotlyjs="cdn")
        return str(path)

    def metric_comparison_figure(self, analyzer: ComparativeAnalyzer):
        """Every metric relative to the best set, sets as colours."""
        matrix = analyzer.generate_comparison_matrix()
        names = matrix['Name'].tolist()
        colors = set_colors(names)
        columns = [m for m in MATRIX_METRICS if m[0] in matrix.columns]
        fig = go.Figure()
        for i, name in enumerate(names):
            relative, text = [], []
            for column, label, higher_better, fmt in columns:
                values = matrix[column].astype(float).to_numpy()
                v = values[i]
                best = np.nanmax(values) if higher_better else np.nanmin(values)
                if higher_better:
                    rel = v / best if best > 0 else 1.0
                else:
                    rel = best / v if v > 0 else 1.0
                relative.append(min(rel, 1.0) if np.isfinite(rel) else 0.0)
                text.append(f"{name}: {fmt.format(v)} (best {fmt.format(best)})")
            fig.add_trace(go.Scatter(
                x=relative, y=[m[1] for m in columns], mode="markers", name=name,
                marker=dict(size=12, color=colors[name], line=dict(color=CHROME["surface"], width=1.5)),
                text=text, hovertemplate="%{text}<extra></extra>"))
        design.plotly_layout(fig, "Metrics relative to the best set", 120 + 48 * len(columns),
                             xaxis=dict(title="Fraction of the best value (1 = best set)", range=[0, 1.08]),
                             yaxis=dict(autorange="reversed"))
        return fig

    def create_interactive_ranking_plot(self, analyzer: ComparativeAnalyzer) -> str:
        """Write the metric comparison figure."""
        if not analyzer.generate_ranking():
            raise ValueError("No ranking data available")
        fig = self.metric_comparison_figure(analyzer)
        path = self.interactive_dir / "metric_comparison.html"
        fig.write_html(str(path), config=PLOTLY_CONFIG, include_plotlyjs="cdn")
        return str(path)

    # ------------------------------------------------------------ driver
    def generate_all_comparative_plots(self, analyzer: ComparativeAnalyzer,
                                       differential_analyzer: Optional[DifferentialAnalyzer] = None) -> Dict[str, str]:
        """Generate all comparative figures; returns {name: path}."""
        plots: Dict[str, str] = {}
        steps = [
            ('dashboard', lambda: self.create_comparison_dashboard(analyzer)),
            ('coverage_distributions', lambda: self.create_coverage_distribution_plots(analyzer.oligo_sets)),
            ('gap_analysis', lambda: self.create_gap_analysis_comparison(analyzer.oligo_sets)),
            ('performance_heatmap', lambda: self.create_performance_heatmap(analyzer)),
            ('metric_comparison', lambda: self.create_interactive_ranking_plot(analyzer)),
        ]
        if differential_analyzer and len(analyzer.oligo_sets) >= 2:
            steps.append(('statistical_comparison',
                          lambda: self.create_statistical_comparison_plot(differential_analyzer, analyzer.oligo_sets)))
        for key, make in steps:
            try:
                plots[key] = make()
            except Exception as e:
                logging.error(f"Failed to create the {key} figure: {e}")
        logging.info(f"Generated {len(plots)} comparative figures")
        return plots
