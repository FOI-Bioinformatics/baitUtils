#!/usr/bin/env python3

"""
interactive_plots.py

Interactive figures of the evaluate command (plotly), styled with the
design system of baitUtils.design:

- coverage_dashboard.html: zoomable depth profile per reference with gaps
  shaded, and the cumulative depth curve
- detailed_gap_analysis.html: gaps by position and length, the largest gaps,
  and gap length against the gap's GC content
- quality_assessment.html: component scores with their weights, observed
  metrics against the design targets, the overall score against the
  category thresholds
- reference_analysis.html: per-reference features against breadth, breadth
  by GC decile, and the challenging windows
"""

import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from baitUtils import design
from baitUtils.design import CHROME, PLOTLY_CONFIG
from baitUtils.quality_scorer import QualityScorer

ACCENT = CHROME["accent"]
GAP = CHROME["gap_fill"]
MAX_PROFILE_REFERENCES = 8


SHORT_COMPONENT = {
    'coverage_breadth': 'Breadth', 'coverage_depth': 'Depth', 'mapping_efficiency': 'Mapping efficiency',
    'gap_characteristics': 'Gaps', 'reference_difficulty': 'Reference difficulty',
}


def _layout(fig, title: str, height: int, **kwargs) -> None:
    design.plotly_layout(fig, title, height, **kwargs)


def _note(fig, text: str, row: int, col: int) -> None:
    """A panel with nothing to plot says so."""
    design.subplot_note(fig, text, row, col)


class InteractivePlotter:
    """Generate interactive plots for coverage analysis."""

    def __init__(
        self,
        coverage_stats: Dict[str, Any],
        gap_analysis: Dict[str, Any],
        quality_scores: Dict[str, Any],
        reference_analysis: Dict[str, Any],
        output_dir: Path,
        coverage_arrays: Optional[Dict[str, np.ndarray]] = None,
    ):
        """
        Args:
            coverage_stats: Output of CoverageAnalyzer.analyze
            gap_analysis: Output of GapAnalyzer.analyze
            quality_scores: QualityScore.to_dict()
            reference_analysis: Output of ReferenceAnalyzer.analyze
            output_dir: Directory that receives the interactive_plots/ folder
            coverage_arrays: Per-reference depth arrays for the profile
        """
        self.coverage_stats = coverage_stats or {}
        self.gap_analysis = gap_analysis or {}
        self.quality_scores = quality_scores or {}
        self.reference_analysis = reference_analysis or {}
        self.coverage_arrays = coverage_arrays or {}
        self.output_dir = Path(output_dir)
        self.plots_dir = self.output_dir / "interactive_plots"
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.default_config = dict(PLOTLY_CONFIG)

    # ----------------------------------------------------------- helpers
    @property
    def _min_coverage(self) -> float:
        return float(self.coverage_stats.get('parameters', {}).get('min_coverage', 1.0))

    @property
    def _target_coverage(self) -> float:
        return float(self.coverage_stats.get('parameters', {}).get('target_coverage', 10.0))

    def _write(self, fig, name: str) -> str:
        path = self.plots_dir / name
        fig.write_html(str(path), config=self.default_config, include_plotlyjs="cdn")
        logging.info(f"Interactive figure saved: {path}")
        return str(path)

    def _gaps_by_reference(self) -> Dict[str, List[Dict[str, Any]]]:
        gaps: Dict[str, List[Dict[str, Any]]] = {}
        for gap in self.gap_analysis.get('gaps', []) or []:
            gaps.setdefault(gap.get('chromosome', ''), []).append(gap)
        return gaps

    # ----------------------------------------------------------- figures
    def coverage_profile_figure(self, max_references: int = MAX_PROFILE_REFERENCES):
        """Depth profile per reference as a plotly figure, or None without arrays."""
        if not self.coverage_arrays:
            return None
        refs = design.select_references({r: {'length': len(a)} for r, a in self.coverage_arrays.items()}, max_references)
        gaps = self._gaps_by_reference()
        per_ref = self.coverage_stats.get('per_reference', {}) or {}
        n = len(refs)
        titles = []
        for ref in refs:
            s = per_ref.get(ref, {})
            titles.append(f"{design.truncate(ref, 40)} - breadth {s.get('coverage_breadth', 0):.1f}%, "
                          f"mean depth {design.fmt_depth(s.get('mean_depth', 0))}, {len(gaps.get(ref, []))} gap(s)")
        fig = make_subplots(rows=n, cols=1, subplot_titles=titles, vertical_spacing=min(0.12, 0.6 / max(n, 1)))
        for i, ref in enumerate(refs, 1):
            array = self.coverage_arrays[ref]
            x, y = design.bin_coverage(array, 1200)
            fig.add_trace(go.Scatter(
                x=x, y=y, mode="lines", line=dict(color=ACCENT, width=1.5), fill="tozeroy",
                fillcolor=design.rgba(ACCENT, 0.12), name=ref, showlegend=False,
                hovertemplate="position %{x:,.0f}<br>depth %{y:.2f}<extra></extra>"), row=i, col=1)
            for gap in gaps.get(ref, []):
                fig.add_vrect(x0=gap['start'], x1=gap['end'], fillcolor=design.rgba(GAP, 0.18), line_width=0,
                              row=i, col=1)
            fig.add_hline(y=self._min_coverage, line=dict(color=CHROME["ink_muted"], width=1, dash="dot"), row=i, col=1)
            fig.update_xaxes(range=[0, len(array)], row=i, col=1)
            fig.update_yaxes(rangemode="tozero", title_text="Depth" if i == 1 else None, row=i, col=1)
        fig.update_xaxes(title_text="Position along reference (bp)", row=n, col=1)
        for annotation in fig.layout.annotations:
            annotation.update(font=dict(size=12, color=CHROME["ink"]), x=0, xanchor="left")
        title = "Coverage profile"
        if len(self.coverage_arrays) > n:
            title += f" ({n} longest of {len(self.coverage_arrays)} references)"
        _layout(fig, title, 170 * n + 120, margin=dict(l=60, r=24, t=70, b=56))
        return fig

    def create_coverage_dashboard(self) -> str:
        """Depth profile per reference and the cumulative depth curve."""
        logging.info("Creating interactive coverage figure...")
        profile = self.coverage_profile_figure()
        if profile is not None:
            fig = profile
            fig.update_layout(title=dict(text=fig.layout.title.text))
        else:
            fig = make_subplots(rows=1, cols=1)
            _layout(fig, "Coverage profile", 400)
            _note(fig, "No per-base coverage available for a profile", 1, 1)
        return self._write(fig, "coverage_dashboard.html")

    def cumulative_depth_figure(self):
        """Fraction of bases with depth at least d."""
        fig = go.Figure()
        if self.coverage_arrays:
            values = np.concatenate([np.asarray(a) for a in self.coverage_arrays.values()])
            upper = int(max(np.percentile(values, 99.5), self._target_coverage, 1))
            depths = np.arange(0, upper + 2)
            fraction = [100.0 * float(np.mean(values >= d)) for d in depths]
            fig.add_trace(go.Scatter(x=depths, y=fraction, mode="lines", line=dict(color=ACCENT, width=2, shape="hv"),
                                     fill="tozeroy", fillcolor=design.rgba(ACCENT, 0.1), name="Bases",
                                     hovertemplate="depth ≥ %{x}: %{y:.1f}% of bases<extra></extra>"))
        else:
            dist = self.coverage_stats.get('depth_distribution', {}) or {}
            if dist:
                thresholds = sorted(dist)
                fig.add_trace(go.Scatter(x=thresholds, y=[dist[t] for t in thresholds], mode="lines+markers",
                                         line=dict(color=ACCENT, width=2), marker=dict(size=8), name="Bases",
                                         hovertemplate="depth ≥ %{x}: %{y:.1f}% of bases<extra></extra>"))
        for value, label, color in ((self._min_coverage, "minimum", CHROME["ink_muted"]),
                                    (self._target_coverage, "target", CHROME["ink_secondary"])):
            fig.add_vline(x=value, line=dict(color=color, width=1, dash="dot"),
                          annotation_text=f"{label} {value:g}x", annotation_position="top right",
                          annotation_font=dict(size=11, color=color))
        fig.update_xaxes(title_text="Depth d")
        fig.update_yaxes(title_text="Bases with depth ≥ d (%)", range=[0, 105])
        _layout(fig, "Cumulative depth", 360, showlegend=False)
        return fig

    def create_detailed_gap_analysis(self) -> str:
        """Write the gap figure."""
        logging.info("Creating interactive gap figure...")
        return self._write(self.gap_figure(), "detailed_gap_analysis.html")

    def gap_figure(self):
        """Gaps by position and length, the largest gaps, and gap length against GC content."""
        gaps = self.gap_analysis.get('gaps', []) or []
        features = self.gap_analysis.get('gap_features', []) or []
        fig = make_subplots(rows=2, cols=2, subplot_titles=[
            "Gaps along the references", "Largest gaps",
            "Gap length against GC content", "Gap length distribution"],
            vertical_spacing=0.16, horizontal_spacing=0.1)

        if gaps:
            refs = design.ordered_unique(g.get('chromosome', '') for g in gaps)
            shown = refs[:8]
            for i, ref in enumerate(shown):
                sel = [g for g in gaps if g.get('chromosome') == ref]
                fig.add_trace(go.Scatter(
                    x=[g['start'] for g in sel], y=[g['length'] for g in sel], mode="markers",
                    marker=dict(size=9, color=design.series_color(i), line=dict(color=CHROME["surface"], width=1.5)),
                    name=design.truncate(ref, 24), legendgroup=ref,
                    text=[f"{ref}:{g['start']:,}-{g['end']:,}" for g in sel],
                    hovertemplate="%{text}<br>%{y:,} bp<extra></extra>"), row=1, col=1)
            if len(refs) > len(shown):
                rest = [g for g in gaps if g.get('chromosome') not in shown]
                fig.add_trace(go.Scatter(
                    x=[g['start'] for g in rest], y=[g['length'] for g in rest], mode="markers",
                    marker=dict(size=8, color=CHROME["emphasis_gray"]), name=f"{len(refs) - len(shown)} other references",
                    text=[f"{g['chromosome']}:{g['start']:,}-{g['end']:,}" for g in rest],
                    hovertemplate="%{text}<br>%{y:,} bp<extra></extra>"), row=1, col=1)
            lengths = [g['length'] for g in gaps]
            log_range = [math.log10(min(lengths)) - 0.6, math.log10(max(lengths)) + 0.6]
            per_ref = self.coverage_stats.get('per_reference', {}) or {}
            longest = max((r.get('length', 0) for r in per_ref.values()), default=0)
            fig.update_xaxes(title_text="Gap start (bp)", range=[0, longest * 1.02] if longest else None, row=1, col=1)
            fig.update_yaxes(title_text="Gap length (bp)", type="log", dtick=1, range=log_range, row=1, col=1)

            largest = sorted(gaps, key=lambda g: g['length'], reverse=True)[:15][::-1]
            fig.add_trace(go.Bar(
                y=[f"{design.truncate(g['chromosome'], 18)}:{g['start']:,}-{g['end']:,}" for g in largest],
                x=[g['length'] for g in largest], orientation="h", marker_color=GAP, opacity=0.85,
                showlegend=False, hovertemplate="%{y}<br>%{x:,} bp<extra></extra>"), row=1, col=2)
            fig.update_xaxes(title_text="Gap length (bp)", row=1, col=2)
            fig.update_yaxes(tickfont=dict(size=10), row=1, col=2)

            if features:
                fig.add_trace(go.Scatter(
                    x=[f.get('gc_content', 0) for f in features], y=[f.get('length', 0) for f in features],
                    mode="markers", marker=dict(size=9, color=ACCENT, line=dict(color=CHROME["surface"], width=1.5)),
                    showlegend=False,
                    text=[f"{f.get('chromosome')}:{f.get('start'):,}-{f.get('end'):,}" for f in features],
                    hovertemplate="%{text}<br>GC %{x:.1f}%<br>%{y:,} bp<extra></extra>"), row=2, col=1)
                fig.update_xaxes(title_text="GC content of the gap (%)", range=[0, 100], row=2, col=1)
                fig.update_yaxes(title_text="Gap length (bp)", type="log", dtick=1, range=log_range, row=2, col=1)
            else:
                _note(fig, "No gap features", 2, 1)

            bins = design.histogram_bins(lengths, min(30, max(5, int(np.sqrt(len(lengths))) + 3)))
            fig.add_trace(go.Histogram(x=lengths, xbins=bins, autobinx=False,
                                       marker_color=ACCENT, showlegend=False,
                                       hovertemplate="%{x} bp: %{y} gaps<extra></extra>"), row=2, col=2)
            fig.update_xaxes(title_text="Gap length (bp)", row=2, col=2)
            fig.update_yaxes(title_text="Gaps", row=2, col=2)
        else:
            for row, col in ((1, 1), (1, 2), (2, 1), (2, 2)):
                _note(fig, "No gaps at or above the minimum gap size", row, col)
        _layout(fig, f"Gap analysis ({len(gaps):,} gaps)", 820, legend_above_subplots=True,
                showlegend=len(design.ordered_unique(g.get('chromosome', '') for g in gaps)) > 1)
        return fig

    def create_quality_assessment_plot(self) -> str:
        """Write the quality figure."""
        logging.info("Creating interactive quality figure...")
        return self._write(self.quality_figure(), "quality_assessment.html")

    def quality_figure(self):
        """Component scores, observed metrics against targets, overall score against thresholds."""
        fig = make_subplots(rows=1, cols=3, column_widths=[0.38, 0.36, 0.26], horizontal_spacing=0.13,
                            subplot_titles=["Component scores", "Observed against target", "Overall score"])
        scores = self.quality_scores.get('component_scores', {}) or {}
        weights = QualityScorer.DEFAULT_WEIGHTS
        if scores:
            names = list(scores)[::-1]
            labels = [f"{SHORT_COMPONENT.get(n, n.replace('_', ' ').capitalize())} ({weights.get(n, 0):.0%})"
                      for n in names]
            fig.add_trace(go.Bar(y=labels, x=[scores[n] for n in names], orientation="h", marker_color=ACCENT,
                                 name="Component score", showlegend=False,
                                 hovertemplate="%{y}: %{x:.2f}<extra></extra>"), row=1, col=1)
            fig.update_xaxes(range=[0, 1.05], title_text="Score (0-1), weight in brackets", row=1, col=1)
            fig.update_yaxes(automargin=True, row=1, col=1)
        else:
            _note(fig, "No quality scores", 1, 1)

        targets = QualityScorer.DEFAULT_BENCHMARKS
        observed = {
            "Breadth": (self.coverage_stats.get('coverage_breadth'), targets['coverage_breadth']),
            "Mean depth": (self.coverage_stats.get('mean_depth'), targets['mean_depth']),
            "Uniformity": (None if self.coverage_stats.get('coverage_gini') is None
                           else 1 - self.coverage_stats['coverage_gini'], targets['depth_uniformity']),
            "Mapping": (self.coverage_stats.get('mapping_efficiency'), targets['mapping_efficiency']),
        }
        rows = [(k, v, t) for k, (v, t) in observed.items() if v is not None]
        if rows:
            ratio = [min(v / t, 1.5) if t else 0 for _, v, t in rows]
            fig.add_trace(go.Scatter(x=[1.0] * len(rows), y=[r[0] for r in rows], mode="markers",
                                     marker=dict(symbol="line-ns", size=18, color=CHROME["ink_muted"],
                                                 line=dict(width=2, color=CHROME["ink_muted"])),
                                     name="Target", hoverinfo="skip"), row=1, col=2)
            fig.add_trace(go.Scatter(x=ratio, y=[r[0] for r in rows], mode="markers",
                                     marker=dict(size=11, color=ACCENT, line=dict(color=CHROME["surface"], width=1.5)),
                                     name="Observed",
                                     text=[f"{v:.2f} of target {t:g}" for _, v, t in rows],
                                     hovertemplate="%{y}: %{text}<extra></extra>"), row=1, col=2)
            fig.update_xaxes(range=[0, 1.6], title_text="Observed / target (1 = target reached)", row=1, col=2)
            fig.update_yaxes(tickfont=dict(size=11), row=1, col=2)
        else:
            _note(fig, "No observed metrics", 1, 2)

        overall = float(self.quality_scores.get('overall_score', 0.0))
        category = str(self.quality_scores.get('category', ''))
        fill = design.category_color(category)
        fig.add_trace(go.Bar(x=[0], y=[overall], width=0.5, marker_color=fill,
                             text=[f"{overall:.2f}"], textposition="inside", insidetextanchor="end",
                             textfont=dict(color=design.label_in_fill(fill), size=13),
                             showlegend=False, hovertemplate=f"{overall:.3f} {category}<extra></extra>"), row=1, col=3)
        for name, threshold in QualityScorer.DEFAULT_THRESHOLDS.items():
            if threshold > 0:
                fig.add_shape(type="line", x0=-0.45, x1=0.95, y0=threshold, y1=threshold,
                              line=dict(dash="dot", color=CHROME["ink_muted"]), row=1, col=3)
                fig.add_annotation(x=0.93, y=threshold, text=name.capitalize(), showarrow=False, xanchor="right",
                                   yanchor="bottom", font=dict(color=CHROME["ink_muted"], size=10), row=1, col=3)
        fig.update_xaxes(range=[-0.45, 0.95], tickvals=[0], ticktext=[category or "Overall"], row=1, col=3)
        fig.update_yaxes(range=[0, 1.05], title_text="Score (0-1)", row=1, col=3)
        _layout(fig, "Quality assessment", 460, legend_above_subplots=True)
        return fig

    def create_reference_analysis_plot(self) -> str:
        """Write the reference figure."""
        logging.info("Creating interactive reference figure...")
        return self._write(self.reference_figure(), "reference_analysis.html")

    def reference_figure(self):
        """Per-reference features against breadth, breadth by GC decile, and challenging windows."""
        features = self.reference_analysis.get('sequence_features', {}) or {}
        per_ref = self.coverage_stats.get('per_reference', {}) or {}
        correlations = self.reference_analysis.get('coverage_correlations', {}) or {}
        windows = self.reference_analysis.get('challenging_windows', {}) or {}

        fig = make_subplots(rows=2, cols=2, subplot_titles=[
            "Breadth against GC content", "Breadth against repeat content",
            "Breadth by GC decile (windows)", "Challenging windows"],
            vertical_spacing=0.16, horizontal_spacing=0.1)

        refs = [r for r in features if r in per_ref] or list(features)
        if refs:
            breadth = [per_ref.get(r, {}).get('coverage_breadth', np.nan) for r in refs]
            lengths = np.array([features[r].get('length', 0) for r in refs], dtype=float)
            size = 8 + 10 * (lengths / lengths.max() if lengths.max() > 0 else 0)
            for col, key, label in ((1, 'gc_content', 'GC content (%)'), (2, 'repeat_content', 'Repeat content (%)')):
                fig.add_trace(go.Scatter(
                    x=[features[r].get(key, np.nan) for r in refs], y=breadth, mode="markers",
                    marker=dict(size=size, color=ACCENT, line=dict(color=CHROME["surface"], width=1.5)),
                    text=refs, showlegend=False,
                    hovertemplate="%{text}<br>" + label + " %{x:.1f}<br>breadth %{y:.1f}%<extra></extra>"), row=1, col=col)
                fig.update_xaxes(title_text=label, row=1, col=col)
                fig.update_yaxes(title_text="Coverage breadth (%)", range=[-3, 105], row=1, col=col)
        else:
            _note(fig, "No reference features", 1, 1)
            _note(fig, "No reference features", 1, 2)

        deciles = correlations.get('breadth_by_gc', []) or []
        if deciles:
            fig.add_trace(go.Bar(
                x=[f"{d['gc_min']}-{d['gc_max']}%" for d in deciles], y=[d['mean_breadth'] for d in deciles],
                marker_color=ACCENT, showlegend=False, text=[f"n={d['windows']:,}" for d in deciles],
                textposition="outside", textfont=dict(size=10, color=CHROME["ink_muted"]),
                hovertemplate="GC %{x}<br>mean breadth %{y:.1f}%<br>%{text} windows<extra></extra>"), row=2, col=1)
            fig.update_xaxes(title_text="GC content of the window", row=2, col=1)
            fig.update_yaxes(title_text="Mean breadth (%)", range=[0, 112], row=2, col=1)
        else:
            _note(fig, "Correlations not applicable" if not correlations.get('applicable') else "No windows", 2, 1)

        rows = [(ref, w) for ref, items in windows.items() for w in items]
        if rows:
            fig.add_trace(go.Scatter(
                x=[w['start'] for _, w in rows], y=[w['breadth'] for _, w in rows], mode="markers",
                marker=dict(size=9, color=GAP, line=dict(color=CHROME["surface"], width=1.5)), showlegend=False,
                text=[f"{ref}:{w['start']:,}-{w['end']:,}<br>{', '.join(w.get('reasons', []))}" for ref, w in rows],
                hovertemplate="%{text}<br>breadth %{y:.1f}%<extra></extra>"), row=2, col=2)
            fig.update_xaxes(title_text="Window start (bp)", row=2, col=2)
            fig.update_yaxes(title_text="Window breadth (%)", row=2, col=2)
        else:
            _note(fig, "No window with low breadth and an extreme feature", 2, 2)
        _layout(fig, "Reference sequence analysis", 820)
        return fig

    def create_all_interactive_plots(self) -> Dict[str, str]:
        """Create all interactive plots and return file paths."""
        plot_files: Dict[str, str] = {}
        for key, method in (("dashboard", self.create_coverage_dashboard),
                            ("gap_analysis", self.create_detailed_gap_analysis),
                            ("quality_assessment", self.create_quality_assessment_plot),
                            ("reference_analysis", self.create_reference_analysis_plot)):
            try:
                plot_files[key] = method()
            except Exception as e:
                logging.error(f"Error creating {key} figure: {e}")
        return plot_files
