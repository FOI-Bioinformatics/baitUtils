#!/usr/bin/env python3

"""
report_generator.py

The HTML report of the evaluate command: a self-contained page with a
table of contents, numbered sections, stat tiles, captioned figures and
tables, all styled by the design system in baitUtils.design. Interactive
figures come from InteractivePlotter; plotly.js is loaded from a CDN or
embedded with offline_plots.
"""

import html
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import plotly.offline as pyo

from baitUtils import design
from baitUtils._version import __version__
from baitUtils.design import status_badge
from baitUtils.interactive_plots import InteractivePlotter
from baitUtils.quality_scorer import QualityScorer
from baitUtils.templates import load_template


def _e(value) -> str:
    """HTML-escape a value that comes from input files (sequence ids)."""
    return html.escape(str(value), quote=True)


# Colours of the quality categories, for the stat tile
CATEGORY_COLORS = design.CATEGORY_COLORS

PLOTLY_CDN = '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'


def _num(value, fmt: str = "{:,}") -> str:
    """Formatted number, or an en dash when missing."""
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "–"
    try:
        return fmt.format(value)
    except (ValueError, TypeError):
        return _e(value)


class InteractiveReportGenerator:
    """Generate the evaluate HTML report."""

    def __init__(
        self,
        coverage_stats: Dict[str, Any],
        gap_analysis: Dict[str, Any],
        quality_scores: Dict[str, Any],
        reference_analysis: Dict[str, Any],
        output_dir: Path,
        offline_plots: bool = False,
        mapper: str = "pblat",
        coverage_arrays: Optional[Dict[str, np.ndarray]] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ):
        """
        Args:
            coverage_stats: Output of CoverageAnalyzer.analyze
            gap_analysis: Output of GapAnalyzer.analyze
            quality_scores: QualityScore.to_dict()
            reference_analysis: Output of ReferenceAnalyzer.analyze
            output_dir: Directory of the report
            offline_plots: Embed plotly.js instead of loading it from a CDN
            mapper: Name of the mapper, for the header and the methodology
            coverage_arrays: Per-reference depth arrays, for the coverage profile
            metadata: Inputs and tool versions shown in the header
                (keys: oligos, reference, command, versions)
        """
        self.coverage_stats = coverage_stats or {}
        self.gap_analysis = gap_analysis or {}
        self.quality_scores = quality_scores or {}
        self.reference_analysis = reference_analysis or {}
        self.output_dir = Path(output_dir)
        self.offline_plots = offline_plots
        self.mapper = mapper
        self.coverage_arrays = coverage_arrays or {}
        self.metadata = metadata or {}

        self.report_title = "Coverage evaluation"
        self.generation_time = datetime.now().isoformat()
        self.html_components: List[str] = []
        self.sections: List[tuple] = []       # (id, title) for the table of contents
        self.figure_count = 0
        self.table_count = 0
        self.plotter = InteractivePlotter(self.coverage_stats, self.gap_analysis, self.quality_scores,
                                          self.reference_analysis, self.output_dir, coverage_arrays=self.coverage_arrays)

    # ------------------------------------------------------------ driver
    def generate_report(self) -> Path:
        """Generate the report and return its path."""
        logging.info("Generating HTML report...")
        self.html_components = []
        self.sections = []
        self.figure_count = 0
        self.table_count = 0

        self._create_header()
        self._create_executive_summary()
        self._create_coverage_overview()
        self._create_per_reference_section()
        self._create_gap_analysis_section()
        self._create_quality_assessment()
        self._create_reference_analysis()
        self._create_recommendations()
        self._create_detailed_statistics()
        self._create_data_exports()
        self._create_methodology()

        report_file = self.output_dir / "coverage_evaluation_report.html"
        report_file.write_text(self._assemble_html(), encoding="utf-8")
        logging.info(f"HTML report written to {report_file}")
        return report_file

    # ----------------------------------------------------------- helpers
    def _section(self, section_id: str, title: str, body: str) -> None:
        self.sections.append((section_id, title))
        number = len(self.sections)
        self.html_components.append(
            f'<section class="section" id="{section_id}">\n'
            f'<h2><span class="num">{number}</span>{_e(title)}</h2>\n{body}\n</section>')

    def _figure(self, plot_html: str, caption: str) -> str:
        self.figure_count += 1
        return (f'<figure><div class="plot">{plot_html}</div>'
                f'<figcaption><strong>Figure {self.figure_count}.</strong> {caption}</figcaption></figure>')

    def _table(self, caption: str, head: List[str], rows: List[List[str]], numeric: Optional[set] = None,
               css: str = "") -> str:
        self.table_count += 1
        numeric = numeric or set()
        head_html = "".join(f'<th class="num">{h}</th>' if i in numeric else f'<th>{h}</th>' for i, h in enumerate(head))
        body = []
        for row in rows:
            body.append("<tr>" + "".join(
                f'<td class="num">{c}</td>' if i in numeric else f'<td>{c}</td>' for i, c in enumerate(row)) + "</tr>")
        return (f'<div class="table-scroll"><table class="{css}"><caption><strong>Table {self.table_count}.</strong> '
                f'{caption}</caption><thead><tr>{head_html}</tr></thead><tbody>{"".join(body)}</tbody></table></div>')

    def _plot(self, fig, div_id: str) -> str:
        # The div id is prefixed so that it never collides with a section id.
        return fig.to_html(full_html=False, include_plotlyjs=False, div_id=f"fig-{div_id}", config=design.PLOTLY_CONFIG)

    @staticmethod
    def _category_for(score: float) -> str:
        """Category from the scorer's thresholds, for quality dicts without one."""
        thresholds = QualityScorer.DEFAULT_THRESHOLDS
        if score >= thresholds['excellent']:
            return "Excellent"
        if score >= thresholds['good']:
            return "Good"
        if score >= thresholds['fair']:
            return "Fair"
        return "Poor"

    @property
    def _params(self) -> Dict[str, Any]:
        return self.coverage_stats.get('parameters', {}) or {}

    # ---------------------------------------------------------- sections
    def _create_header(self) -> None:
        params = self._params
        versions = self.metadata.get('versions', {}) or {}
        meta = [
            ("Baits", _e(Path(str(self.metadata.get('oligos', ''))).name) if self.metadata.get('oligos') else None),
            ("Reference", _e(Path(str(self.metadata.get('reference', ''))).name) if self.metadata.get('reference') else None),
            ("Mapper", _e(self.mapper) + (f" {_e(versions['pblat'] if self.mapper == 'pblat' else versions.get('minimap2', ''))}"
                                         if versions.get(self.mapper) else "")),
            ("Minimum identity", f"{params.get('min_identity', 'n/a')}%"),
            ("Minimum aligned length", f"{params.get('min_length', 'n/a')} bp"),
            ("Minimum depth", f"{params.get('min_coverage', 'n/a')}x"),
            ("Target depth", f"{params.get('target_coverage', 'n/a')}x"),
            ("Generated", datetime.now().strftime('%Y-%m-%d %H:%M')),
        ]
        items = "".join(f"<div><dt>{k}</dt><dd>{v}</dd></div>" for k, v in meta if v)
        self.html_components.append(f"""
<header class="report-header">
  <p class="kicker">baitUtils {_e(__version__)} evaluate</p>
  <h1>Coverage evaluation of a bait set</h1>
  <p class="subtitle">Mapping, coverage, gaps and quality of one bait set against its reference.</p>
  <dl class="meta">{items}</dl>
</header>""")

    def _create_executive_summary(self) -> None:
        cs, ga = self.coverage_stats, self.gap_analysis
        breadth = cs.get('coverage_breadth', 0.0)
        quality_score = float(self.quality_scores.get('overall_score', 0) or 0)
        category = self.quality_scores.get('category') or self._category_for(quality_score)
        tiles = [
            ("Coverage breadth", f"{breadth:.1f}<small>%</small>",
             f"{cs.get('coverage_breadth_non_n', breadth):.1f}% of non-N bases" if 'coverage_breadth_non_n' in cs else
             f"{_num(cs.get('covered_bases'))} of {_num(cs.get('total_bases') or cs.get('reference_length'))} bp"),
            ("Mean depth", f"{cs.get('mean_depth', 0):.1f}<small>x</small>",
             f"median {cs.get('median_depth', 0):.1f}x, target {self._params.get('target_coverage', 'n/a')}x"),
            ("Mapped baits", f"{_num(cs.get('mapped_oligos', 0))}<small> of {_num(cs.get('total_oligos', 0))}</small>",
             f"{cs.get('mapping_efficiency', 0):.1f}% mapping efficiency"),
            ("Gaps", f"{_num(ga.get('total_gaps', 0))}",
             f"{design.fmt_bp(ga.get('total_gap_length', 0))} uncovered, largest {design.fmt_bp(ga.get('max_gap_size', 0))}"
             if ga.get('total_gaps') else "none at or above the minimum gap size"),
        ]
        tiles_html = "".join(
            f'<div class="kpi"><div class="label">{label}</div><div class="value">{value}</div><div class="sub">{sub}</div></div>'
            for label, value, sub in tiles)
        tiles_html += (f'<div class="kpi"><div class="label">Overall quality</div>'
                       f'<div class="value">{quality_score:.2f}<small> / 1</small></div>'
                       f'<div class="sub">{status_badge(category)}</div></div>')
        self._section("summary", "Summary", f'<div class="kpis">{tiles_html}</div>')

    def _create_coverage_overview(self) -> None:
        parts = []
        profile = self.plotter.coverage_profile_figure()
        if profile is not None:
            n = len(self.coverage_arrays)
            shown = min(n, 8)
            parts.append(self._figure(self._plot(profile, "coverage-profile"),
                                      f"Depth along {'each reference' if n <= shown else f'the {shown} longest of {n} references'}: "
                                      f"mean number of aligned baits per bin, gaps shaded, the dotted line at the minimum depth "
                                      f"of {self._params.get('min_coverage', 1)}x. Drag to zoom; double-click to reset."))
        else:
            parts.append('<p class="note">No per-base coverage was available for a coverage profile.</p>')
        parts.append(self._figure(self._plot(self.plotter.cumulative_depth_figure(), "cumulative-depth"),
                                  "Share of reference bases with at least a given depth. The curve reaches the breadth at the "
                                  "minimum depth and falls to the share of bases at or above the target depth."))
        self._section("coverage", "Coverage", "".join(parts))

    def _create_per_reference_section(self) -> None:
        per_ref = self.coverage_stats.get('per_reference', {}) or {}
        if not per_ref:
            self._section("references", "Per-reference results", '<p class="note">No per-reference statistics.</p>')
            return
        features = self.reference_analysis.get('sequence_features', {}) or {}
        difficulty = self.reference_analysis.get('challenging_regions', {}) or {}
        refs = design.select_references(per_ref, 50)
        rows = []
        for ref in refs:
            s = per_ref[ref]
            f = features.get(ref, {})
            d = difficulty.get(ref, {})
            rows.append([
                _e(ref), _num(s.get('length', 0)), _num(f.get('gc_content'), "{:.1f}"),
                _num(s.get('n_bases', 0)), _num(s.get('coverage_breadth'), "{:.1f}"),
                _num(s.get('coverage_breadth_non_n'), "{:.1f}"), _num(s.get('mean_depth'), "{:.2f}"),
                _num(s.get('max_depth')), _num(s.get('gaps', 0)),
                f"{_num(s.get('plus_hits', 0))} / {_num(s.get('minus_hits', 0))}",
                _e(d.get('difficulty_level', '–')),
            ])
        note = f" The {len(refs)} longest of {len(per_ref)} references are listed." if len(per_ref) > len(refs) else ""
        table = self._table(
            f"Coverage by reference sequence.{note} Gaps are runs below the minimum depth of any length; "
            f"hits are counted on the plus and minus strand.",
            ["Reference", "Length (bp)", "GC (%)", "N bases", "Breadth (%)", "Breadth, non-N (%)", "Mean depth",
             "Max depth", "Gaps", "Hits + / −", "Difficulty"],
            rows, numeric={1, 2, 3, 4, 5, 6, 7, 8, 9}, css="compact")
        self._section("references", "Per-reference results", table)

    def _create_gap_analysis_section(self) -> None:
        ga = self.gap_analysis
        parts = []
        total = ga.get('total_gaps', 0)
        if total:
            parts.append(f'<p class="lead">{_num(total)} gap{"s" if total != 1 else ""} of at least the minimum gap size, '
                         f'{design.fmt_bp(ga.get("total_gap_length", 0))} in total '
                         f'({ga.get("gap_percentage", 0):.1f}% of the reference); mean {design.fmt_bp(ga.get("mean_gap_size", 0))}, '
                         f'median {design.fmt_bp(ga.get("median_gap_size", 0))}, largest {design.fmt_bp(ga.get("max_gap_size", 0))}.</p>')
            parts.append(self._figure(self._plot(self.plotter.gap_figure(), "gaps"),
                                      "Gaps by start position and length (colour: reference), the largest gaps, gap length "
                                      "against the GC content of the gap sequence, and the length distribution."))
            parts.append(self._create_largest_gaps_table())
            dist = ga.get('size_distribution', {}) or {}
            if dist:
                rows = [[_e(k), _num(v)] for k, v in dist.items()]
                parts.append(self._table("Gaps by length class (bp).", ["Length class", "Gaps"], rows, numeric={1}, css="kv compact"))
        else:
            parts.append('<p class="lead">No gap at or above the minimum gap size.</p>')
        self._section("gaps", "Gaps", "".join(parts))

    def _create_largest_gaps_table(self) -> str:
        largest = (self.gap_analysis.get('largest_gaps', []) or [])[:15]
        if not largest:
            return ""
        features = {(f.get('chromosome'), f.get('start')): f for f in self.gap_analysis.get('gap_features', []) or []}
        rows = []
        for i, gap in enumerate(largest, 1):
            f = features.get((gap.get('chromosome'), gap.get('start')), {})
            rows.append([str(i), _e(gap.get('chromosome', 'Unknown')), _num(gap.get('start', 0)), _num(gap.get('end', 0)),
                         _num(gap.get('length', 0)), _num(f.get('gc_content'), "{:.1f}"),
                         _num(f.get('n_content'), "{:.1f}")])
        return self._table("Largest gaps (0-based, half-open coordinates).",
                           ["#", "Reference", "Start", "End", "Length (bp)", "GC (%)", "N (%)"], rows,
                           numeric={0, 2, 3, 4, 5, 6}, css="compact")

    def _create_quality_assessment(self) -> None:
        scores = self.quality_scores.get('component_scores', {}) or {}
        if not scores:
            self._section("quality", "Quality score", '<p class="note">No quality score was computed.</p>')
            return
        overall = float(self.quality_scores.get('overall_score', 0) or 0)
        category = self.quality_scores.get('category') or self._category_for(overall)
        weights = QualityScorer.DEFAULT_WEIGHTS
        rows = [[_e(name.replace('_', ' ').capitalize()), f"{weights.get(name, 0):.0%}", f"{score:.2f}",
                 f"{self.quality_scores.get('weighted_scores', {}).get(name, score * weights.get(name, 0)):.3f}"]
                for name, score in scores.items()]
        rows.append(["<strong>Overall</strong>", "", "", f"<strong>{overall:.3f}</strong>"])
        body = (f'<p class="lead">Overall score {overall:.2f} of 1: {status_badge(category)}.</p>'
                + self._figure(self._plot(self.plotter.quality_figure(), "quality"),
                               "Component scores with their weights, the observed breadth, depth, uniformity and mapping "
                               "efficiency relative to the design targets, and the overall score against the category thresholds.")
                + self._table("Quality score components.", ["Component", "Weight", "Score (0-1)", "Weighted"], rows,
                              numeric={1, 2, 3}, css="kv compact"))
        self._section("quality", "Quality score", body)

    def _create_reference_analysis(self) -> None:
        if not self.reference_analysis:
            self._section("reference", "Reference sequence analysis", '<p class="note">No reference analysis.</p>')
            return
        parts = [self._figure(self._plot(self.plotter.reference_figure(), "reference"),
                              "Per-reference breadth against GC and repeat content (marker size: length), mean window "
                              "breadth by GC decile, and windows with low breadth and an extreme sequence feature.")]
        parts.append(self._create_correlation_tables())
        parts.append(self._create_challenging_windows_table())
        self._section("reference", "Reference sequence analysis", "".join(parts))

    def _create_correlation_tables(self) -> str:
        corr = self.reference_analysis.get('coverage_correlations') or {}
        if not corr.get('applicable'):
            return f'<p class="note">Feature-coverage correlations not applicable: {_e(corr.get("reason", "too few windows"))}.</p>'
        rows = []
        for name, entry in (corr.get('features') or {}).items():
            cells = [_e(name.replace('_', ' '))]
            for target in ('breadth', 'mean_depth'):
                rho = entry.get(f'rho_{target}', float('nan'))
                p = entry.get(f'p_{target}', float('nan'))
                cells += ["–", "–"] if rho != rho else [f"{rho:+.3f}", f"{p:.3g}"]
            rows.append(cells)
        table = self._table(
            f"Spearman rank correlation between window features and coverage over {corr.get('windows', 0):,} windows of "
            f"{corr.get('window_size', 0):,} bp. A negative rho with breadth means that windows high in the feature are "
            f"covered less.",
            ["Feature", "rho (breadth)", "p", "rho (mean depth)", "p"], rows, numeric={1, 2, 3, 4}, css="compact")
        gc_rows = [[f"{b['gc_min']}-{b['gc_max']}%", _num(b['windows']), f"{b['mean_breadth']:.1f}", f"{b['mean_depth']:.2f}"]
                   for b in corr.get('breadth_by_gc', []) or []]
        gc_table = self._table("Window breadth and depth by GC decile.", ["GC", "Windows", "Mean breadth (%)", "Mean depth"],
                               gc_rows, numeric={1, 2, 3}, css="compact kv") if gc_rows else ""
        return table + gc_table

    def _create_challenging_windows_table(self) -> str:
        windows = self.reference_analysis.get('challenging_windows') or {}
        rows = [[_e(ref), _num(w['start']), _num(w['end']), f"{w['breadth']:.1f}", _e(', '.join(w.get('reasons', [])))]
                for ref, entries in windows.items() for w in entries]
        if not rows:
            return '<p class="note">No window with low breadth and an extreme sequence feature.</p>'
        return self._table("Windows with low breadth and an extreme sequence feature (first 200).",
                           ["Reference", "Start", "End", "Breadth (%)", "Features"], rows[:200], numeric={1, 2, 3}, css="compact")

    def _create_recommendations(self) -> None:
        items = list(self.gap_analysis.get('suggestions', []) or []) + list(self.quality_scores.get('recommendations', []) or [])
        if not items:
            items = ["Coverage analysis complete. No specific recommendations at this time."]
        self._section("recommendations", "Recommendations",
                      '<div class="recommendations"><ol>' + "".join(f"<li>{_e(r)}</li>" for r in items) + "</ol></div>")

    def _create_detailed_statistics(self) -> None:
        body = f'<details><summary>All coverage statistics</summary><div class="details-body">{self._create_statistics_table()}</div></details>'
        self._section("statistics", "Detailed statistics", body)

    def _create_statistics_table(self) -> str:
        cs, ga = self.coverage_stats, self.gap_analysis
        groups = [
            ("Coverage", [
                ("Coverage breadth", _num(cs.get('coverage_breadth'), "{:.2f}%")),
                ("Coverage breadth, non-N bases", _num(cs.get('coverage_breadth_non_n'), "{:.2f}%")),
                ("Breadth at the target depth", _num(cs.get('target_coverage_breadth'), "{:.2f}%")),
                ("Mean depth", _num(cs.get('mean_depth'), "{:.2f}x")),
                ("Median depth", _num(cs.get('median_depth'), "{:.2f}x")),
                ("Depth standard deviation", _num(cs.get('depth_std'), "{:.2f}x")),
                ("Maximum depth", _num(cs.get('max_depth'), "{:,}x")),
                ("Coefficient of variation (covered bases)", _num(cs.get('coverage_cv'), "{:.3f}")),
                ("Gini coefficient (covered bases)", _num(cs.get('coverage_gini'), "{:.3f}")),
            ]),
            ("Mapping", [
                ("Baits in the input", _num(cs.get('total_oligos'))),
                ("Baits with a passing hit", _num(cs.get('mapped_oligos'))),
                ("Mapping efficiency", _num(cs.get('mapping_efficiency'), "{:.1f}%")),
                ("Passing hits", _num(cs.get('total_mappings'))),
                ("Hits on the plus / minus strand",
                 f"{_num((cs.get('strand_counts') or {}).get('+'))} / {_num((cs.get('strand_counts') or {}).get('-'))}"),
            ]),
            ("Reference", [
                ("Reference sequences", _num(cs.get('reference_sequences'))),
                ("Reference length", _num(cs.get('reference_length'), "{:,} bp")),
                ("N bases", _num(cs.get('n_bases'), "{:,} bp")),
                ("Covered bases", _num(cs.get('covered_bases'), "{:,} bp")),
                ("Uncovered bases", _num(cs.get('uncovered_bases'), "{:,} bp")),
            ]),
            ("Gaps", [
                ("Gaps", _num(ga.get('total_gaps'))),
                ("Total gap length", _num(ga.get('total_gap_length'), "{:,} bp")),
                ("Share of the reference in gaps", _num(ga.get('gap_percentage'), "{:.2f}%")),
                ("Mean gap length", _num(ga.get('mean_gap_size'), "{:,.0f} bp")),
                ("Median gap length", _num(ga.get('median_gap_size'), "{:,.0f} bp")),
                ("Largest gap", _num(ga.get('max_gap_size'), "{:,} bp")),
            ]),
        ]
        rows = []
        for group, items in groups:
            rows.append([f"<strong>{group}</strong>", ""])
            rows.extend([[label, value] for label, value in items])
        return self._table("All statistics of the evaluation.", ["Metric", "Value"], rows, numeric={1}, css="kv compact")

    def _create_data_exports(self) -> None:
        body = """
<div class="exports">
  <a href="evaluation.json"><strong>evaluation.json</strong><span>every result, with tool versions and arguments</span></a>
  <a href="data/coverage.bedgraph"><strong>coverage.bedgraph</strong><span>run-length depth, 0-based half-open</span></a>
  <a href="data/gap_regions.bed"><strong>gap_regions.bed</strong><span>runs below the minimum depth</span></a>
  <a href="hits.tsv"><strong>hits.tsv</strong><span>per-bait hit table</span></a>
  <a href="plots/"><strong>plots/</strong><span>static figures</span></a>
</div>"""
        self._section("exports", "Data files", body)

    def _create_methodology(self) -> None:
        params = self._params
        weights = QualityScorer.DEFAULT_WEIGHTS
        targets = QualityScorer.DEFAULT_BENCHMARKS
        thresholds = QualityScorer.DEFAULT_THRESHOLDS
        weight_items = "".join(
            f"<li><strong>{name.replace('_', ' ').title()}</strong>: {weight:.0%}</li>" for name, weight in weights.items())
        body = f"""
<h4>Mapping</h4>
<p>Oligos were mapped to the reference with {_e(self.mapper)}. Hits count when their BLAT identity is
at least {params.get('min_identity', 'N/A')}% and their aligned length is at least {params.get('min_length', 'N/A')} bp.</p>
<h4>Coverage</h4>
<p>Depth is the number of aligned blocks over each reference position. Coverage breadth is the percentage of
positions with depth at least {params.get('min_coverage', 'N/A')}; it is also reported over the bases that are not N.</p>
<h4>Gaps</h4>
<p>Gaps are runs of positions below the minimum depth, at or above the minimum gap size. Each gap is characterised by
GC content, Shannon entropy, duplicated 6-mers, homopolymer runs and N content.</p>
<h4>Quality score</h4>
<p>The overall score is a weighted sum of five component scores in 0-1:</p>
<ul>{weight_items}</ul>
<p>Targets: breadth {targets['coverage_breadth']:.0f}%, mean depth {targets['mean_depth']:.0f}x, uniformity (1 - Gini)
{targets['depth_uniformity']:.1f}, mapping efficiency {targets['mapping_efficiency']:.0f}%, gaps at most
{targets['gap_percentage']:.0f}% of the reference and {targets['largest_gap']:,} bp for the largest gap. The reference
difficulty component is breadth relative to the breadth expected for the reference's difficulty.
Categories: Excellent &gt;= {thresholds['excellent']}, Good &gt;= {thresholds['good']}, Fair &gt;= {thresholds['fair']},
otherwise Poor.</p>
"""
        self._section("methods", "Methods", body)

    # ---------------------------------------------------------- assembly
    def _assemble_html(self) -> str:
        plotly_script = f"<script>{pyo.get_plotlyjs()}</script>" if self.offline_plots else PLOTLY_CDN
        toc = "".join(f'<li><a href="#{sid}">{i}. {_e(title)}</a></li>' for i, (sid, title) in enumerate(self.sections, 1))
        versions = self.metadata.get('versions', {}) or {}
        footer_bits = [f"baitUtils {_e(__version__)}"]
        for tool in ("python", "numpy", "scipy", "biopython"):
            if versions.get(tool):
                footer_bits.append(f"{tool} {_e(versions[tool])}")
        command = self.metadata.get('command')
        command_html = f"<br><code>{_e(command)}</code>" if command else ""
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Coverage evaluation - baitUtils</title>
{plotly_script}
<style>
{design.css_variables()}
{self._get_css_styles()}
</style>
</head>
<body>
<div class="layout">
  <nav class="toc" aria-label="Contents">
    <div class="brand">baitUtils<small>evaluate report</small></div>
    <ol>{toc}</ol>
  </nav>
  <main class="report">
    {''.join(self.html_components)}
    <footer class="report-footer">Generated {datetime.now().strftime('%Y-%m-%d %H:%M')} with {' · '.join(footer_bits)}.{command_html}</footer>
  </main>
</div>
<script>{self._get_js_scripts()}</script>
</body>
</html>
"""

    def _get_css_styles(self) -> str:
        return load_template("report.css")

    def _get_js_scripts(self) -> str:
        return load_template("report.js")
