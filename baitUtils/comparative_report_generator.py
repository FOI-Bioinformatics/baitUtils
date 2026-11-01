#!/usr/bin/env python3

"""
comparative_report_generator.py

The HTML report of the compare command, on the same design system as the
evaluate report: a comparison table with the best value per metric marked,
the ranking, the statistical tests as tables with their sample sizes, gap
overlap, the figures, recommendations and the methods.
"""

import base64
import html
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from baitUtils import design
from baitUtils._version import __version__
from baitUtils.comparative_analyzer import ComparativeAnalyzer
from baitUtils.comparative_visualizations import ComparativeVisualizer, set_colors, MATRIX_METRICS
from baitUtils.differential_analysis import DifferentialAnalyzer
from baitUtils.quality_scorer import QualityScorer
from baitUtils.templates import load_template


def _e(value) -> str:
    """HTML-escape a value that comes from user input (set names, sequence ids)."""
    return html.escape(str(value), quote=True)


IMAGE_TYPES = {'.png': 'image/png', '.svg': 'image/svg+xml', '.jpg': 'image/jpeg', '.jpeg': 'image/jpeg'}
SMALL_N = 5   # below this many paired observations a test is reported with a caveat
PLOTLY_CDN = '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'


class BestValue(str):
    """Table cell text marked as the best value of its column."""


class ComparativeReportGenerator:
    """Generate the compare HTML report."""

    def __init__(self, analyzer: ComparativeAnalyzer, output_dir: Path,
                 differential_analyzer: Optional[DifferentialAnalyzer] = None,
                 visualizer: Optional[ComparativeVisualizer] = None,
                 plots: Optional[Dict[str, str]] = None,
                 metadata: Optional[Dict[str, Any]] = None):
        """
        Args:
            analyzer: Comparative analyzer with the oligo set results
            output_dir: Directory of the report
            differential_analyzer: Statistical tests, or None when disabled
            visualizer: Visualizer configured with the user's plot format and dpi
            plots: Figures already generated with that visualizer ({name: path});
                when given, they are not generated again
            metadata: Inputs, parameters and versions shown in the header
                (keys: reference, mapper, parameters, versions, command)
        """
        self.analyzer = analyzer
        self.differential_analyzer = differential_analyzer
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.report_file = self.output_dir / "comparative_analysis_report.html"
        self.visualizer = visualizer or ComparativeVisualizer(output_dir)
        self.plots = plots
        self.metadata = metadata or {}
        self.sections: List[tuple] = []
        self.figure_count = 0
        self.table_count = 0

    # ------------------------------------------------------------ driver
    def generate_report(self) -> str:
        logging.info("Generating comparative analysis report...")
        plots = self.plots
        if plots is None:
            plots = self.visualizer.generate_all_comparative_plots(self.analyzer, self.differential_analyzer)
        self.sections, self.figure_count, self.table_count = [], 0, 0
        content = self._generate_html_content(plots)
        self.report_file.write_text(content, encoding="utf-8")
        logging.info(f"Comparative analysis report saved to {self.report_file}")
        return str(self.report_file)

    # ----------------------------------------------------------- helpers
    @property
    def _names(self) -> List[str]:
        return [s.name for s in self.analyzer.oligo_sets]

    @property
    def _colors(self) -> Dict[str, str]:
        return set_colors(self._names)

    def _set_label(self, name: str) -> str:
        return f'<span class="badge">{design.swatch(self._colors.get(name, design.CHROME["ink_muted"]))}{_e(name)}</span>'

    def _section(self, section_id: str, title: str, body: str) -> str:
        self.sections.append((section_id, title))
        return (f'<section class="section" id="{section_id}">\n'
                f'<h2><span class="num">{len(self.sections)}</span>{_e(title)}</h2>\n{body}\n</section>')

    def _figure(self, inner: str, caption: str) -> str:
        self.figure_count += 1
        return (f'<figure><div class="plot">{inner}</div>'
                f'<figcaption><strong>Figure {self.figure_count}.</strong> {caption}</figcaption></figure>')

    def _table(self, caption: str, head: List[str], rows: List[List[str]], numeric: Optional[set] = None,
               css: str = "compact", numbered: bool = True) -> str:
        numeric = numeric or set()
        if numbered:
            self.table_count += 1
            caption = f'<strong>Table {self.table_count}.</strong> {caption}'
        head_html = "".join(f'<th class="num">{h}</th>' if i in numeric else f'<th>{h}</th>' for i, h in enumerate(head))

        def cell(i, c):
            classes = []
            if i in numeric:
                classes.append("num")
            if isinstance(c, BestValue):
                classes.append("best")
            attr = f' class="{" ".join(classes)}"' if classes else ""
            return f'<td{attr}>{c}</td>'

        body = "".join("<tr>" + "".join(cell(i, c) for i, c in enumerate(row)) + "</tr>" for row in rows)
        return (f'<div class="table-scroll"><table class="{css}"><caption>{caption}</caption>'
                f'<thead><tr>{head_html}</tr></thead><tbody>{body}</tbody></table></div>')

    def _relative(self, path) -> str:
        """Path of an output file relative to the report, for links."""
        try:
            return Path(path).resolve().relative_to(self.output_dir.resolve()).as_posix()
        except ValueError:
            return Path(path).name

    # ----------------------------------------------------------- content
    def _generate_html_content(self, plots: Dict[str, str]) -> str:
        comparison_matrix = self.analyzer.generate_comparison_matrix()
        ranking = self.analyzer.generate_ranking()
        gap_analysis = self.analyzer.analyze_gap_overlap()
        sections = [
            self._generate_executive_summary(),
            self._generate_comparison_overview(comparison_matrix),
            self._generate_performance_ranking(ranking),
            self._generate_coverage_profiles(),
            self._generate_statistical_analysis() if self.differential_analyzer else "",
            self._generate_gap_overlap_analysis(gap_analysis),
            self._generate_visualizations_section(plots),
            self._generate_recommendations(),
            self._generate_detailed_analysis(),
            self._generate_methodology_section(),
        ]
        body = "\n".join(s for s in sections if s)
        toc = "".join(f'<li><a href="#{sid}">{i}. {_e(title)}</a></li>' for i, (sid, title) in enumerate(self.sections, 1))
        return (self._get_html_header()
                + f'<nav class="toc" aria-label="Contents"><div class="brand">baitUtils<small>compare report</small></div><ol>{toc}</ol></nav>\n'
                + '<main class="report">\n' + self._header_block() + body + self._footer_block() + '\n</main>\n'
                + self._get_html_footer())

    def _get_html_header(self) -> str:
        """Document head and the opening of the layout."""
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Comparison of bait sets - baitUtils</title>
{PLOTLY_CDN}
<style>
{design.css_variables()}
{load_template("report.css")}
</style>
</head>
<body>
<div class="layout">
"""

    def _header_block(self) -> str:
        params = self.metadata.get('parameters', {}) or {}
        versions = self.metadata.get('versions', {}) or {}
        mapper = self.metadata.get('mapper', getattr(self.analyzer, 'mapper', 'pblat'))
        meta = [
            ("Reference", _e(Path(str(self.metadata.get('reference', ''))).name) if self.metadata.get('reference') else None),
            ("Sets", ", ".join(_e(n) for n in self._names)),
            ("Mapper", _e(mapper) + (f" {_e(versions[mapper])}" if versions.get(mapper) else "")),
            ("Minimum identity", f"{params['min_identity']}%" if 'min_identity' in params else None),
            ("Minimum aligned length", f"{params['min_length']} bp" if 'min_length' in params else None),
            ("Minimum depth", f"{params['min_coverage']}x" if 'min_coverage' in params else None),
            ("Minimum gap size", f"{params['min_gap_size']} bp" if 'min_gap_size' in params else None),
            ("Generated", datetime.now().strftime('%Y-%m-%d %H:%M')),
        ]
        items = "".join(f"<div><dt>{k}</dt><dd>{v}</dd></div>" for k, v in meta if v)
        sets_rows = [[self._set_label(s.name), f"<code>{_e(Path(str(s.file_path)).name)}</code>",
                      f"{s.coverage_stats.get('total_oligos', 0):,}", f"{s.coverage_stats.get('mapped_oligos', 0):,}"]
                     for s in self.analyzer.oligo_sets]
        sets_table = self._table("The bait sets compared; each set keeps its colour in every figure.",
                                 ["Set", "File", "Baits", "Mapped"], sets_rows, numeric={2, 3}, css="compact kv",
                                 numbered=False)
        return f"""
<header class="report-header">
  <p class="kicker">baitUtils {_e(__version__)} compare</p>
  <h1>Comparison of {len(self._names)} bait sets</h1>
  <p class="subtitle">Coverage, gaps, quality and statistical tests of several bait sets against one reference.</p>
  <dl class="meta">{items}</dl>
</header>
{sets_table}
"""

    def _footer_block(self) -> str:
        versions = self.metadata.get('versions', {}) or {}
        bits = [f"baitUtils {_e(__version__)}"] + [f"{t} {_e(versions[t])}" for t in ("python", "scipy") if versions.get(t)]
        command = self.metadata.get('command')
        command_html = f"<br><code>{_e(command)}</code>" if command else ""
        return (f'<footer class="report-footer">Generated {datetime.now().strftime("%Y-%m-%d %H:%M")} with '
                f'{" · ".join(bits)}.{command_html}</footer>')

    def _get_html_footer(self) -> str:
        return f"</div>\n<script>{load_template('report.js')}</script>\n</body>\n</html>\n"

    # ----------------------------------------------------------- sections
    def _generate_executive_summary(self) -> str:
        best = self.analyzer.identify_best_performer('quality_score')
        ranking = self.analyzer.generate_ranking()
        sets = self.analyzer.oligo_sets
        breadths = [s.coverage_stats.get('coverage_breadth', 0) for s in sets]
        gaps = [s.gap_analysis.get('total_gaps', 0) for s in sets]
        effs = [s.coverage_stats.get('mapping_efficiency', 0) for s in sets]
        tiles = [
            ("Sets compared", f"{len(sets)}", ""),
            ("Best quality score", f"{best.quality_score.overall_score:.2f}<small> / 1</small>",
             f"{self._set_label(best.name)} {design.status_badge(best.quality_score.category.value)}"),
            ("Coverage breadth", f"{min(breadths):.1f}<small>%</small> – {max(breadths):.1f}<small>%</small>", "range across sets"),
            ("Gaps", f"{min(gaps):,} – {max(gaps):,}", "range across sets"),
            ("Mapping efficiency", f"{min(effs):.1f}<small>%</small> – {max(effs):.1f}<small>%</small>", "range across sets"),
        ]
        tiles_html = "".join(
            f'<div class="kpi"><div class="label">{label}</div><div class="value">{value}</div><div class="sub">{sub}</div></div>'
            for label, value, sub in tiles)
        items = "".join(
            f'<li>{self._set_label(name)}: composite score {score:.2f}, quality '
            f'{next(s for s in sets if s.name == name).quality_score.overall_score:.2f} '
            f'({_e(next(s for s in sets if s.name == name).quality_score.category.value)})</li>'
            for name, score in ranking)
        body = (f'<div class="kpis">{tiles_html}</div>'
                f'<p><strong>{_e(best.name)}</strong> has the highest quality score, {best.quality_score.overall_score:.2f} of 1 '
                f'(Grade {_e(best.quality_score.category.value)}); breadth ranges from {min(breadths):.1f}% to {max(breadths):.1f}%, '
                f'gap counts from {min(gaps):,} to {max(gaps):,}.</p>'
                f'<h3>Ranking by composite score</h3><ol>{items}</ol>')
        return self._section("summary", "Summary", body)

    def _generate_comparison_overview(self, comparison_matrix) -> str:
        if comparison_matrix is None:
            return ""
        columns = [m for m in MATRIX_METRICS if m[0] in comparison_matrix.columns]
        head = ["Set", "Grade"] + [label for _, label, _, _ in columns]
        rows = []
        best_per_column = {}
        for column, _label, higher_better, _fmt in columns:
            values = comparison_matrix[column].astype(float)
            best_per_column[column] = values.max() if higher_better else values.min()
        for _, row in comparison_matrix.iterrows():
            cells = [self._set_label(row['Name']), design.status_badge(row.get('Quality_Grade', ''))]
            for column, _label, _hb, fmt in columns:
                value = float(row[column])
                text = fmt.format(value)
                cells.append(BestValue(text) if np.isclose(value, best_per_column[column]) else text)
            rows.append(cells)
        table = self._table("Metrics of every set; the best value of each column is bold.", head, rows,
                            numeric=set(range(2, len(head))))
        return self._section("overview", "Comparison", table)

    def _generate_performance_ranking(self, ranking) -> str:
        if not ranking:
            return ""
        top = max(score for _, score in ranking) or 1.0
        rows = []
        for i, (name, score) in enumerate(ranking, 1):
            oligo_set = next(s for s in self.analyzer.oligo_sets if s.name == name)
            bar = f'<span class="inline-bar"><span style="width:{max(0.0, min(1.0, score / top)) * 100:.0f}%"></span></span>'
            rows.append([str(i), self._set_label(name), f"{bar}{score:.2f}",
                         f"{oligo_set.quality_score.overall_score:.2f}", design.status_badge(oligo_set.quality_score.category.value)])
        table = self._table("Composite score in 0-1: quality score 40 %, coverage breadth 25 %, mapping efficiency 20 %, "
                            "gap term 15 %.", ["#", "Set", "Composite score", "Quality score", "Grade"], rows, numeric={0, 3})
        return self._section("ranking", "Ranking", table)

    def _generate_coverage_profiles(self) -> str:
        fig = self.visualizer.coverage_profiles_figure(self.analyzer.oligo_sets)
        if fig is None:
            return ""
        inner = fig.to_html(full_html=False, include_plotlyjs=False, div_id="fig-profiles", config=design.PLOTLY_CONFIG)
        body = self._figure(inner, "Depth along the longest references for every set (mean depth per bin). "
                                   "Drag to zoom; double-click to reset.")
        return self._section("profiles", "Coverage profiles", body)

    def _generate_detailed_analysis(self) -> str:
        pairwise = self.analyzer.calculate_pairwise_comparisons()
        if not pairwise:
            return ""
        rows = [[f"{self._set_label(a)} − {self._set_label(b)}", f"{m.coverage_breadth_diff:+.1f}", f"{m.mean_depth_diff:+.2f}",
                 f"{m.gap_count_diff:+d}", f"{m.quality_score_diff:+.2f}", f"{m.mapping_efficiency_diff:+.1f}"]
                for (a, b), m in pairwise.items()]
        table = self._table("Pairwise differences (first set minus second).",
                            ["Pair", "Breadth (pp)", "Mean depth (x)", "Gaps", "Quality score", "Mapping efficiency (pp)"],
                            rows, numeric={1, 2, 3, 4, 5})
        return self._section("pairwise", "Pairwise differences", table)

    def _generate_statistical_analysis(self) -> str:
        da = self.differential_analyzer
        if not da:
            return ""
        alpha, method = da.significance_level, da.correction_method
        sets = self.analyzer.oligo_sets
        quality_tests = da.compare_quality_metrics(sets)
        parts = [f'<p>Paired tests across reference sequences (each reference is one observation); p-values adjusted for '
                 f'multiple comparisons ({_e(method)}) within each family; alpha = {alpha:g}. '
                 f'Labels: {_e(da.significance_legend())}.</p>']
        rows = [self._test_row(metric.replace('_', ' ').capitalize(), test) for metric, test in quality_tests.items()]
        parts.append(self._table("Per-reference metrics.", self._test_head(), rows, numeric={2, 3, 4, 5, 6}))
        small = [t for t in quality_tests.values() if t.applicable and t.n < SMALL_N]
        if small:
            parts.append(f'<p class="note">With fewer than {SMALL_N} references the paired tests have little power and the '
                         f'effect sizes are unstable; read them as descriptive.</p>')
        if len(sets) == 2:
            dist = da.compare_coverage_distributions(sets[0], sets[1])
            rows = [self._test_row(label, test) for label, test in (
                ("Depth distribution", dist.ks_test), ("Depth median", dist.mann_whitney_test),
                ("Depth variability", dist.levene_test))]
            parts.append(self._table(f"Mean depth per {dist.window_size:,} bp window, {_e(sets[0].name)} against {_e(sets[1].name)}.",
                                     self._test_head(), rows, numeric={2, 3, 4, 5, 6}))
            summary_rows = [[self._set_label(name), f"{s['mean']:.2f}", f"{s['median']:.2f}", f"{s['std']:.2f}",
                             f"{s['q25']:.2f}", f"{s['q75']:.2f}", f"{s['windows']:,}"] for name, s in dist.summary_stats.items()]
            parts.append(self._table("Window depth summary.", ["Set", "Mean", "Median", "SD", "Q1", "Q3", "Windows"],
                                     summary_rows, numeric={1, 2, 3, 4, 5, 6}))
            identity = da.compare_oligo_identity(sets[0], sets[1])
            gap_test = da.analyze_gap_patterns(sets[0], sets[1])['gap_size_distribution']
            parts.append(self._table("Best-hit identity per bait and gap lengths.", self._test_head(),
                                     [self._test_row("Best-hit identity per bait", identity), self._test_row("Gap lengths", gap_test)],
                                     numeric={2, 3, 4, 5, 6}))
        return self._section("statistics", "Statistical tests", "".join(parts))

    @staticmethod
    def _test_head() -> List[str]:
        return ["Quantity", "Test", "n", "Statistic", "p", "Adjusted p", "Effect size", "Label", "Result"]

    @staticmethod
    def _test_row(title: str, test) -> List[str]:
        if not test.applicable:
            return [_e(title), _e(test.test_name), f"{test.n}", "–", "–", "–", "–", "n/a", _e(test.interpretation)]
        effect = "–" if test.effect_size != test.effect_size else f"{test.effect_size:.3f}"
        adjusted = "–" if test.p_adjusted is None else f"{test.p_adjusted:.4f}"
        return [_e(title), f"{_e(test.test_name)}", f"{test.n}", f"{test.statistic:.3f}", f"{test.p_value:.4f}", adjusted,
                effect, _e(test.significance_level), _e(test.interpretation)]

    @staticmethod
    def _render_test(title: str, test, alpha: float) -> str:
        """One test as a short HTML block (kept for callers of the previous API)."""
        if not test.applicable:
            return (f'<div class="callout muted"><h4>{_e(title)}</h4><p><strong>Test:</strong> {_e(test.test_name)}</p>'
                    f'<p><strong>Result:</strong> {_e(test.interpretation)}</p></div>')
        adjusted = f" (adjusted {test.p_adjusted:.4f})" if test.p_adjusted is not None else ""
        effect = "n/a" if test.effect_size != test.effect_size else f"{test.effect_size:.3f}"
        significant = test.p_reported <= alpha
        return (f'<div class="callout{"" if significant else " muted"}"><h4>{_e(title)}</h4>'
                f'<p><strong>Test:</strong> {_e(test.test_name)} (n = {test.n})</p>'
                f'<p><strong>P-value:</strong> {test.p_value:.4f}{adjusted} {_e(test.significance_level)}</p>'
                f'<p><strong>Effect size:</strong> {effect}</p>'
                f'<p><strong>Interpretation:</strong> {_e(test.interpretation)}</p></div>')

    def _generate_gap_overlap_analysis(self, gap_analysis) -> str:
        if not gap_analysis:
            return ""
        rows = []
        for set_name, gaps in (gap_analysis.get('unique_gaps', {}) or {}).items():
            examples = "; ".join(f"{_e(g['chromosome'])}:{g['start']:,}-{g['end']:,}" for g in gaps[:3])
            if len(gaps) > 3:
                examples += f" (+{len(gaps) - 3} more)"
            rows.append([self._set_label(set_name), f"{len(gaps):,}", examples])
        overlap_rows = []
        for pair, overlaps in (gap_analysis.get('overlapping_gaps', {}) or {}).items():
            total = sum(o['length'] for o in overlaps)
            overlap_rows.append([_e(pair.replace('_vs_', ' and ')), f"{len(overlaps):,}", design.fmt_bp(total)])
        body = (f'<p>{gap_analysis.get("total_unique_regions", 0):,} distinct gap regions across all sets.</p>'
                + self._table("Gaps present in one set only.", ["Set", "Unique gaps", "Examples"], rows, numeric={1})
                + (self._table("Gap regions shared by two sets.", ["Sets", "Shared gaps", "Shared length"], overlap_rows, numeric={1, 2})
                   if overlap_rows else ""))
        return self._section("gaps", "Gap overlap", body)

    def _generate_recommendations(self) -> str:
        best = self.analyzer.identify_best_performer('quality_score')
        ranking = self.analyzer.generate_ranking()
        sets = self.analyzer.oligo_sets
        items = [f"<strong>{_e(best.name)}</strong> has the highest quality score ({best.quality_score.overall_score:.2f}, "
                 f"{_e(best.quality_score.category.value)}): breadth {best.coverage_stats.get('coverage_breadth', 0):.1f}%, "
                 f"{best.gap_analysis.get('total_gaps', 0):,} gaps, mapping efficiency {best.coverage_stats.get('mapping_efficiency', 0):.1f}%."]
        if len(ranking) > 1:
            second = next(s for s in sets if s.name == ranking[1][0])
            items.append(f"<strong>{_e(second.name)}</strong> ranks second by composite score (quality "
                         f"{second.quality_score.overall_score:.2f}, {_e(second.quality_score.category.value)}).")
        breadths = [s.coverage_stats.get('coverage_breadth', 0) for s in sets]
        if max(breadths) - min(breadths) > 20:
            items.append("Breadth differs by more than 20 percentage points between sets; the sets with low breadth lack "
                         "baits in regions the others cover.")
        if max(s.gap_analysis.get('total_gaps', 0) for s in sets) > 100:
            items.append("Sets with many gaps are candidates for baitUtils fill, which selects additional baits in uncovered regions.")
        if min(s.coverage_stats.get('mapping_efficiency', 0) for s in sets) < 80:
            items.append("Mapping efficiency below 80% in a set points to baits that do not match the reference; review their "
                         "design or the reference.")
        items.append("Gap regions shared by every set mark reference regions that none of the designs reaches.")
        return self._section("recommendations", "Recommendations",
                             '<div class="recommendations"><ol>' + "".join(f"<li>{i}</li>" for i in items) + "</ol></div>")

    def _generate_visualizations_section(self, plots: Dict[str, str]) -> str:
        titles = {
            'coverage_distributions': ("Coverage metrics by set",
                                       "Breadth, mean depth, uniformity (1 - Gini) and mapping efficiency of every set."),
            'gap_analysis': ("Gap metrics by set", "Gap count, share of the reference in gaps, mean and largest gap of every set."),
            'performance_heatmap': ("Comparison matrix", "Every metric of every set; darker cells are better within a column."),
            'statistical_comparison': ("Statistical comparison",
                                       "Adjusted p-values against alpha, effect sizes, the distribution of window depths and the "
                                       "paired breadth per reference."),
        }
        parts = []
        for key, (title, caption) in titles.items():
            plot_file = plots.get(key)
            if not plot_file or not Path(plot_file).exists():
                continue
            suffix = Path(plot_file).suffix.lower()
            mime = IMAGE_TYPES.get(suffix)
            if mime is None:
                inner = f'<p><a href="{self._relative(plot_file)}" target="_blank">Open {_e(title)} ({suffix[1:]})</a></p>'
            else:
                try:
                    data = base64.b64encode(Path(plot_file).read_bytes()).decode('utf-8')
                    inner = f'<img src="data:{mime};base64,{data}" alt="{_e(title)}">'
                except Exception as e:
                    logging.warning(f"Failed to embed plot {plot_file}: {e}")
                    continue
            parts.append(self._figure(inner, f"{_e(title)}. {caption}"))
        interactive = {k: v for k, v in plots.items() if str(v).lower().endswith('.html')}
        if interactive:
            labels = {'dashboard': 'Coverage profiles by set', 'metric_comparison': 'Metrics relative to the best set'}
            links = "".join(f'<a href="{self._relative(path)}" target="_blank"><strong>{_e(labels.get(key, key.replace("_", " ").title()))}'
                            f'</strong><span>{_e(Path(path).name)}</span></a>' for key, path in interactive.items())
            parts.append(f'<h3>Interactive figures</h3><div class="exports">{links}</div>')
        exports = ('<h3>Data files</h3><div class="exports">'
                   '<a href="comparison.json"><strong>comparison.json</strong><span>every result, with versions and arguments</span></a>'
                   '<a href="comparison_matrix.csv"><strong>comparison_matrix.csv</strong><span>the comparison table</span></a>'
                   '<a href="pairwise_comparisons.csv"><strong>pairwise_comparisons.csv</strong><span>pairwise differences</span></a>'
                   '<a href="gap_overlap_analysis.txt"><strong>gap_overlap_analysis.txt</strong><span>gap overlap</span></a></div>')
        return self._section("figures", "Figures and data files", "".join(parts) + exports)

    def _generate_methodology_section(self) -> str:
        weights = QualityScorer.DEFAULT_WEIGHTS
        thresholds = QualityScorer.DEFAULT_THRESHOLDS
        weight_items = "".join(f"<li><strong>{name.replace('_', ' ').title()}</strong>: {weight:.0%}</li>"
                               for name, weight in weights.items())
        legend = self.differential_analyzer.significance_legend() if self.differential_analyzer else "not run"
        mapper = self.metadata.get('mapper', getattr(self.analyzer, 'mapper', 'pblat'))
        body = f"""
<h4>Pipeline</h4>
<ol>
<li><strong>Mapping:</strong> oligos mapped to the reference with {_e(mapper)} and filtered on identity and aligned length</li>
<li><strong>Coverage:</strong> per-position depth from the aligned blocks; breadth, depth and uniformity statistics</li>
<li><strong>Gaps:</strong> runs below the minimum depth, at or above the minimum gap size, characterised by sequence features</li>
<li><strong>Quality score:</strong> weighted sum of five component scores in 0-1 (below)</li>
<li><strong>Statistical tests:</strong> non-parametric tests on observed window depths, per-reference metrics, bait identities and gap lengths</li>
<li><strong>Benchmarks:</strong> observed metrics against the design targets of QualityScorer.DEFAULT_BENCHMARKS</li>
<li><strong>Ranking:</strong> composite score in 0-1: quality score 40 %, coverage breadth 25 %, mapping efficiency 20 %, gap term 15 %</li>
</ol>
<h4>Quality score</h4>
<ul>{weight_items}</ul>
<p>Categories on the overall score (0-1): Excellent &gt;= {thresholds['excellent']}, Good &gt;= {thresholds['good']},
Fair &gt;= {thresholds['fair']}, otherwise Poor.</p>
<h4>Statistical tests</h4>
<ul>
<li><strong>Kolmogorov-Smirnov, Mann-Whitney U and Levene</strong> on mean depth per window (two sets)</li>
<li><strong>Wilcoxon signed-rank (two sets) or Friedman (more)</strong> on per-reference breadth, mean depth and gap count,
each reference being one paired observation</li>
<li><strong>Mann-Whitney U</strong> on best-hit identity per bait and on gap lengths (two sets)</li>
<li><strong>Effect sizes:</strong> rank-biserial correlation for Mann-Whitney, paired standardised mean difference for
Wilcoxon, Kendall's W for Friedman</li>
</ul>
<p>Significance labels: {_e(legend)}.</p>
"""
        return self._section("methods", "Methods", body)
