"""
Regression tests for the phase 1 fixes of the 2026-10 code audit
(docs/code_audit_2026-10.md): crash guards, the 0-1 quality scale in
compare, reference analysis in the outputs, consistent statistics in the
compare report, and HTML escaping.
"""

import json
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

from baitUtils.quality_scorer import QualityScorer, QualityScore, QualityCategory
from baitUtils.coverage_viz import CoverageVisualizer
from baitUtils.sequence_analysis import SequenceAnalyzer
from baitUtils.json_export import to_jsonable
from baitUtils.comparative_analyzer import ComparativeAnalyzer, OligoSetResult
from baitUtils.differential_analysis import DifferentialAnalyzer, StatisticalTest
from baitUtils.comparative_report_generator import ComparativeReportGenerator
from baitUtils.report_generator import InteractiveReportGenerator
from baitUtils.coverage_analysis import SequenceProcessor


def _quality(score: float) -> QualityScore:
    comps = {k: score for k in QualityScorer.DEFAULT_WEIGHTS}
    return QualityScore(score, QualityCategory.GOOD, comps, comps, {}, [])


def _oligo_set(name, breadth, efficiency, gaps, per_reference=None, coverage_arrays=None, score=0.8):
    return OligoSetResult(
        name=name, file_path=f"{name}.fa",
        coverage_stats={'coverage_breadth': breadth, 'mapping_efficiency': efficiency,
                        'mean_depth': 5.0, 'coverage_gini': 0.2,
                        'per_reference': per_reference or {}},
        gap_analysis={'total_gaps': len(gaps), 'gaps': gaps, 'max_gap_size': max([g['length'] for g in gaps], default=0)},
        quality_score=_quality(score), benchmark_results={},
        coverage_arrays=coverage_arrays or {}, per_oligo=None)


# --------------------------------------------------------------- crash guards

class TestCrashGuards:
    def test_reference_difficulty_with_expected_difficulty_two(self):
        """F-1: one extreme-GC, challenging reference used to divide by zero."""
        scorer = QualityScorer(
            coverage_stats={'coverage_breadth': 60.0, 'mean_depth': 5, 'coverage_cv': 0.4, 'mapping_efficiency': 90},
            gap_analysis={'gap_percentage': 5.0, 'max_gap_size': 500, 'total_gaps': 3},
            reference_analysis={'analysis_summary': {
                'total_sequences': 1, 'challenging_sequences': 1,
                'sequences_with_extreme_gc': 1, 'sequences_with_high_repeats': 0}})
        score = scorer.calculate_quality_score()
        assert 0.0 <= score.component_scores['reference_difficulty'] <= 1.0
        assert 0.0 <= score.overall_score <= 1.0

    def test_reference_difficulty_with_every_fraction_one(self):
        scorer = QualityScorer(
            coverage_stats={'coverage_breadth': 90.0},
            gap_analysis={},
            reference_analysis={'analysis_summary': {
                'total_sequences': 2, 'challenging_sequences': 2,
                'sequences_with_extreme_gc': 2, 'sequences_with_high_repeats': 2}})
        assert 0.0 <= scorer._score_reference_difficulty() <= 1.0

    def test_coverage_heatmap_with_no_mapped_hits(self, tmp_path):
        """F-2: all-zero per-reference statistics used to divide by zero."""
        stats = {'per_reference': {
            'chrA': {'coverage_breadth': 0.0, 'mean_depth': 0.0, 'max_depth': 0, 'gaps': 0},
            'chrB': {'coverage_breadth': 0.0, 'mean_depth': 0.0, 'max_depth': 0, 'gaps': 0}}}
        viz = CoverageVisualizer(coverage_stats=stats, gap_analysis={}, output_dir=tmp_path)
        viz.plot_coverage_heatmap()
        assert (tmp_path / "plots" / "coverage_heatmap.png").exists() or any(tmp_path.rglob("coverage_heatmap.*"))

    def test_empty_sequence_does_not_reach_viennarna(self):
        """F-3: fold_compound('&') segfaults; empty input must return NA."""
        analyzer = SequenceAnalyzer()
        assert analyzer.calculate_self_dimer_dg("") == 'NA'
        assert analyzer.calculate_hairpin_dg("") == 'NA'
        stats = analyzer.analyze_sequence("", "empty")
        assert stats['hairpin_dg'] == 'NA' and stats['self_dimer_dg'] == 'NA'


# ----------------------------------------------------------------- JSON / CLI

class TestSmallFixes:
    def test_infinite_values_become_null(self):
        """F-28: Infinity is not JSON."""
        out = to_jsonable({'cv': float('inf'), 'neg': float('-inf'), 'nan': float('nan'),
                           'np': np.float64('inf'), 'ok': 1.5})
        assert out == {'cv': None, 'neg': None, 'nan': None, 'np': None, 'ok': 1.5}
        json.dumps(out)

    def test_stats_missing_input_exits_nonzero(self, tmp_path, run_cli):
        """F-27"""
        with pytest.raises(SystemExit) as exc:
            run_cli(["stats", "-i", tmp_path / "missing.fa", "-o", tmp_path / "out"])
        assert exc.value.code == 1

    def test_uncovered_fasta_without_n_split(self, tmp_path):
        """F-24: the uncovered FASTA is written even without an N-split path."""
        ref = tmp_path / "ref.fa"
        ref.write_text(">chrA\n" + "ACGT" * 100 + "\n")
        out = tmp_path / "uncovered.fa"
        SequenceProcessor.export_uncovered_fasta({'chrA': [(10, 60)]}, ref, out, None)
        assert out.read_text().startswith(">chrA_uncovered_1")


# --------------------------------------------------------- compare: 0-1 scale

class TestCompareScale:
    def test_ranking_composite_is_in_0_1_and_quality_weight_counts(self):
        """F-4: the composite used to mix a 0-1 score with 0-10 terms."""
        analyzer = ComparativeAnalyzer(reference_file="ref.fa", output_dir=Path("."))
        high_q = _oligo_set("Q", breadth=80.0, efficiency=80.0, gaps=[], score=0.95)
        low_q = _oligo_set("L", breadth=80.0, efficiency=80.0, gaps=[], score=0.50)
        analyzer.oligo_sets = [low_q, high_q]
        ranking = dict(analyzer.generate_ranking())
        assert all(0.0 <= v <= 1.0 for v in ranking.values())
        # Same breadth, efficiency and gaps: the difference is 0.45 * 0.4
        assert ranking["Q"] - ranking["L"] == pytest.approx(0.45 * 0.4)

    def test_report_prints_score_on_its_own_scale(self, tmp_path):
        best = _oligo_set("Best", 90.0, 95.0, [], score=0.84)
        analyzer = MagicMock()
        analyzer.oligo_sets = [best, _oligo_set("Other", 80.0, 90.0, [], score=0.70)]
        analyzer.identify_best_performer.return_value = best
        analyzer.generate_ranking.return_value = [("Best", 0.86), ("Other", 0.74)]
        gen = ComparativeReportGenerator(analyzer, tmp_path)
        summary = gen._generate_executive_summary()
        assert "0.84" in summary and "/10" not in summary

    def test_evaluate_report_uses_scorer_category(self, tmp_path):
        """F-10: 0.843 is Good for the scorer; the report said Excellent."""
        gen = InteractiveReportGenerator(
            coverage_stats={}, gap_analysis={}, reference_analysis={}, output_dir=tmp_path,
            quality_scores={'overall_score': 0.843, 'category': 'Good'})
        gen._create_executive_summary()
        html = gen.html_components[-1]
        assert ">Good<" in html and "Excellent" not in html
        assert InteractiveReportGenerator._category_for(0.843) == "Good"
        assert InteractiveReportGenerator._category_for(0.85) == "Excellent"


# ------------------------------------------------------ compare: statistics

class TestDifferentialConsistency:
    def _sets(self):
        refs = {f"r{i}": None for i in range(8)}
        per_a = {r: {'coverage_breadth': 90 + i, 'mean_depth': 5 + i, 'gaps': 2} for i, r in enumerate(refs)}
        per_b = {r: {'coverage_breadth': 80 + i, 'mean_depth': 4 + i, 'gaps': 2} for i, r in enumerate(refs)}
        gaps_a = [{'chromosome': 'r0', 'start': 0, 'end': 150, 'length': 150}]
        gaps_b = [{'chromosome': 'r0', 'start': 0, 'end': 150, 'length': 150},
                  {'chromosome': 'r1', 'start': 0, 'end': 300, 'length': 300}]
        return _oligo_set("A", 90, 90, gaps_a, per_a), _oligo_set("B", 80, 90, gaps_b, per_b)

    def test_interpretation_follows_adjusted_p_value(self):
        """F-12: the text used the raw p, the stars the adjusted one."""
        a, b = self._sets()
        tests = DifferentialAnalyzer(significance_level=0.05, correction_method='bonferroni').compare_quality_metrics([a, b])
        for test in tests.values():
            if not test.applicable:
                continue
            significant = test.p_reported <= 0.05
            assert test.interpretation.startswith("Difference" if significant else "No difference")
            assert (test.significance_level != 'ns') == significant

    def test_star_labels_scale_with_alpha(self):
        """F-13"""
        strict = DifferentialAnalyzer(significance_level=0.01)
        assert strict._get_significance_level(0.03) == 'ns'
        assert strict._get_significance_level(0.009) == '*'
        assert strict._get_significance_level(0.0001) == '***'
        assert "p<=0.01" in strict.significance_legend()
        default = DifferentialAnalyzer()
        assert default._get_significance_level(0.03) == '*'
        assert default._get_significance_level(0.0005) == '***'

    def test_gap_count_uses_gap_analysis_definition(self):
        """F-14: per-reference gap counts come from the gap analysis, not from every run below min depth."""
        a, b = self._sets()
        analyzer = DifferentialAnalyzer()
        counts_a = analyzer._gaps_per_reference(a)
        counts_b = analyzer._gaps_per_reference(b)
        assert counts_a == {'r0': 1} and counts_b == {'r0': 1, 'r1': 1}
        test = analyzer.compare_quality_metrics([a, b])['gap_count']
        # Only r1 differs (0 vs 1), so the Wilcoxon test runs on that difference
        assert test.applicable and test.n == 8

    def test_touching_gaps_do_not_overlap(self):
        """F-15: half-open intervals that touch share no base."""
        analyzer = ComparativeAnalyzer(reference_file="ref.fa", output_dir=Path("."))
        g1 = [{'chromosome': 'c', 'start': 0, 'end': 100, 'length': 100}]
        g2 = [{'chromosome': 'c', 'start': 100, 'end': 200, 'length': 100}]
        assert analyzer._find_overlapping_gaps(g1, g2) == []
        assert analyzer._find_unique_gaps('x', {'x': g1, 'y': g2}) == g1

    def test_interactive_links_are_relative_to_the_report(self, tmp_path):
        """F-16"""
        analyzer = MagicMock()
        gen = ComparativeReportGenerator(analyzer, tmp_path)
        plot = tmp_path / "interactive_comparative_plots" / "comparison_dashboard.html"
        plot.parent.mkdir(exist_ok=True)
        plot.write_text("<html></html>")
        section = gen._generate_visualizations_section({'dashboard': str(plot)})
        assert 'href="interactive_comparative_plots/comparison_dashboard.html"' in section


# --------------------------------------------------------------- HTML escaping

class TestHtmlEscaping:
    def test_compare_report_escapes_set_names_and_test_text(self, tmp_path):
        """F-11"""
        bad = "A<b>x</b>"
        best = _oligo_set(bad, 90.0, 95.0, [], score=0.84)
        analyzer = MagicMock()
        analyzer.oligo_sets = [best, _oligo_set("B", 80.0, 90.0, [], score=0.70)]
        analyzer.identify_best_performer.return_value = best
        analyzer.generate_ranking.return_value = [(bad, 0.86), ("B", 0.74)]
        gen = ComparativeReportGenerator(analyzer, tmp_path)
        summary = gen._generate_executive_summary()
        assert "<b>x</b>" not in summary and "A&lt;b&gt;x&lt;/b&gt;" in summary
        test = StatisticalTest("T<i>", 1.0, 0.5, 0.1, 'ns', "<script>", None, True, 3)
        rendered = ComparativeReportGenerator._render_test("Title", test, 0.05)
        assert "<script>" not in rendered and "&lt;script&gt;" in rendered

    def test_evaluate_report_escapes_reference_ids(self, tmp_path):
        gen = InteractiveReportGenerator(
            coverage_stats={}, quality_scores={}, reference_analysis={}, output_dir=tmp_path,
            gap_analysis={'largest_gaps': [{'chromosome': '<img src=x onerror=alert(1)>', 'start': 1, 'end': 2, 'length': 1}]})
        table = gen._create_largest_gaps_table()
        assert "<img" not in table and "&lt;img" in table


# ------------------------------------------------------ evaluate: JSON content

class TestEvaluateJson:
    def test_reference_analysis_is_written(self, dataset, tmp_path, run_cli, fake_pblat):
        """F-5 and F-22: the correlations and challenging windows reach the outputs."""
        out = tmp_path / "eval"
        run_cli(["evaluate", "-i", dataset["baits"], "-r", dataset["reference"], "-o", out])
        data = json.loads((out / "evaluation.json").read_text())
        assert data["reference_summary"] is not None
        assert data["reference_summary"]["total_sequences"] == 2
        ref = data["reference_analysis"]
        assert set(ref) == {"sequence_features", "challenging_regions", "coverage_correlations", "challenging_windows"}
        assert ref["coverage_correlations"]["applicable"] is True
        assert "gc_content" in ref["coverage_correlations"]["features"]
        assert ref["coverage_correlations"]["breadth_by_gc"]
        html = (out / "coverage_evaluation_report.html").read_text()
        assert "Feature-Coverage Correlations" in html and "Breadth by GC decile" in html
        assert "Quality Trends" not in (out / "interactive_plots" / "quality_assessment.html").read_text()
        # gap features are produced (used by the interactive gap panels) but kept out of the JSON
        assert "gap_features" not in data["gap_analysis"]


# --------------------------------------------------- compare: plot settings

class TestComparePlotSettings:
    def test_plot_format_and_dpi_are_kept_by_the_report(self, dataset, tmp_path, run_cli, fake_pblat):
        """F-6: the report used to regenerate every plot with default settings."""
        out = tmp_path / "cmp"
        run_cli(["compare", "-r", dataset["reference"], "-o", out,
                 "--sets", f"A:{dataset['baits']}", f"B:{dataset['baits']}",
                 "--plot-format", "svg", "--plot-dpi", "72", "--threads", "1"])
        static = {p.suffix for p in (out / "comparative_plots").iterdir()}
        assert static == {".svg"}
        html = (out / "comparative_analysis_report.html").read_text()
        assert "data:image/svg+xml;base64" in html and "data:image/png" not in html
        assert 'href="interactive_comparative_plots/comparison_dashboard.html"' in html
        assert "/10" not in html
        assert "Significance labels: *** p<=0.001, ** p<=0.01, * p<=0.05, ns p>0.05" in html
        data = json.loads((out / "comparison.json").read_text())
        assert data["reference_analysis"]["coverage_correlations"]["applicable"] is True
        assert all(0.0 <= score <= 1.0 for _, score in data["ranking"])
