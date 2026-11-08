"""
Regression tests for phase 3 of the second audit (docs/code_audit_2026-11.md):
the compare statistics, ranking and report. Each test failed before its fix.
"""

import csv
import json
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest

from baitUtils.comparative_analyzer import ComparativeAnalyzer, OligoSetResult
from baitUtils.differential_analysis import DifferentialAnalyzer


def _set(name, arrays=None, per_oligo=None, gaps=None, **stats):
    coverage_stats = {'per_reference': {}, **stats}
    return OligoSetResult(name, f"{name}.fa", coverage_stats, {'gaps': gaps or []}, MagicMock(),
                          coverage_arrays=arrays, per_oligo=per_oligo)


# ------------------------------------------------------------------ N-41

class TestIdentityTest:
    def test_baits_excluded_by_max_hits_are_not_compared(self):
        """N-41: with --max-hits 0 nothing was mapped, yet the identity test reported p=1e-8 '***'."""
        rng = np.random.default_rng(0)
        table_a = pd.DataFrame({'best_identity': rng.uniform(95, 100, 40), 'multi_mapping_excluded': True})
        table_b = pd.DataFrame({'best_identity': rng.uniform(85, 90, 40), 'multi_mapping_excluded': True})
        result = DifferentialAnalyzer().compare_oligo_identity(_set("A", per_oligo=table_a),
                                                               _set("B", per_oligo=table_b))
        assert not result.applicable

    def test_only_mapped_baits_count(self):
        rng = np.random.default_rng(0)
        table_a = pd.DataFrame({'best_identity': rng.uniform(95, 100, 40),
                                'multi_mapping_excluded': [i < 10 for i in range(40)]})
        table_b = pd.DataFrame({'best_identity': rng.uniform(85, 90, 40), 'multi_mapping_excluded': False})
        result = DifferentialAnalyzer().compare_oligo_identity(_set("A", per_oligo=table_a),
                                                               _set("B", per_oligo=table_b))
        assert result.applicable and result.n == 30 + 40


# ------------------------------------------------------------ N-42, N-45

class TestWindowDepthTests:
    def test_paired_test_sees_a_shift_in_every_window(self):
        """N-42: B = A + 1x in every window; KS, Mann-Whitney and Levene all said 'no difference'."""
        rng = np.random.default_rng(0)
        base = np.repeat(rng.integers(5, 60, size=20), 1000)
        a, b = _set("A", {'chr1': base}), _set("B", {'chr1': base + 1})
        result = DifferentialAnalyzer().compare_coverage_distributions(b, a)
        assert result.wilcoxon_test.applicable
        assert result.wilcoxon_test.p_adjusted < 0.001
        assert result.wilcoxon_test.effect_size == pytest.approx(1.0)      # B deeper in every window
        assert result.wilcoxon_test.n == 20
        assert result.ks_test.significance_level == "ns"                   # shape unchanged

    def test_windows_are_paired_by_reference_name(self):
        """Pairs must not depend on the order in which each set stores its references."""
        rng = np.random.default_rng(1)
        x, y = rng.integers(0, 20, 3000), rng.integers(40, 60, 3000)
        a = _set("A", {'chrX': x, 'chrY': y})
        b = _set("B", {'chrY': y, 'chrX': x})      # same data, other order
        result = DifferentialAnalyzer().compare_coverage_distributions(a, b)
        assert not result.wilcoxon_test.applicable    # identical depth in every window

    def test_all_tied_mann_whitney_is_not_applicable(self):
        """N-45: all values tied gave 'applicable', p null and adjusted p 1.0."""
        gaps = [{'chromosome': 'c', 'start': i * 1000, 'end': i * 1000 + 300, 'length': 300} for i in range(3)]
        result = DifferentialAnalyzer().analyze_gap_patterns(_set("A", gaps=gaps), _set("B", gaps=gaps))
        assert not result['gap_size_distribution'].applicable

    def test_wilcoxon_n_counts_references_that_differ(self):
        per_ref_a = {f"r{i}": {'coverage_breadth': 90.0, 'mean_depth': 5.0, 'gaps': 1} for i in range(6)}
        per_ref_b = {f"r{i}": {'coverage_breadth': 90.0 + (i < 3), 'mean_depth': 5.0, 'gaps': 1} for i in range(6)}
        a, b = _set("A"), _set("B")
        a.coverage_stats['per_reference'], b.coverage_stats['per_reference'] = per_ref_a, per_ref_b
        breadth = DifferentialAnalyzer().compare_quality_metrics([a, b])['coverage_breadth']
        assert breadth.n == 3 and "3 of 6 references differ" in breadth.interpretation

    def test_half_integer_u_is_printed_in_full(self):
        gaps_a = [{'chromosome': 'c', 'start': 0, 'end': n, 'length': n} for n in (100, 200, 300)]
        gaps_b = [{'chromosome': 'c', 'start': 0, 'end': n, 'length': n} for n in (200, 250)]
        test = DifferentialAnalyzer().analyze_gap_patterns(_set("A", gaps=gaps_a), _set("B", gaps=gaps_b))
        test = test['gap_size_distribution']
        assert f"U={test.statistic:g}" in test.interpretation


# ------------------------------------------------------------ N-43, N-44

def _fixture_compare(run_cli, dataset, out, *sets):
    run_cli(["compare", "-r", dataset["reference"], "-o", out, "--threads", "1",
             "--sets", *[f"{name}:{dataset[key]}" for name, key in sets]])
    return json.loads((out / "comparison.json").read_text())


class TestComparisonJson:
    def test_reference_coverage_is_per_set_and_order_independent(self, dataset, tmp_path, run_cli, fake_pblat):
        """N-43: the correlations came from the first set only and changed with the order of --sets."""
        first = _fixture_compare(run_cli, dataset, tmp_path / "a", ("Full", "baits"), ("Sub", "baits_subset"))
        second = _fixture_compare(run_cli, dataset, tmp_path / "b", ("Sub", "baits_subset"), ("Full", "baits"))

        def per_set(data):
            return {s['name']: s['reference_coverage']['coverage_correlations']['breadth_by_gc'] for s in data['sets']}
        assert per_set(first) == per_set(second)
        assert per_set(first)['Full'] != per_set(first)['Sub']
        assert set(first['reference_analysis']) == {'sequence_features', 'challenging_regions'}

    def test_gap_term_follows_the_share_in_gaps(self):
        """N-44: a set covering nothing (one gap per reference) got almost the full 15 % gap weight."""
        analyzer = ComparativeAnalyzer(reference_file="ref.fa", output_dir=Path("."))

        def result(name, breadth, gap_pct, gaps, quality):
            r = _set(name, coverage_breadth=breadth, mapping_efficiency=breadth)
            r.gap_analysis = {'gap_percentage': gap_pct, 'total_gaps': gaps}
            r.quality_score = MagicMock(overall_score=quality)
            return r
        analyzer.oligo_sets = [result("none", 0.0, 100.0, 2, 0.05), result("fragmented", 99.0, 1.0, 50, 0.8)]
        scores = dict(analyzer.generate_ranking())
        assert scores["none"] == pytest.approx(0.05 * 0.4)                  # no gap credit
        assert scores["fragmented"] > 0.9


# ------------------------------------------------------- N-46 to N-50

class TestCompareDetails:
    def test_uniformity_diff_is_positive_for_the_more_uniform_first_set(self, dataset, tmp_path, run_cli, fake_pblat):
        """N-46: Uniformity_Diff was the Gini difference, the opposite sign of its name."""
        data = _fixture_compare(run_cli, dataset, tmp_path / "c", ("Full", "baits"), ("Sub", "baits_subset"))
        gini = {s['name']: s['coverage_stats']['coverage_gini'] for s in data['sets']}
        rows = list(csv.DictReader(open(tmp_path / "c" / "pairwise_comparisons.csv")))
        assert float(rows[0]['Uniformity_Diff']) == pytest.approx(gini['Sub'] - gini['Full'])

    def test_recommendations_name_the_ranking_winner_once(self):
        """N-47: the best set by quality and ranking[1] could be the same set, and the winner was not named."""
        from baitUtils.comparative_report_generator import ComparativeReportGenerator
        from baitUtils.quality_scorer import QualityCategory

        def result(name, quality):
            r = _set(name, coverage_breadth=90.0, mapping_efficiency=95.0)
            r.gap_analysis = {'total_gaps': 1}
            r.quality_score = MagicMock(overall_score=quality, category=QualityCategory.GOOD)
            return r
        deep, plain = result("Deep", 0.844), result("Plain", 0.837)
        analyzer = MagicMock(oligo_sets=[deep, plain])
        analyzer.identify_best_performer.return_value = deep
        analyzer.generate_ranking.return_value = [("Plain", 0.893), ("Deep", 0.858)]
        html = ComparativeReportGenerator(analyzer, Path("."))._generate_recommendations()
        assert "<strong>Plain</strong> ranks first" in html
        assert html.count("<strong>Deep</strong>") == 2       # ranks second, and has the best quality score

    def test_overlapping_gaps_of_different_sets_are_one_region(self):
        """N-48: the same hole with different edges counted as two 'distinct' regions."""
        gaps = [{'chromosome': 'chrA', 'start': 1080, 'end': 1380}, {'chromosome': 'chrA', 'start': 1080, 'end': 1860},
                {'chromosome': 'chrB', 'start': 0, 'end': 100}, {'chromosome': 'chrB', 'start': 100, 'end': 200}]
        assert ComparativeAnalyzer._count_merged_regions(gaps) == 3     # touching intervals share no base

    def test_effect_panel_shows_only_signed_effects(self):
        """N-49: KS D, the Levene ratio and Kendall's W were drawn as 'first set larger'."""
        from baitUtils.comparative_visualizations import signed_effects

        def t(name, effect):
            return MagicMock(test_name=name, effect_size=effect)
        tests = [("ks", t("Kolmogorov-Smirnov", 0.4)), ("lev", t("Levene", 0.18)), ("fr", t("Friedman", 0.5)),
                 ("wx", t("Wilcoxon signed-rank (paired windows)", -0.6)), ("mw", t("Mann-Whitney U (gap sizes)", 0.2))]
        assert [label for label, _ in signed_effects(tests)] == ["wx", "mw"]

    def test_empty_set_name_is_an_error(self, dataset, tmp_path, run_cli, fake_pblat):
        """N-50: ':baits.fa' was accepted and shown as a blank set."""
        with pytest.raises(SystemExit) as excinfo:
            run_cli(["compare", "-r", dataset["reference"], "-o", tmp_path / "o",
                     "--sets", f":{dataset['baits']}", f" :{dataset['baits_subset']}"])
        assert excinfo.value.code == 1
