"""
Regression tests for phase 2 of the second audit (docs/code_audit_2026-11.md):
statistics over A, C, G and T only, undefined without them, and the gap,
score and coverage edge cases. Each test failed before its fix.
"""

import json
import math

import numpy as np
import pytest

from baitUtils import sequence_features as sf
from baitUtils.gap_analysis import GapAnalyzer
from baitUtils.reference_analyzer import ReferenceAnalyzer
from baitUtils.quality_scorer import QualityScorer


def _random(rng, n):
    return "".join(rng.choice(list("ACGT"), n))


# ------------------------------------------------------------------ N-28

class TestNWindows:
    def _analyze(self):
        """Four references: 3000 random bp at depth 5, then 1000 N (uncoverable)."""
        rng = np.random.default_rng(1)
        seqs, cov = {}, {}
        for i in range(4):
            seqs[f"r{i}"] = _random(rng, 3000) + "N" * 1000
            cov[f"r{i}"] = np.r_[np.full(3000, 5), np.zeros(1000)].astype(int)
        return ReferenceAnalyzer("x.fa", {}, window_size=1000, coverage_arrays=cov,
                                 reference_sequences=seqs).analyze()

    def test_n_windows_have_no_composition(self):
        features = sf.window_features(sf.encode("ACGT" * 250 + "N" * 1000), 1000, 1000)
        assert features["gc_content"][0] == pytest.approx(50.0)
        for key in ("gc_content", "entropy", "homopolymer_content", "repeat_density"):
            assert np.isnan(features[key][1]), key
        assert features["n_content"][1] == 100.0
        assert math.isnan(sf.gc_content(sf.encode("NNNN")))

    def test_n_windows_create_no_correlation_decile_or_flag(self):
        """N-28: N windows gave rho 0.76 between GC and breadth on data where coverage ignores GC."""
        result = self._analyze()
        corr = result["coverage_correlations"]
        assert corr["windows"] == 12                       # the 4 N-only windows are not assessable
        assert math.isnan(corr["features"]["gc_content"]["rho_breadth"])   # breadth is constant
        assert all(b["gc_min"] >= 10 for b in corr["breadth_by_gc"])       # no "0-10 %" bin
        assert sum(len(w) for w in result["challenging_windows"].values()) == 0

    def test_window_breadth_is_over_non_n_bases(self):
        coverage = np.r_[np.ones(500), np.zeros(500)]
        n_mask = np.r_[np.zeros(500, bool), np.ones(400, bool), np.zeros(100, bool)]
        cov = sf.window_coverage(coverage, 1000, 1000, 1.0, n_mask=n_mask)
        assert cov["breadth"][0] == pytest.approx(500 / 600 * 100)    # 100 real uncovered bases


# ------------------------------------------------------------------ N-29

class TestGapFeaturesOverAcgt:
    def test_n_gap_is_neither_gc_extreme_nor_low_complexity(self):
        """N-29: a gap over 500 N had GC 0, a 500 bp homopolymer and 'extreme composition'."""
        features = GapAnalyzer({}, "x")._calculate_sequence_features("N" * 500)
        assert math.isnan(features["gc_content"])
        assert features["max_homopolymer"] == 0
        assert not features["extreme_composition"] and not features["at_rich"]
        assert features["n_content"] == 100.0

    def test_suggestions_for_an_n_gap_name_only_n(self):
        rng = np.random.default_rng(0)
        seq = _random(rng, 1000) + "N" * 500 + _random(rng, 1000)
        cov = np.ones(len(seq), int)
        cov[900:1600] = 0                       # real bases on both sides of the N run uncovered
        result = GapAnalyzer({'reference_length': len(seq), 'assessable_bases': len(seq) - 500},
                             "x", min_gap_size=50, coverage_arrays={'r': cov},
                             reference_sequences={'r': seq}).analyze()
        text = " ".join(result["suggestions"])
        assert "extreme GC" not in text and "low-complexity" not in text


# ------------------------------------------------------- N-30, N-31, N-33

class TestGapsOverRealSequence:
    def test_n_run_is_not_a_gap_and_splits_gaps(self):
        """N-30: a 5000 bp N run made a 5040 bp gap and pulled the gap score to 0.47."""
        rng = np.random.default_rng(0)
        seq = _random(rng, 8000) + "N" * 5000 + _random(rng, 7000)
        cov = np.ones(len(seq), int)
        cov[7980:13060] = 0                     # 20 + 5000 N + 60 uncovered
        result = GapAnalyzer({'reference_length': len(seq), 'assessable_bases': 15000},
                             "x", min_gap_size=50, coverage_arrays={'r': cov},
                             reference_sequences={'r': seq}).analyze()
        assert [(g['start'], g['end']) for g in result['gaps']] == [(13000, 13060)]
        assert result['gap_percentage'] == pytest.approx(60 / 15000 * 100)

    def test_score_uses_breadth_over_non_n_bases(self):
        stats = {'coverage_breadth': 74.7, 'coverage_breadth_non_n': 99.6, 'mean_depth': 20.0,
                 'coverage_gini': 0.1, 'mapping_efficiency': 100.0}
        gaps = {'gap_percentage': 0.0, 'max_gap_size': 0, 'total_gaps': 0}
        score = QualityScorer(stats, gaps).calculate_quality_score()
        assert score.component_scores['coverage_breadth'] == 1.0

    def test_all_n_contig_is_not_extreme_gc(self):
        rng = np.random.default_rng(0)
        seqs = {'good': _random(rng, 3000), 'allN': "N" * 2000}
        cov = {'good': np.ones(3000, int), 'allN': np.zeros(2000, int)}
        result = ReferenceAnalyzer("x.fa", {}, window_size=1000, coverage_arrays=cov,
                                   reference_sequences=seqs).analyze()
        reasons = result["challenging_regions"]["allN"]["reasons"]
        assert not any("GC" in r for r in reasons)
        assert any("N content" in r for r in reasons)

    def test_repeat_content_does_not_grow_with_gap_length(self):
        """N-31: 6-mer duplication was 15 % at 500 bp and 70 % at 5 kb on random sequence."""
        rng = np.random.default_rng(0)
        analyzer = GapAnalyzer({}, "x")
        assert analyzer._calculate_repeat_content(_random(rng, 5000)) < 1.0
        assert analyzer._calculate_repeat_content("ACG" * 400) > 90.0

    def test_negative_extension_larger_than_the_gap_does_not_crash(self):
        """N-33: an empty slice returned a feature dict without 'extreme_composition'."""
        rng = np.random.default_rng(0)
        seq = _random(rng, 3000)
        cov = np.ones(3000, int)
        cov[1000:1300] = 0
        result = GapAnalyzer({'reference_length': 3000}, "x", min_gap_size=100, extend_bp=-200,
                             coverage_arrays={'r': cov}, reference_sequences={'r': seq}).analyze()
        assert result['total_gaps'] == 1


# ------------------------------------------------------------ N-26, N-27

class TestNothingCovered:
    def test_uniformity_is_undefined_and_not_recommended_on(self, dataset, tmp_path):
        """N-27: Gini 1.0 and CV inf gave 'Poor' uniformity and 'CV=inf' advice with nothing covered."""
        from baitUtils.coverage_stats import CoverageAnalyzer
        from baitUtils.benchmark import BenchmarkAnalyzer
        empty = tmp_path / "empty.psl"
        empty.write_text("")
        stats = CoverageAnalyzer(empty, dataset["reference"]).analyze()
        assert math.isnan(stats['coverage_gini']) and math.isnan(stats['coverage_cv'])
        result = BenchmarkAnalyzer(stats, {'gap_percentage': 100.0}).benchmark_depth_uniformity()
        assert result.category == "Not applicable"
        score = QualityScorer(stats, {'gap_percentage': 100.0, 'max_gap_size': 3000, 'total_gaps': 2})
        recommendations = score.calculate_quality_score().recommendations
        assert not any("CV=" in r for r in recommendations)

    def test_min_coverage_zero_is_rejected(self, run_cli, capsys):
        """N-26: --min-coverage 0 counted every base as covered."""
        for command in ("evaluate", "compare", "check"):
            with pytest.raises(SystemExit) as excinfo:
                run_cli([command, "--min-coverage", "0", "--help"])
            assert excinfo.value.code == 2, command
        assert "must be above 0" in capsys.readouterr().err


# ------------------------------------------------------------ N-34, N-35

class TestGapReporting:
    def test_size_classes_agree(self):
        """N-34: [lo, hi) in one place and (lo, hi] in the other put 500 and 1000 bp in different classes."""
        rng = np.random.default_rng(0)
        seq = _random(rng, 6000)
        cov = np.ones(6000, int)
        cov[100:600] = 0                        # 500 bp
        cov[1000:2000] = 0                      # 1000 bp
        analyzer = GapAnalyzer({'reference_length': 6000}, "x", min_gap_size=100,
                               coverage_arrays={'r': cov}, reference_sequences={'r': seq})
        result = analyzer.analyze()
        assert result['size_distribution']['500-1000'] == 1
        assert result['size_distribution']['1000-5000'] == 1
        assert set(result['size_stratified_analysis']) == {'500-1000', '1000-5000'}

    def test_first_class_starts_at_the_minimum_gap_size(self):
        labels = [label for label, _, _ in GapAnalyzer({}, "x", min_gap_size=1000)._size_classes()]
        assert labels[0] == "1000-5000"

    def test_text_report_uses_one_based_coordinates(self, dataset, tmp_path, run_cli, fake_pblat):
        """N-35: 'chrA:1080-1380' was a 0-based start in samtools notation."""
        out = tmp_path / "eval"
        run_cli(["evaluate", "-i", dataset["baits"], "-r", dataset["reference"], "-o", out,
                 "--no-interactive-plots", "--no-html-report"])
        text = (out / "gap_analysis.txt").read_text()
        assert "chrA:1081-1380" in text
        data = json.loads((out / "evaluation.json").read_text())
        assert {'start': 1080, 'end': 1380} in [{'start': g['start'], 'end': g['end']}
                                                for g in data['gap_analysis']['gaps']]


# ----------------------------------------------------------- N-11, N-12, F-52

class TestStatsOverAcgt:
    def test_tm_entropy_and_homopolymers_ignore_n(self):
        from baitUtils.sequence_analysis import SequenceAnalyzer
        rng = np.random.default_rng(3)
        bait = _random(rng, 120)
        analyzer = SequenceAnalyzer()
        with_n = bait[:60] + "N" * 10 + bait[60:]
        assert analyzer.calculate_melting_temperature(with_n) == 'NA'           # N-11
        result = analyzer.analyze_sequence(with_n, "x")
        assert result['max_homopolymer_length'] < 10                            # N-12
        assert analyzer.calculate_entropy("ACGTNNNN") == 2.0                    # F-52
        assert analyzer.calculate_entropy("NNNN") == 'NA'
