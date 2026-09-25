"""
Value-level tests for sequence_features and ReferenceAnalyzer.
"""

import numpy as np
import pytest

from baitUtils import sequence_features as sf
from baitUtils.coverage_stats import CoverageAnalyzer
from baitUtils.reference_analyzer import ReferenceAnalyzer


class TestSequenceFeatures:
    def test_encode_and_counts(self):
        codes = sf.encode("ACGTNacgtX")
        assert list(codes) == [0, 1, 2, 3, 4, 0, 1, 2, 3, 4]
        assert list(sf.base_counts(codes)) == [2, 2, 2, 2, 2]
        assert list(sf.masked_mask("ACgtN")) == [False, False, True, True, False]

    def test_gc_and_n_content(self):
        codes = sf.encode("GGCCAANN")
        assert sf.gc_content(codes) == pytest.approx(4 / 6 * 100)
        assert sf.n_content(codes) == pytest.approx(25.0)

    def test_entropy(self):
        assert sf.shannon_entropy(sf.encode("ACGT" * 10)) == pytest.approx(2.0)
        assert sf.shannon_entropy(sf.encode("AAAA")) == 0.0
        assert sf.shannon_entropy(sf.encode("NNNN")) == 0.0

    def test_homopolymers(self):
        codes = sf.encode("ACGTTTTTTGCAAAAANNNNNNNN")
        assert sf.max_homopolymer(codes) == 6
        mask = sf.homopolymer_mask(codes, min_run=5)
        assert mask.sum() == 6 + 5   # the N run is not counted

    def test_repeat_content_from_duplicated_kmers(self):
        unique = "ACGGTCAAGTCCGATTGCAGGCTATCGATCCGTAAGCTTAGCATGCCGTAGGATCCAGTTGACCGTACGGATCGAAGTCTAGC"
        assert sf.repeat_content(sf.encode(unique), k=12) == 0.0
        repeated = "ACGGTCAAGTCCGATTGCAG" * 4
        assert sf.repeat_content(sf.encode(repeated), k=12) > 90.0

    def test_linguistic_complexity_bounds(self):
        rng = np.random.default_rng(0)
        random_seq = "".join(rng.choice(list("ACGT"), 2000))
        assert sf.linguistic_complexity(sf.encode(random_seq)) > 0.8
        assert sf.linguistic_complexity(sf.encode("ATATATATATATATATATAT")) < 0.2

    def test_window_starts_cover_tail(self):
        assert list(sf.window_starts(25, 10, 10)) == [0, 10, 20]
        assert list(sf.window_starts(30, 10, 10)) == [0, 10, 20]
        assert list(sf.window_starts(5, 10, 10)) == [0]
        assert list(sf.window_starts(0, 10, 10)) == []

    def test_window_features_match_direct_computation(self):
        seq = "GGGGGGGGGGAAAAAAAAAACATGCATGCANNNNNNNNNN"
        codes = sf.encode(seq)
        w = sf.window_features(codes, 10, 10, homopolymer_min_run=5, repeat_k=4)
        assert list(w["start"]) == [0, 10, 20, 30]
        assert w["gc_content"][0] == pytest.approx(100.0)
        assert w["gc_content"][1] == pytest.approx(0.0)
        assert w["gc_content"][2] == pytest.approx(50.0)
        assert w["n_content"][3] == pytest.approx(100.0)
        assert w["homopolymer_content"][0] == pytest.approx(100.0)
        assert w["homopolymer_content"][2] == pytest.approx(0.0)
        assert w["entropy"][2] > 1.9

    def test_window_coverage(self):
        cov = np.array([0, 0, 1, 2, 3, 0, 0, 0, 0, 5])
        w = sf.window_coverage(cov, 5, 5, min_coverage=1)
        assert list(w["breadth"]) == [60.0, 20.0]
        assert list(w["mean_depth"]) == [1.2, 1.0]


@pytest.fixture(scope="module")
def coverage(dataset):
    analyzer = CoverageAnalyzer(psl_file=dataset["psl"], reference_file=dataset["reference"],
                                min_coverage=1.0, oligos_file=dataset["baits"])
    analyzer.analyze()
    return analyzer


class TestReferenceAnalyzer:
    def test_sequence_features_and_summary(self, dataset, coverage):
        results = ReferenceAnalyzer(dataset["reference"], coverage.stats, window_size=500,
                                    coverage_arrays=coverage.coverage_arrays).analyze()
        feats = results["sequence_features"]
        assert set(feats) == {"chrA", "chrB"}
        assert feats["chrA"]["length"] == 3000
        assert 40 < feats["chrA"]["gc_content"] < 60
        assert feats["chrA"]["windows"] == 6
        assert feats["chrB"]["windows"] == 4
        summary = results["analysis_summary"]
        assert summary["total_sequences"] == 2
        assert summary["total_length"] == 5000
        assert summary["correlations_applicable"] is True
        for ref in ("chrA", "chrB"):
            assert results["challenging_regions"][ref]["difficulty_level"] == "Easy"

    def test_window_features_carry_coverage(self, dataset, coverage):
        results = ReferenceAnalyzer(dataset["reference"], coverage.stats, window_size=500,
                                    coverage_arrays=coverage.coverage_arrays).analyze()
        windows = results["window_features"]["chrA"]
        assert [w["start"] for w in windows] == [0, 500, 1000, 1500, 2000, 2500]
        # The 300 bp hole at 1080-1380 lowers breadth in windows 1000-1500
        assert windows[2]["breadth"] == pytest.approx((500 - 300) / 500 * 100)
        assert windows[0]["breadth"] == pytest.approx(100.0)
        assert windows[5]["breadth"] == pytest.approx((500 - 60) / 500 * 100)

    def test_correlations_are_real_spearman_values(self, dataset, coverage):
        results = ReferenceAnalyzer(dataset["reference"], coverage.stats, window_size=200,
                                    coverage_arrays=coverage.coverage_arrays).analyze()
        corr = results["coverage_correlations"]
        assert corr["applicable"] is True
        assert corr["windows"] == 15 + 10
        gc = corr["features"]["gc_content"]
        assert -1.0 <= gc["rho_breadth"] <= 1.0
        assert 0.0 <= gc["p_breadth"] <= 1.0
        assert sum(b["windows"] for b in corr["breadth_by_gc"]) == 25

    def test_without_arrays_correlations_not_applicable(self, dataset, coverage):
        results = ReferenceAnalyzer(dataset["reference"], coverage.stats).analyze()
        assert results["coverage_correlations"]["applicable"] is False
        assert results["challenging_windows"] == {}
        assert results["analysis_summary"]["correlations_applicable"] is False

    def test_challenging_reference_and_windows(self, tmp_path):
        ref = tmp_path / "hard.fa"
        ref.write_text(">hard\n" + "A" * 400 + "ACGT" * 100 + "N" * 200 + "\n")
        cov = {"hard": np.concatenate([np.zeros(400), np.ones(400), np.zeros(200)]).astype(np.int32)}
        results = ReferenceAnalyzer(ref, {}, window_size=200, coverage_arrays=cov).analyze()
        region = results["challenging_regions"]["hard"]
        assert region["challenging_score"] >= 3
        assert any("homopolymer" in r.lower() for r in region["reasons"])
        assert any("N content" in r for r in region["reasons"])
        windows = results["challenging_windows"]["hard"]
        assert [w["start"] for w in windows] == [0, 200, 800]
        assert "homopolymer-rich" in windows[0]["reasons"]
        assert "N-rich" in windows[2]["reasons"]
