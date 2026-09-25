"""
Coverage, gap and mapping-efficiency values computed on the hand-derived
fixture dataset (see conftest.py for the expected numbers).
"""

import numpy as np
import pytest

from baitUtils.coverage_stats import CoverageAnalyzer, find_gap_intervals
from baitUtils.gap_analysis import GapAnalyzer


class TestFindGapIntervals:
    def test_empty_array(self):
        assert find_gap_intervals(np.array([]), 1) == []

    def test_no_gaps(self):
        assert find_gap_intervals(np.ones(10), 1) == []

    def test_gaps_at_both_ends_and_middle(self):
        arr = np.array([0, 0, 1, 1, 0, 1, 0, 0])
        assert find_gap_intervals(arr, 1) == [(0, 2), (4, 5), (6, 8)]

    def test_threshold_is_inclusive_of_min_coverage(self):
        arr = np.array([2, 1, 0, 3])
        assert find_gap_intervals(arr, 2) == [(1, 3)]


@pytest.fixture(scope="module")
def analyzer(dataset):
    analyzer = CoverageAnalyzer(
        psl_file=dataset["psl"],
        reference_file=dataset["reference"],
        min_coverage=1.0,
        min_identity=90.0,
        min_length=100,
        oligos_file=dataset["baits"],
    )
    analyzer.analyze()
    return analyzer


class TestCoverageAnalyzerValues:
    def test_mapping_efficiency_uses_input_fasta(self, analyzer, dataset):
        exp = dataset["expected"]
        assert analyzer.stats["total_oligos"] == exp["n_baits"]
        assert analyzer.stats["mapped_oligos"] == exp["n_mapped"]
        assert analyzer.stats["mapping_efficiency"] == pytest.approx(
            100.0 * exp["n_mapped"] / exp["n_baits"])

    def test_covered_bases_per_reference(self, analyzer, dataset):
        per_ref = analyzer.stats["per_reference"]
        for chrom, covered in dataset["expected"]["covered"].items():
            assert per_ref[chrom]["covered_bases"] == covered

    def test_gap_count_per_reference(self, analyzer):
        per_ref = analyzer.stats["per_reference"]
        assert per_ref["chrA"]["gaps"] == 2   # 300 bp hole and 60 bp tail
        assert per_ref["chrB"]["gaps"] == 1   # 80 bp tail

    def test_low_identity_hit_is_excluded(self, analyzer):
        assert "L_00" not in {m["query_name"] for m in analyzer.mappings}


class TestGapAnalyzerRealGaps:
    def test_only_the_hole_passes_min_gap_size(self, analyzer, dataset):
        gap_analyzer = GapAnalyzer(
            coverage_data=analyzer.stats,
            reference_file=dataset["reference"],
            min_gap_size=100,
            coverage_arrays=analyzer.coverage_arrays,
            min_coverage=1.0,
        )
        results = gap_analyzer.analyze()
        chrom, start, end = dataset["expected"]["hole"]
        assert results["total_gaps"] == 1
        assert gap_analyzer.gaps[0]["chromosome"] == chrom
        assert gap_analyzer.gaps[0]["start"] == start
        assert gap_analyzer.gaps[0]["end"] == end
        assert results["total_gap_length"] == end - start

    def test_small_min_gap_size_finds_tails(self, analyzer, dataset):
        gap_analyzer = GapAnalyzer(
            coverage_data=analyzer.stats,
            reference_file=dataset["reference"],
            min_gap_size=50,
            coverage_arrays=analyzer.coverage_arrays,
        )
        results = gap_analyzer.analyze()
        assert results["total_gaps"] == 3

    def test_without_arrays_no_gaps_are_invented(self, analyzer, dataset):
        gap_analyzer = GapAnalyzer(coverage_data=analyzer.stats, reference_file=dataset["reference"])
        results = gap_analyzer.analyze()
        assert results["total_gaps"] == 0


class TestMergeUncoveredIntervals:
    def test_runs_are_merged_and_tails_kept_across_chromosomes(self):
        from baitUtils.coverage_analysis import merge_uncovered_intervals
        intervals = [
            ("chrA", 0, 100, 1), ("chrA", 100, 150, 0), ("chrA", 150, 200, 0),
            ("chrA", 200, 250, 2), ("chrA", 250, 300, 0),
            ("chrB", 0, 40, 3), ("chrB", 40, 60, 0), ("chrB", 60, 80, 1),
            ("genome", 0, 380, 0),
        ]
        regions, total = merge_uncovered_intervals(intervals, 1)
        assert regions == {"chrA": [(100, 200), (250, 300)], "chrB": [(40, 60)]}
        assert total == 170

    def test_no_uncovered(self):
        from baitUtils.coverage_analysis import merge_uncovered_intervals
        regions, total = merge_uncovered_intervals([("chrA", 0, 10, 5)], 1)
        assert regions == {} and total == 0


class TestNAwareBreadth:
    def test_breadth_excluding_n(self, tmp_path):
        # 1000 bp reference: 200 bp N run at 400-600, one bait covering 0-120
        ref = tmp_path / "ref.fa"
        ref.write_text(">chrN\n" + "ACGT" * 100 + "N" * 200 + "ACGT" * 100 + "\n")
        psl = tmp_path / "hits.psl"
        psl.write_text("\t".join(str(c) for c in [120, 0, 0, 0, 0, 0, 0, 0, "+", "b", 120, 0, 120,
                                                 "chrN", 1000, 0, 120, 1, "120,", "0,", "0,"]) + "\n")
        analyzer = CoverageAnalyzer(psl, ref, min_coverage=1.0)
        analyzer.analyze()
        stats = analyzer.stats
        assert stats["n_bases"] == 200
        assert stats["assessable_bases"] == 800
        assert stats["coverage_breadth"] == pytest.approx(12.0)
        assert stats["coverage_breadth_non_n"] == pytest.approx(15.0)
        per_ref = stats["per_reference"]["chrN"]
        assert per_ref["n_bases"] == 200
        assert per_ref["coverage_breadth_non_n"] == pytest.approx(15.0)

    def test_fixture_without_n_has_equal_breadths(self, analyzer):
        assert analyzer.stats["n_bases"] == 0
        assert analyzer.stats["coverage_breadth_non_n"] == pytest.approx(analyzer.stats["coverage_breadth"])


class TestCoverageExport:
    def test_export_writes_csv_bedgraph_and_bed(self, analyzer, tmp_path, dataset):
        coverage_df, gaps = analyzer.export_coverage_data(tmp_path)
        assert list(coverage_df.columns) == ["chromosome", "position", "coverage"]
        assert len(coverage_df) == 5000
        assert coverage_df["position"].iloc[0] == 1
        assert coverage_df.loc[coverage_df["chromosome"] == "chrA", "coverage"].sum() == int(analyzer.coverage_arrays["chrA"].sum())

        bedgraph = (tmp_path / "data" / "coverage.bedgraph").read_text().splitlines()
        rows = [line.split("\t") for line in bedgraph]
        # Runs are contiguous and cover each reference exactly once
        for chrom, length in (("chrA", 3000), ("chrB", 2000)):
            runs = [(int(s), int(e), int(d)) for c, s, e, d in rows if c == chrom]
            assert runs[0][0] == 0 and runs[-1][1] == length
            assert all(a[1] == b[0] for a, b in zip(runs, runs[1:]))
            assert all(a[2] != b[2] for a, b in zip(runs, runs[1:]))
        chrom, start, end = dataset["expected"]["hole"]
        assert (chrom, str(start), str(end), "0") in {tuple(r) for r in rows}

        bed = (tmp_path / "data" / "gap_regions.bed").read_text().splitlines()
        assert f"{chrom}\t{start}\t{end}" in bed
        assert len(gaps) == 3
