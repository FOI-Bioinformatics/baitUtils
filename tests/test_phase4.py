"""
Tests for phase 4 of the 2026-10 code audit: one-pass strand counts, the
shared reference loader, the fake minimap2 and the two-set compare
fixture, and the modules that had no tests (gap filling algorithm, plot
command edge cases).
"""

import json
import re

import numpy as np
import pandas as pd
import pytest

from baitUtils.coverage_stats import CoverageAnalyzer
from baitUtils.gap_analysis import GapAnalyzer
from baitUtils.reference_analyzer import ReferenceAnalyzer
from baitUtils.gap_filling_algorithm import (
    OligoMapping, OligoScorer, SpacingConstraintChecker, GreedySelector, MultiPassSelector,
    CoveragePatternAnalyzer,
)


# ------------------------------------------------------------ coverage_stats

class TestStrandCountsAndSharedLoader:
    def _mixed_psl(self, dataset, tmp_path):
        psl = tmp_path / "mixed.psl"
        rows = [dataset["psl"].read_text()]
        for i, (chrom, start) in enumerate([("chrA", 100), ("chrB", 300), ("chrB", 700)]):
            rows.append("\t".join(str(c) for c in [120, 0, 0, 0, 0, 0, 0, 0, "-", f"M_{i}", 120, 0, 120,
                                                    chrom, 3000 if chrom == "chrA" else 2000, start, start + 120,
                                                    1, "120,", "0,", f"{start},"]) + "\n")
        psl.write_text("".join(rows))
        return psl

    def test_per_reference_strand_counts_match_a_direct_count(self, dataset, tmp_path):
        """F-29: the counts are now gathered in one pass; they must not change."""
        analyzer = CoverageAnalyzer(self._mixed_psl(dataset, tmp_path), dataset["reference"], min_length=0)
        stats = analyzer.analyze()
        for ref_id, ref in stats["per_reference"].items():
            plus = sum(1 for m in analyzer.mappings if m["target_name"] == ref_id and m["strand"] == "+")
            minus = sum(1 for m in analyzer.mappings if m["target_name"] == ref_id and m["strand"] == "-")
            assert (ref["plus_hits"], ref["minus_hits"]) == (plus, minus)
        assert stats["per_reference"]["chrA"]["minus_hits"] == 1
        assert stats["per_reference"]["chrB"]["minus_hits"] == 2
        assert stats["strand_counts"] == {"+": 39, "-": 3}

    def test_analyzers_accept_loaded_sequences(self, dataset, tmp_path):
        """F-30: the gap and reference analyzers do not reread the FASTA when given the sequences."""
        analyzer = CoverageAnalyzer(dataset["psl"], dataset["reference"])
        stats = analyzer.analyze()
        sequences = analyzer.sequences()
        assert set(sequences) == {"chrA", "chrB"} and len(sequences["chrA"]) == 3000
        missing = tmp_path / "does_not_exist.fa"
        gaps = GapAnalyzer(stats, missing, coverage_arrays=analyzer.coverage_arrays,
                           reference_sequences=sequences).analyze()
        assert gaps["total_gaps"] == 1 and gaps["gaps"][0]["start"] == 1080
        ref = ReferenceAnalyzer(missing, stats, coverage_arrays=analyzer.coverage_arrays,
                                reference_sequences=sequences).analyze()
        assert ref["analysis_summary"]["total_sequences"] == 2
        assert ref["coverage_correlations"]["applicable"] is True

    def test_gap_features_are_case_insensitive(self, dataset, tmp_path):
        """The reference is kept as loaded; gap slices are upper-cased for the features."""
        analyzer = CoverageAnalyzer(dataset["psl"], dataset["reference"])
        stats = analyzer.analyze()
        lower = {k: v.lower() for k, v in analyzer.sequences().items()}
        upper = {k: v.upper() for k, v in analyzer.sequences().items()}
        f_lower = GapAnalyzer(stats, dataset["reference"], coverage_arrays=analyzer.coverage_arrays,
                              reference_sequences=lower).analyze()["gap_features"][0]
        f_upper = GapAnalyzer(stats, dataset["reference"], coverage_arrays=analyzer.coverage_arrays,
                              reference_sequences=upper).analyze()["gap_features"][0]
        assert f_lower == f_upper


# ------------------------------------------------------------------ fixtures

class TestFakeMappers:
    def test_fake_pblat_honours_the_query_fasta(self, dataset, tmp_path, fake_pblat, run_cli):
        out = tmp_path / "map"
        run_cli(["map", "-i", dataset["baits_subset"], "-r", dataset["reference"], "-o", out, "--prefix", "s"])
        mapped = (out / "s-mapped-sequence-ids.txt").read_text().split()
        assert len(mapped) == dataset["expected"]["subset"]["n_mapped"]
        assert "A_09" not in mapped and "A_08" in mapped

    def test_map_with_fake_minimap2_names_the_filtered_file_by_format(self, dataset, tmp_path, fake_minimap2, run_cli):
        """F-25"""
        out = tmp_path / "map"
        run_cli(["map", "-i", dataset["baits_subset"], "-r", dataset["reference"], "-o", out,
                 "--prefix", "s", "--mapper", "minimap2"])
        assert (out / "s-mapping.paf").exists() and (out / "s-mapping_filtered.paf").exists()
        assert not (out / "s-mapping_filtered.psl").exists()
        hits = pd.read_csv(out / "s-hits.tsv", sep="\t")
        assert len(hits) == dataset["expected"]["subset"]["n_mapped"]

    def test_evaluate_subset_has_the_wider_hole(self, dataset, tmp_path, fake_pblat, run_cli):
        out = tmp_path / "eval"
        run_cli(["evaluate", "-i", dataset["baits_subset"], "-r", dataset["reference"], "-o", out,
                 "--no-html-report", "--no-interactive-plots"])
        data = json.loads((out / "evaluation.json").read_text())
        exp = dataset["expected"]["subset"]
        assert data["coverage_stats"]["total_oligos"] == exp["n_baits"]
        assert data["coverage_stats"]["mapped_oligos"] == exp["n_mapped"]
        chrom, start, end = exp["hole"]
        gaps = {(g["chromosome"], g["start"], g["end"]) for g in data["gap_analysis"]["gaps"]}
        assert (chrom, start, end) in gaps
        assert data["coverage_stats"]["per_reference"]["chrA"]["covered_bases"] == exp["covered"]["chrA"]
        assert data["coverage_stats"]["per_reference"]["chrB"]["covered_bases"] == exp["covered"]["chrB"]


class TestCompareTwoDifferentSets:
    @pytest.fixture(scope="class")
    def comparison(self, dataset, tmp_path_factory):
        """compare on the full set and the subset; statistics then have something to test."""
        import os, stat, sys
        from baitUtils.__main__ import main
        out = tmp_path_factory.mktemp("cmp")
        bin_dir = out / "bin"
        bin_dir.mkdir()
        script = bin_dir / "pblat"
        awk = ('FNR==NR { if (substr($0,1,1)==">") { split(substr($0,2), f, " "); ids[f[1]]=1 }; next } '
               '($10 in ids)')
        script.write_text("#!/usr/bin/env bash\nif [ $# -eq 0 ]; then exit 0; fi\n"
                          'out="${@: -1}"\nquery="${@: -2:1}"\n'
                          f"awk '{awk}' \"$query\" '{dataset['psl']}' > \"$out\"\n")
        script.chmod(script.stat().st_mode | stat.S_IEXEC)
        old_path, old_argv = os.environ["PATH"], sys.argv
        os.environ["PATH"] = f"{bin_dir}{os.pathsep}{old_path}"
        sys.argv = ["baitUtils", "compare", "-r", str(dataset["reference"]), "-o", str(out / "result"),
                    "--sets", f"Full:{dataset['baits']}", f"Subset:{dataset['baits_subset']}",
                    f"Again:{dataset['baits']}", "--threads", "1"]
        try:
            main()
        finally:
            os.environ["PATH"], sys.argv = old_path, old_argv
        return out / "result"

    def test_statistics_are_applicable_and_rendered(self, comparison):
        data = json.loads((comparison / "comparison.json").read_text())
        per_ref = data["statistics"]["per_reference_metrics"]
        # three sets, two references that both differ between Full and Subset: Friedman runs
        assert per_ref["coverage_breadth"]["test_name"] == "Friedman"
        assert per_ref["coverage_breadth"]["applicable"] is True
        assert per_ref["coverage_breadth"]["p_adjusted"] is not None
        assert per_ref["coverage_breadth"]["significance_level"] in {"***", "**", "*", "ns"}
        assert per_ref["coverage_breadth"]["interpretation"].startswith(("Difference", "No difference"))
        html = (comparison / "comparative_analysis_report.html").read_text()
        assert "Friedman" in html
        assert re.search(r"P-value:</strong> \d\.\d{4} \(adjusted \d\.\d{4}\)", html)
        assert "Full" in html and "Subset" in html

    def test_ranking_prefers_the_full_set(self, comparison):
        data = json.loads((comparison / "comparison.json").read_text())
        ranking = dict(data["ranking"])
        assert ranking["Full"] > ranking["Subset"]
        assert ranking["Full"] == pytest.approx(ranking["Again"])
        sets = {s["name"]: s for s in data["sets"]}
        assert sets["Subset"]["coverage_stats"]["mapped_oligos"] == 33
        assert sets["Full"]["gap_analysis"]["total_gaps"] < sets["Subset"]["gap_analysis"]["total_gaps"]

    def test_gap_overlap_reports_the_subset_only_gap(self, comparison):
        text = (comparison / "gap_overlap_analysis.txt").read_text()
        # the chrA hole is shared by every set; only the subset's chrB tail is unique
        assert "Total unique gap regions: 3" in text
        assert "Subset: 1 unique gaps" in text and "chrB:1680-2000" in text
        assert "Full: 0 unique gaps" in text


# ------------------------------------------------------ gap filling algorithm

def _mapping(oligo, ref, start, end, hit="0"):
    return OligoMapping(oligo, ref, start, end, hit_id=hit)


def _coverage_calculator(ref_lengths):
    """Pure-Python stand-in for the bedtools-based calculator used by fill."""
    def calculate(selected, all_mappings, min_coverage):
        uncovered = {}
        total = 0
        for ref, length in ref_lengths.items():
            depth = np.zeros(length, dtype=int)
            for m in all_mappings.get(ref, []):
                if m.key in selected:
                    for s, e in m.blocks:
                        depth[s:e] += 1
            below = depth < min_coverage
            runs = []
            pos = 0
            while pos < length:
                if below[pos]:
                    start = pos
                    while pos < length and below[pos]:
                        pos += 1
                    runs.append((start, pos))
                    total += pos - start
                else:
                    pos += 1
            if runs:
                uncovered[ref] = runs
        return uncovered, total
    return calculate


class TestGapFillingAlgorithm:
    def test_spacing_constraint(self):
        positions = {"chr": [(0, 120), (300, 420)]}
        checker = SpacingConstraintChecker()
        assert checker.can_select_oligo(positions, _mapping("x", "chr", 150, 270), 30)
        assert not checker.can_select_oligo(positions, _mapping("x", "chr", 20, 140), 30)
        assert not checker.can_select_oligo(positions, _mapping("x", "chr", 290, 410), 30)
        assert checker.can_select_oligo({}, _mapping("x", "chr", 0, 120), 30)

    def test_scorer_rewards_overlap_with_uncovered_regions(self):
        uncovered = {"chr": [(100, 400)]}
        inside = OligoScorer.score_oligo(_mapping("a", "chr", 150, 270), uncovered, min_contribution=5)
        edge = OligoScorer.score_oligo(_mapping("b", "chr", 380, 500), uncovered, min_contribution=5)
        outside = OligoScorer.score_oligo(_mapping("c", "chr", 600, 720), uncovered, min_contribution=5)
        assert inside > edge > outside == 0.0
        assert OligoScorer.score_oligo(_mapping("b", "chr", 398, 518), uncovered, min_contribution=5) == 0.0

    def test_difficulty_is_higher_with_fewer_options(self):
        uncovered = {"chr": [(0, 100), (500, 600)]}
        mappings = {"chr": [_mapping("a", "chr", 0, 120), _mapping("b", "chr", 10, 130), _mapping("c", "chr", 500, 620)]}
        scores = CoveragePatternAnalyzer.analyze_coverage_patterns(uncovered, mappings)
        assert scores["chr:500-600"] > scores["chr:0-100"]

    def test_multi_pass_selection_closes_the_gaps(self):
        length = 1000
        mappings = {"chr": [_mapping(f"o{i}", "chr", s, s + 120) for i, s in enumerate(range(0, 900, 100))]}
        selector = MultiPassSelector()
        selected = selector.multi_pass_selection(
            mappings, forced_oligos=set(), coverage_calculator=_coverage_calculator({"chr": length}),
            min_coverage=1, spacing_distance=30, min_contribution=5, max_passes=10)
        uncovered, total = _coverage_calculator({"chr": length})(selected, mappings, 1)
        # every base before the last oligo's end is covered; only the tail after 920 remains
        assert total == length - 920
        assert uncovered == {"chr": [(920, 1000)]}
        assert {k[0] for k in selected} == {f"o{i}" for i in range(9)}

    def test_forced_oligos_are_selected_first_and_spacing_is_respected(self):
        mappings = {"chr": [_mapping("a", "chr", 0, 120), _mapping("b", "chr", 10, 130), _mapping("c", "chr", 300, 420)]}
        selector = MultiPassSelector()
        selected = selector.multi_pass_selection(
            mappings, forced_oligos={"b"}, coverage_calculator=_coverage_calculator({"chr": 500}),
            min_coverage=1, spacing_distance=30, min_contribution=5, max_passes=5)
        names = {k[0] for k in selected}
        assert "b" in names and "c" in names
        assert "a" not in names          # starts within 30 bp of the forced oligo b

    def test_single_pass_respects_max_oligos_per_pass(self):
        mappings = {"chr": [_mapping(f"o{i}", "chr", s, s + 120) for i, s in enumerate(range(0, 900, 200))]}
        uncovered = {"chr": [(0, 1000)]}
        chosen = GreedySelector().single_pass_greedy(mappings, uncovered, set(), 30, 5, max_oligos_per_pass=2)
        assert len(chosen) == 2


# --------------------------------------------------------------------- plot

class TestPlotCommand:
    def test_pca_survives_missing_values_and_svg_output(self, dataset, tmp_path, run_cli):
        """F-39: a NaN used to abort the whole plot command."""
        stats_out = tmp_path / "stats"
        run_cli(["stats", "-i", dataset["baits"], "-o", stats_out, "--threads", "1"])
        table = stats_out / "sequence_statistics.tsv"
        df = pd.read_csv(table, sep="\t")
        df.loc[0, "gc_content"] = np.nan
        df.to_csv(table, sep="\t", index=False)
        plots = tmp_path / "plots"
        run_cli(["plot", "-i", table, "-o", plots, "--plot-type", "pca", "histogram", "--format", "svg",
                 "--columns", "gc_content", "melting_temperature", "entropy", "--color", "kept"])
        assert (plots / "pca_scores_plot.svg").exists()
        assert (plots / "histogram_gc_content.svg").exists()
        assert not list(plots.glob("*.png"))
