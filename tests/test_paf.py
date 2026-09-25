"""
Tests for PAF (minimap2) input: parsing, equivalence with PSL on the fixture
dataset, and the mapper dispatch.
"""

import shutil

import pytest

from baitUtils.coverage_stats import CoverageAnalyzer
from baitUtils.mapping_utils import (
    parse_paf_line, parse_paf, parse_alignments, alignment_format, parse_psl,
)


def paf(*cols):
    return "\t".join(str(c) for c in cols)


PERFECT = paf("b0", 120, 0, 120, "+", "chrA", 3000, 0, 120, 120, 120, 60, "NM:i:0", "tp:A:P", "cg:Z:120M")
MISMATCH = paf("L", 120, 0, 120, "+", "chrA", 3000, 500, 620, 100, 120, 60, "NM:i:20", "cg:Z:120M")
T_INSERT = paf("G", 120, 0, 120, "+", "chrB", 2000, 500, 630, 120, 130, 60, "NM:i:10", "cg:Z:60M10D60M")
Q_INSERT = paf("Q", 120, 0, 120, "-", "chrB", 2000, 100, 210, 110, 120, 60, "NM:i:10", "cg:Z:50M10I60M")
NO_CIGAR = paf("b1", 120, 2, 116, "+", "chrA", 3000, 202, 316, 114, 114, 60, "tp:A:P")


class TestParsePafLine:
    def test_simple_hit(self):
        hit = parse_paf_line(PERFECT)
        assert (hit.q_name, hit.t_name, hit.t_start, hit.t_end, hit.strand) == ("b0", "chrA", 0, 120, "+")
        assert hit.matches == 120 and hit.mismatches == 0
        assert hit.block_sizes == [120] and hit.t_starts == [0] and hit.q_starts == [0]
        assert hit.identity == pytest.approx(100.0)

    def test_mismatches_from_nm(self):
        hit = parse_paf_line(MISMATCH)
        assert hit.mismatches == 20
        assert hit.identity == pytest.approx(100 - 1000 * 20 / 120 / 10)

    def test_target_insert_gives_two_blocks(self):
        hit = parse_paf_line(T_INSERT)
        assert hit.target_blocks == [(500, 560), (570, 630)]
        assert hit.q_starts == [0, 60]
        assert (hit.t_num_insert, hit.t_base_insert, hit.mismatches) == (1, 10, 0)
        assert hit.aligned_length == 120 and hit.target_span == 130
        assert hit.identity == pytest.approx(100.0)

    def test_query_insert_on_minus_strand(self):
        hit = parse_paf_line(Q_INSERT)
        assert hit.strand == "-"
        assert hit.target_blocks == [(100, 150), (150, 210)]
        assert (hit.q_num_insert, hit.q_base_insert, hit.mismatches) == (1, 10, 0)
        assert hit.identity < 100.0

    def test_without_cigar_single_block_and_estimated_mismatches(self):
        hit = parse_paf_line(NO_CIGAR)
        assert hit.target_blocks == [(202, 316)]
        assert hit.matches == 114 and hit.mismatches == 0

    def test_bad_lines(self):
        assert parse_paf_line("b0\t120\t0") is None
        assert parse_paf_line(PERFECT.replace("\t+\t", "\t?\t")) is None
        assert parse_paf_line(PERFECT.replace("\t120\t0\t120\t", "\tx\t0\t120\t", 1)) is None


class TestDispatch:
    def test_alignment_format(self):
        assert alignment_format("a/b.paf") == "paf"
        assert alignment_format("a/b.PAF.gz") == "paf"
        assert alignment_format("a/b.psl") == "psl"
        assert alignment_format("mapping") == "psl"

    def test_parse_alignments_reads_both(self, dataset):
        psl_hits = list(parse_alignments(dataset["psl"]))
        paf_hits = list(parse_alignments(dataset["paf"]))
        assert [h.q_name for h in psl_hits] == [h.q_name for h in paf_hits]
        for a, b in zip(psl_hits, paf_hits):
            assert a.target_blocks == b.target_blocks
            assert a.identity == pytest.approx(b.identity)

    def test_paf_comments_skipped(self, tmp_path, caplog):
        path = tmp_path / "x.paf"
        path.write_text("# comment\n" + PERFECT + "\nbad\n")
        assert [h.q_name for h in parse_paf(path)] == ["b0"]
        assert "Skipped 1 malformed" in caplog.text


class TestCoverageFromPaf:
    def test_same_statistics_as_psl(self, dataset):
        def stats(path):
            analyzer = CoverageAnalyzer(psl_file=path, reference_file=dataset["reference"],
                                        min_coverage=1.0, oligos_file=dataset["baits"])
            analyzer.analyze()
            return analyzer.stats

        psl_stats, paf_stats = stats(dataset["psl"]), stats(dataset["paf"])
        assert paf_stats["mapped_oligos"] == psl_stats["mapped_oligos"] == dataset["expected"]["n_mapped"]
        for chrom in ("chrA", "chrB"):
            assert paf_stats["per_reference"][chrom]["covered_bases"] == psl_stats["per_reference"][chrom]["covered_bases"]


class TestMapperCli:
    def test_map_with_fake_minimap2(self, dataset, tmp_path, fake_minimap2, run_cli):
        import pandas as pd
        outdir = tmp_path / "map"
        outdir.mkdir()
        run_cli(["map", "-i", dataset["baits"], "-q", dataset["reference"], "-o", outdir,
                 "--prefix", "run", "--mapper", "minimap2", "--filterIdentity", "90"])
        assert (outdir / "run-mapping.paf").exists()
        mapped = (outdir / "run-mapped-sequence-ids.txt").read_text().split()
        assert len(mapped) == dataset["expected"]["n_mapped"]
        hits = pd.read_csv(outdir / "run-hits.tsv", sep="\t")
        assert hits.set_index("oligo_id").loc["G_00", "best_end"] == 630

    def test_evaluate_with_fake_minimap2(self, dataset, tmp_path, fake_minimap2, run_cli):
        import json
        out = tmp_path / "eval"
        run_cli(["evaluate", "-i", dataset["baits"], "-r", dataset["reference"], "-o", out,
                 "--mapper", "minimap2", "--no-html-report", "--no-interactive-plots"])
        data = json.loads((out / "evaluation.json").read_text())
        assert data["coverage_stats"]["mapped_oligos"] == dataset["expected"]["n_mapped"]
        assert data["gap_analysis"]["total_gaps"] == 1

    def test_check_accepts_paf(self, dataset, tmp_path, run_cli, monkeypatch):
        pytest.importorskip("pybedtools")
        if shutil.which("bedtools") is None:
            pytest.skip("requires bedtools")
        monkeypatch.chdir(tmp_path)
        uncovered = tmp_path / "unc.tsv"
        run_cli(["check", "--alignments", dataset["paf"], "--reference", dataset["reference"],
                 "--min_coverage", "1", "--longest_uncovered_out", uncovered])
        rows = {tuple(line.split()[:3]) for line in uncovered.read_text().splitlines()[1:]}
        chrom, start, end = dataset["expected"]["hole"]
        assert (chrom, str(start), str(end)) in rows


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="requires minimap2")
class TestRealMinimap2:
    def test_evaluate_with_real_minimap2(self, dataset, tmp_path, run_cli):
        import json
        out = tmp_path / "eval_mm2"
        run_cli(["evaluate", "-i", dataset["baits"], "-r", dataset["reference"], "-o", out,
                 "--mapper", "minimap2", "--min-identity", "90", "--no-html-report", "--no-interactive-plots"])
        data = json.loads((out / "evaluation.json").read_text())
        # Tiles map exactly; the 20-mismatch bait and the random baits do not pass
        assert data["coverage_stats"]["mapped_oligos"] >= dataset["expected"]["n_mapped"] - 1
        assert data["coverage_stats"]["mapped_oligos"] <= dataset["expected"]["n_mapped"]
        assert data["gap_analysis"]["total_gaps"] == 1
