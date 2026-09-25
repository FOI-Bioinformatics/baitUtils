"""
End-to-end tests that drive every subcommand through the CLI dispatcher.

pblat is replaced by a shell script that writes the fixture PSL, so these
tests exercise the real subprocess seam without the external tool. check and
fill require pybedtools and the bedtools binary and are skipped otherwise.
"""

import shutil

import pandas as pd
import pytest


def _has_bedtools() -> bool:
    try:
        import pybedtools  # noqa: F401
    except ImportError:
        return False
    return shutil.which("bedtools") is not None


class TestStats:
    def test_stats_writes_tsv_summary_and_filtered_fasta(self, dataset, tmp_path, run_cli):
        out = tmp_path / "stats"
        run_cli(["stats", "-i", dataset["baits"], "-o", out, "--filter"])

        tsv = out / "sequence_statistics.tsv"
        assert tsv.exists()
        df = pd.read_csv(tsv, sep="\t")
        assert len(df) == dataset["expected"]["n_baits"]
        assert (df["length"] == 120).all()
        assert {"hairpin_dg", "self_dimer_dg", "melting_temperature"} <= set(df.columns)
        assert "mfe" not in df.columns
        assert (out / "summary.txt").exists()
        assert (out / "filtered_sequences.fasta").read_text().count(">") == len(df)

    def test_stats_dimer_filter_removes_palindromes(self, dataset, tmp_path, run_cli):
        pytest.importorskip("RNA")
        fasta = tmp_path / "baits.fa"
        # A/C-only sequence cannot pair with itself; the palindrome forms a full duplex
        fasta.write_text(">plain\n" + "AACCCAACAC" * 12 + "\n>palindrome\n" + "GAATTCGGATCCGAATTC" * 6 + "GAATTCGGATCC\n")
        out = tmp_path / "stats"
        run_cli(["stats", "-i", fasta, "-o", out, "--filter", "--min-dimer-dg", "-15"])
        df = pd.read_csv(out / "sequence_statistics.tsv", sep="\t")
        kept = dict(zip(df["sequence_id"], df["kept"]))
        assert kept["palindrome"] == False  # noqa: E712
        assert kept["plain"] == True  # noqa: E712

    def test_stats_reports_melting_temperature(self, dataset, tmp_path, run_cli):
        out = tmp_path / "stats"
        run_cli(["stats", "-i", dataset["baits"], "-o", out])
        df = pd.read_csv(out / "sequence_statistics.tsv", sep="\t")
        assert pd.to_numeric(df["melting_temperature"], errors="coerce").notna().all()


class TestPlot:
    def test_plot_from_stats_tsv(self, dataset, tmp_path, run_cli):
        stats_out = tmp_path / "stats"
        run_cli(["stats", "-i", dataset["baits"], "-o", stats_out])
        plots = tmp_path / "plots"
        run_cli(["plot", "-i", stats_out / "sequence_statistics.tsv", "-o", plots,
                 "--columns", "gc_content", "entropy",
                 "--plot_type", "histogram", "scatterplot"])
        assert any(plots.glob("*.png"))


class TestMap:
    def test_map_splits_mapped_and_unmapped(self, dataset, tmp_path, fake_pblat, run_cli):
        outdir = tmp_path / "map"
        outdir.mkdir()
        run_cli(["map", "-i", dataset["baits"], "-q", dataset["reference"],
                 "-o", outdir, "--prefix", "run", "--filterIdentity", "90"])

        mapped = (outdir / "run-mapped-sequence-ids.txt").read_text().split()
        unmapped = (outdir / "run-unmapped-sequence-ids.txt").read_text().split()
        assert len(mapped) == dataset["expected"]["n_mapped"]
        assert len(unmapped) == dataset["expected"]["n_unmapped"]
        assert "L_00" in unmapped
        assert all(u in unmapped for u in ("U_00", "U_01", "U_02"))

        hits = pd.read_csv(outdir / "run-hits.tsv", sep="\t")
        assert len(hits) == dataset["expected"]["n_mapped"]
        assert (hits["n_hits"] == 1).all()
        assert hits.set_index("oligo_id").loc["G_00", "best_target"] == "chrB"

    def test_map_max_hits_excludes_multi_mapping_baits(self, dataset, tmp_path, fake_pblat, run_cli):
        outdir = tmp_path / "map"
        outdir.mkdir()
        run_cli(["map", "-i", dataset["baits"], "-q", dataset["reference"],
                 "-o", outdir, "--prefix", "run", "--max-hits", "0"])
        mapped = (outdir / "run-mapped-sequence-ids.txt").read_text().split()
        assert mapped == []


class TestEvaluate:
    def test_evaluate_runs_end_to_end(self, dataset, tmp_path, fake_pblat, run_cli):
        out = tmp_path / "eval"
        run_cli(["evaluate", "-i", dataset["baits"], "-r", dataset["reference"], "-o", out])
        assert (out / "coverage_statistics.txt").exists()
        assert (out / "gap_analysis.txt").exists()
        assert (out / "coverage_evaluation_report.html").exists()

        import json
        data = json.loads((out / "evaluation.json").read_text())
        exp = dataset["expected"]
        assert data["coverage_stats"]["total_oligos"] == exp["n_baits"]
        assert data["coverage_stats"]["mapped_oligos"] == exp["n_mapped"]
        assert data["gap_analysis"]["total_gaps"] == 1
        chrom, start, end = exp["hole"]
        gap = data["gap_analysis"]["gaps"][0]
        assert (gap["chromosome"], gap["start"], gap["end"]) == (chrom, start, end)
        assert 0.0 <= data["quality_score"]["overall_score"] <= 1.0
        assert set(data["benchmarks"]) == {"coverage_breadth", "depth_uniformity", "gap_reduction"}


class TestCompare:
    def test_compare_runs_end_to_end(self, dataset, tmp_path, fake_pblat, run_cli):
        out = tmp_path / "cmp"
        run_cli(["compare", "-r", dataset["reference"], "-o", out,
                 "--sets", f"A:{dataset['baits']}", f"B:{dataset['baits']}"])
        assert (out / "comparison_matrix.csv").exists()
        report = (out / "comparative_analysis_report.html").read_text()
        assert "P-value:</strong> nan" not in report
        assert "Effect Size:</strong> nan" not in report

        import json
        data = json.loads((out / "comparison.json").read_text())
        assert [s_["name"] for s_ in data["sets"]] == ["A", "B"]
        assert data["statistics"]["per_reference_metrics"]["coverage_breadth"]["applicable"] is False
        assert data["statistics"]["bait_identity"]["applicable"] is False  # identical sets
        assert data["parameters"]["multiple_comparison_correction"] == "fdr"
        assert "Not applicable" in report  # identical sets: paired tests cannot run


@pytest.mark.skipif(not _has_bedtools(), reason="requires pybedtools and bedtools")
class TestCheckAndFill:
    def test_check_reports_hole(self, dataset, tmp_path, run_cli, monkeypatch):
        monkeypatch.chdir(tmp_path)
        uncovered = tmp_path / "uncovered.txt"
        before = set(tmp_path.iterdir())
        run_cli(["check", "--psl", dataset["psl"], "--reference", dataset["reference"],
                 "--min_coverage", "1", "--longest_uncovered_out", uncovered])
        assert set(tmp_path.iterdir()) - before == {uncovered}
        rows = [line.split() for line in uncovered.read_text().splitlines()
                if line.strip() and not line.lower().startswith(("reference", "chrom", "#"))]
        regions = {(r[0], int(r[1]), int(r[2])) for r in rows}
        assert dataset["expected"]["hole"] in regions
        # Tails beyond the last mapped bait are reported because sizes come from the FASTA
        assert ("chrA", 2940, 3000) in regions
        assert ("chrB", 1920, 2000) in regions

    def test_fill_selects_oligos(self, dataset, tmp_path, run_cli, monkeypatch):
        monkeypatch.chdir(tmp_path)
        selected = tmp_path / "selected.txt"
        before = set(tmp_path.iterdir())
        run_cli(["fill", "--psl", dataset["psl"], "--reference", dataset["reference"],
                 "--output", selected, "--min_contribution", "1"])
        assert selected.exists()
        assert len(selected.read_text().split()) > 0
        mappings = pd.read_csv(tmp_path / "selected_mappings.tsv", sep="\t")
        assert set(mappings.columns) == {"oligo_id", "reference", "start", "end"}
        assert set(mappings["oligo_id"]) == set(selected.read_text().split())
        # No temporary files are left in the working directory
        assert set(tmp_path.iterdir()) - before == {selected, tmp_path / "selected_mappings.tsv"}


class TestMissingDependencies:
    @pytest.mark.skipif(_has_bedtools(), reason="only meaningful without pybedtools")
    def test_check_without_bedtools_exits_with_message(self, dataset, tmp_path, run_cli, capsys, monkeypatch):
        monkeypatch.chdir(tmp_path)
        with pytest.raises(SystemExit) as exc:
            run_cli(["check", "--psl", dataset["psl"], "--reference", dataset["reference"]])
        assert exc.value.code != 0
        assert "bedtools" in str(exc.value.code).lower()
