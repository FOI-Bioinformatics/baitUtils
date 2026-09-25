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
        assert (out / "summary.txt").exists()
        assert (out / "filtered_sequences.fasta").read_text().count(">") == len(df)

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
                 "-o", "run", "-Z", outdir, "--filterIdentity", "90"])

        mapped = (outdir / "run-mapped-sequence-ids.txt").read_text().split()
        unmapped = (outdir / "run-unmapped-sequence-ids.txt").read_text().split()
        assert len(mapped) == dataset["expected"]["n_mapped"]
        assert len(unmapped) == dataset["expected"]["n_unmapped"]
        assert "L_00" in unmapped
        assert all(u in unmapped for u in ("U_00", "U_01", "U_02"))


class TestEvaluate:
    def test_evaluate_runs_end_to_end(self, dataset, tmp_path, fake_pblat, run_cli):
        out = tmp_path / "eval"
        run_cli(["evaluate", "-i", dataset["baits"], "-r", dataset["reference"], "-o", out])
        assert (out / "coverage_statistics.txt").exists()
        assert (out / "gap_analysis.txt").exists()
        assert (out / "coverage_evaluation_report.html").exists()


class TestCompare:
    def test_compare_runs_end_to_end(self, dataset, tmp_path, fake_pblat, run_cli):
        out = tmp_path / "cmp"
        run_cli(["compare", "-r", dataset["reference"], "-o", out,
                 "--sets", f"A:{dataset['baits']}", f"B:{dataset['baits']}"])
        assert (out / "comparison_matrix.csv").exists()
        assert (out / "comparative_analysis_report.html").exists()


@pytest.mark.skipif(not _has_bedtools(), reason="requires pybedtools and bedtools")
class TestCheckAndFill:
    def test_check_reports_hole(self, dataset, tmp_path, run_cli, monkeypatch):
        monkeypatch.chdir(tmp_path)
        uncovered = tmp_path / "uncovered.txt"
        run_cli(["check", "--psl", dataset["psl"], "--reference", dataset["reference"],
                 "--min_coverage", "1", "--longest_uncovered_out", uncovered])
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
        run_cli(["fill", "--psl", dataset["psl"], "--reference", dataset["reference"],
                 "--output", selected, "--min_contribution", "1"])
        assert selected.exists()
        assert len(selected.read_text().split()) > 0


class TestMissingDependencies:
    @pytest.mark.skipif(_has_bedtools(), reason="only meaningful without pybedtools")
    def test_check_without_bedtools_exits_with_message(self, dataset, tmp_path, run_cli, capsys, monkeypatch):
        monkeypatch.chdir(tmp_path)
        with pytest.raises(SystemExit) as exc:
            run_cli(["check", "--psl", dataset["psl"], "--reference", dataset["reference"]])
        assert exc.value.code != 0
        assert "bedtools" in str(exc.value.code).lower()
