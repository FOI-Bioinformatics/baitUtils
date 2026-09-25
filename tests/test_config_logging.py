"""
Tests for the configuration file, unified logging, the version record and
the strand options.
"""

import json
import logging

import pandas as pd
import pytest

from baitUtils.config import load_config, apply_config
from baitUtils.json_export import tool_versions, arguments_record
from baitUtils.logging_utils import configure_logging
from baitUtils.mapping_utils import filter_hits, parse_psl_line
from baitUtils.coverage_stats import CoverageAnalyzer


class TestConfigFile:
    def test_config_sets_defaults_and_lifts_required(self, dataset, tmp_path, run_cli):
        config = tmp_path / "cfg.json"
        config.write_text(json.dumps({
            "common": {"log_level": "WARNING"},
            "stats": {"input": str(dataset["baits"]), "outdir": str(tmp_path / "out"), "filter": True,
                      "unknown_key": 1},
        }))
        run_cli(["--config", config, "stats"])
        assert (tmp_path / "out" / "filtered_sequences.fasta").exists()

    def test_command_line_overrides_config(self, dataset, tmp_path, run_cli):
        config = tmp_path / "cfg.json"
        config.write_text(json.dumps({"stats": {"input": str(dataset["baits"]), "outdir": str(tmp_path / "a")}}))
        run_cli(["--config", config, "stats", "-o", tmp_path / "b"])
        assert (tmp_path / "b" / "sequence_statistics.tsv").exists()
        assert not (tmp_path / "a").exists()

    def test_apply_config_warns_on_unknown_keys(self, caplog):
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("--alpha", type=float, default=1.0)
        parser.add_argument("--needed", required=True)
        apply_config({"cmd": parser}, {"cmd": {"alpha": 2.5, "needed": "x", "bogus": 1}})
        args = parser.parse_args([])
        assert args.alpha == 2.5 and args.needed == "x"
        assert "bogus" in caplog.text

    def test_load_config_rejects_non_object(self, tmp_path):
        path = tmp_path / "bad.json"
        path.write_text("[1, 2]")
        with pytest.raises(ValueError):
            load_config(path)


class TestLogging:
    def test_ensure_logging_does_not_override_dispatcher(self):
        from baitUtils import logging_utils
        logging_utils._configured = False
        logging_utils.ensure_logging(level="ERROR")
        assert logging.getLogger().level == logging.ERROR
        logging_utils.ensure_logging(level="DEBUG")
        assert logging.getLogger().level == logging.ERROR
        configure_logging()

    def test_levels(self):
        configure_logging(level="ERROR")
        assert logging.getLogger().level == logging.ERROR
        configure_logging(level="ERROR", verbose=True)
        assert logging.getLogger().level == logging.DEBUG
        configure_logging(level="DEBUG", quiet=True)
        assert logging.getLogger().level == logging.WARNING
        configure_logging()
        assert logging.getLogger().level == logging.INFO

    def test_global_log_level_reaches_commands(self, dataset, tmp_path, run_cli):
        run_cli(["--log-level", "ERROR", "stats", "-i", dataset["baits"], "-o", tmp_path / "s"])
        assert logging.getLogger().level == logging.ERROR
        configure_logging()


class TestVersionsRecord:
    def test_tool_versions_has_core_entries(self):
        versions = tool_versions()
        for key in ("baitUtils", "python", "numpy", "pandas", "scipy", "biopython", "pblat", "minimap2", "bedtools"):
            assert key in versions
        assert versions["baitUtils"]

    def test_arguments_record_drops_callables(self):
        import argparse
        args = argparse.Namespace(a=1, func=lambda: None, path=None)
        record = arguments_record(args)
        assert record == {"a": 1, "path": None}

    def test_evaluation_json_carries_versions_and_arguments(self, dataset, tmp_path, fake_pblat, run_cli):
        out = tmp_path / "eval"
        run_cli(["evaluate", "-i", dataset["baits"], "-r", dataset["reference"], "-o", out,
                 "--no-html-report", "--no-interactive-plots", "--min-identity", "92"])
        data = json.loads((out / "evaluation.json").read_text())
        assert data["versions"]["baitUtils"]
        assert data["arguments"]["min_identity"] == 92.0
        assert "func" not in data["arguments"]


class TestStrand:
    def _hit(self, strand):
        cols = [120, 0, 0, 0, 0, 0, 0, 0, strand, "q", 120, 0, 120, "chrA", 3000, 0, 120, 1, "120,", "0,", "0,"]
        return parse_psl_line("\t".join(str(c) for c in cols))

    def test_filter_hits_by_strand(self):
        hits = [self._hit("+"), self._hit("-"), self._hit("+")]
        assert len(list(filter_hits(hits, strand="both"))) == 3
        assert [h.strand for h in filter_hits(hits, strand="plus")] == ["+", "+"]
        assert [h.strand for h in filter_hits(hits, strand="minus")] == ["-"]

    def test_coverage_analyzer_strand_counts(self, dataset, tmp_path):
        # Fixture hits are all on the plus strand; add one minus-strand hit
        psl = tmp_path / "mixed.psl"
        minus = "\t".join(str(c) for c in [120, 0, 0, 0, 0, 0, 0, 0, "-", "M_00", 120, 0, 120,
                                           "chrB", 2000, 1900, 2000, 1, "100,", "0,", "1900,"])
        psl.write_text(dataset["psl"].read_text() + minus + "\n")
        both = CoverageAnalyzer(psl, dataset["reference"], strand="both")
        both.analyze()
        assert both.stats["strand_counts"] == {"+": dataset["expected"]["n_mapped"], "-": 1}
        assert both.stats["per_reference"]["chrB"]["minus_hits"] == 1
        assert both.stats["per_reference"]["chrB"]["covered_bases"] == dataset["expected"]["covered"]["chrB"] + 80

        plus = CoverageAnalyzer(psl, dataset["reference"], strand="plus")
        plus.analyze()
        assert plus.stats["strand_counts"] == {"+": dataset["expected"]["n_mapped"], "-": 0}
        assert plus.stats["per_reference"]["chrB"]["covered_bases"] == dataset["expected"]["covered"]["chrB"]

        minus_only = CoverageAnalyzer(psl, dataset["reference"], strand="minus")
        minus_only.analyze()
        assert minus_only.stats["total_mappings"] == 1

    def test_map_orients_minus_strand_baits(self, dataset, tmp_path, monkeypatch, run_cli):
        from Bio import SeqIO
        from Bio.Seq import Seq
        # Bait file with one reverse-complemented tile, and a matching PSL on the minus strand
        records = {r.id: r for r in SeqIO.parse(str(dataset["baits"]), "fasta")}
        baits = tmp_path / "rc.fa"
        with open(baits, "w") as fh:
            fh.write(f">A_00\n{records['A_00'].seq}\n>A_01rc\n{Seq(str(records['A_01'].seq)).reverse_complement()}\n")
        psl = tmp_path / "rc.psl"
        rows = [
            [120, 0, 0, 0, 0, 0, 0, 0, "+", "A_00", 120, 0, 120, "chrA", 3000, 0, 120, 1, "120,", "0,", "0,"],
            [120, 0, 0, 0, 0, 0, 0, 0, "-", "A_01rc", 120, 0, 120, "chrA", 3000, 120, 240, 1, "120,", "0,", "120,"],
        ]
        psl.write_text("\n".join("\t".join(str(c) for c in r) for r in rows) + "\n")
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        script = bin_dir / "pblat"
        script.write_text(f"#!/bin/sh\nfor last; do :; done\ncp '{psl}' \"$last\"\n")
        script.chmod(0o755)
        monkeypatch.setenv("PATH", f"{bin_dir}:{__import__('os').environ['PATH']}")

        outdir = tmp_path / "map"
        outdir.mkdir()
        run_cli(["map", "-i", baits, "-q", dataset["reference"], "-o", outdir, "--prefix", "run",
                 "--orient-to-reference"])
        written = {r.id: str(r.seq) for r in SeqIO.parse(str(outdir / "run-mapped-sequences.fa"), "fasta")}
        assert written["A_00"] == str(records["A_00"].seq)
        assert written["A_01rc"] == str(records["A_01"].seq)   # back in reference orientation
        hits = pd.read_csv(outdir / "run-hits.tsv", sep="\t").set_index("oligo_id")
        assert hits.loc["A_01rc", "best_strand"] == "-"
