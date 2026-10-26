"""
Tests for phase 2 of the 2026-10 code audit: configuration file parsing,
kebab-case option names with hidden aliases for the old spellings, aligned
defaults, mapper-specific option warnings, and the new evaluate/compare
switches (--max-hits, --temp-dir, analysis parameters).
"""

import argparse
import json
import logging

import pandas as pd
import pytest

from baitUtils.__main__ import build_parser, _extract_config
from baitUtils.config import apply_config
from baitUtils.cli_utils import add_argument, register_alias
from baitUtils.sequence_mapping import resolve_mapper_options
from baitUtils.coverage_stats import CoverageAnalyzer


# ------------------------------------------------------------------ config

class TestConfigParsing:
    def test_equals_form_and_position_are_accepted(self, tmp_path):
        """F-8: --config=file and --config after the command were ignored."""
        cfg = tmp_path / "c.json"
        cfg.write_text("{}")
        assert _extract_config([f"--config={cfg}", "stats"]) == (cfg, ["stats"])
        assert _extract_config(["stats", "-i", "x", "--config", str(cfg)]) == (cfg, ["stats", "-i", "x"])
        assert _extract_config(["stats", "--min-dimer-dg", "-15"]) == (None, ["stats", "--min-dimer-dg", "-15"])

    def test_trailing_config_is_a_usage_error(self, capsys):
        with pytest.raises(SystemExit) as exc:
            _extract_config(["stats", "--config"])
        assert exc.value.code == 2
        assert "expected one argument" in capsys.readouterr().err

    def test_common_keys_warn_only_when_no_command_has_them(self, caplog):
        """F-9: the README example used to emit seven warnings."""
        _, commands = build_parser()
        with caplog.at_level(logging.WARNING):
            apply_config(commands, {"common": {"threads": 4, "log_level": "WARNING"},
                                    "evaluate": {"min_identity": 95.0}})
        assert "Configuration" not in caplog.text
        assert commands["stats"].parse_args(["-i", "a", "-o", "b"]).threads == 4
        assert commands["check"].parse_args(["--alignments", "a", "--reference", "b"]).log_level == "WARNING"
        assert commands["evaluate"].parse_args(["-i", "a", "-r", "b", "-o", "c"]).min_identity == 95.0

    def test_unknown_section_and_unused_common_key_warn(self, caplog):
        _, commands = build_parser()
        with caplog.at_level(logging.WARNING):
            apply_config(commands, {"evaluat": {"min_identity": 95}, "common": {"nothing": 1}})
        assert "sections that are not commands: evaluat" in caplog.text
        assert "used by no command: nothing" in caplog.text

    def test_config_in_equals_form_end_to_end(self, dataset, tmp_path, run_cli):
        cfg = tmp_path / "c.json"
        cfg.write_text(json.dumps({"stats": {"outdir": str(tmp_path / "out"), "filter": True}}))
        run_cli([f"--config={cfg}", "stats", "-i", dataset["baits"]])
        assert (tmp_path / "out" / "filtered_sequences.fasta").exists()


# ----------------------------------------------------------------- aliases

class TestAliases:
    def test_register_alias_resolves_to_the_same_action(self):
        parser = argparse.ArgumentParser()
        action = add_argument(parser, "--min-coverage", type=float, default=1.0, aliases=("--min_coverage",))
        assert parser.parse_args(["--min_coverage", "3"]).min_coverage == 3.0
        assert parser.parse_args(["--min-coverage", "4"]).min_coverage == 4.0
        assert "--min_coverage" not in parser.format_help()
        with pytest.raises(ValueError):
            register_alias(parser, action, "--min-coverage")

    @pytest.mark.parametrize("command, old, new", [
        ("stats", ["-i", "a", "-o", "b", "--mingc", "40", "--maxgc", "60", "--mint", "60", "--maxt", "80",
                   "--maxmask", "5", "--noNs", "--processes", "2", "--outdir", "c"],
                  ["-i", "a", "-o", "b", "--min-gc", "40", "--max-gc", "60", "--min-tm", "60", "--max-tm", "80",
                   "--max-mask", "5", "--no-ns", "--threads", "2", "--output", "c"]),
        ("plot", ["-i", "a", "--plot_type", "pca", "-l"], ["-i", "a", "--plot-type", "pca", "-v"]),
        ("map", ["-i", "a", "-q", "r", "-X", "2", "--minIdentity", "85", "--filterIdentity", "90",
                 "--minMatchCount", "3", "--minMatch", "1", "--minScore", "20", "--log"],
                ["-i", "a", "-r", "r", "--threads", "2", "--min-identity", "85", "--filter-identity", "90",
                 "--min-match-count", "3", "--min-match", "1", "--min-score", "20", "--verbose"]),
        ("check", ["--psl", "a", "--reference", "r", "--min_coverage", "2", "--min_similarity", "95",
                   "--min_length", "50", "--forced_oligos", "f", "--uncovered_length_cutoff", "10",
                   "--coverage_out", "c", "--longest_uncovered_out", "l", "--uncovered_fasta", "u",
                   "--n_split_fasta", "n", "--extend_region", "5", "--min_oligo_length", "100",
                   "--temp_dir", "t", "--log_level", "DEBUG"],
                  ["--alignments", "a", "--reference", "r", "--min-coverage", "2", "--min-identity", "95",
                   "--min-length", "50", "--forced-oligos", "f", "--uncovered-length-cutoff", "10",
                   "--coverage-out", "c", "--longest-uncovered-out", "l", "--uncovered-fasta", "u",
                   "--n-split-fasta", "n", "--extend-region", "5", "--min-oligo-length", "100",
                   "--temp-dir", "t", "--log-level", "DEBUG"]),
        ("fill", ["--psl", "a", "--reference", "r", "--spacing_distance", "10", "--min_contribution", "2",
                  "--max_passes", "3", "--max_oligos_per_pass", "7", "--stall_rounds", "2"],
                 ["--alignments", "a", "--reference", "r", "--spacing-distance", "10", "--min-contribution", "2",
                  "--max-passes", "3", "--max-oligos-per-pass", "7", "--stall-rounds", "2"]),
        ("evaluate", ["-i", "a", "-r", "r", "-o", "o", "--log"], ["-i", "a", "-r", "r", "-o", "o", "--verbose"]),
        ("compare", ["-r", "r", "-o", "o", "--sets", "A:a", "-l"], ["-r", "r", "-o", "o", "--sets", "A:a", "-v"]),
    ])
    def test_old_spellings_parse_like_the_new_ones(self, command, old, new):
        _, commands = build_parser()
        assert vars(commands[command].parse_args(old)) == vars(commands[command].parse_args(new))

    def test_help_shows_only_kebab_case(self):
        _, commands = build_parser()
        for name, parser in commands.items():
            help_text = parser.format_help()
            assert "_" not in "".join(line.split()[0] for line in help_text.splitlines()
                                      if line.strip().startswith("--")), name
            assert "--minIdentity" not in help_text and "--psl" not in help_text


# ---------------------------------------------------------------- defaults

class TestAlignedDefaults:
    def test_check_and_fill_share_evaluate_defaults(self):
        _, commands = build_parser()
        check = commands["check"].parse_args(["--alignments", "a", "--reference", "r"])
        fill = commands["fill"].parse_args(["--alignments", "a", "--reference", "r"])
        evaluate = commands["evaluate"].parse_args(["-i", "a", "-r", "r", "-o", "o"])
        for args in (check, fill):
            assert args.min_coverage == evaluate.min_coverage == 1.0
            assert args.min_identity == evaluate.min_identity == 90.0
            assert args.min_length == evaluate.min_length == 100
            assert args.longest_uncovered_out is None
            assert args.log_level is None and args.verbose is False and args.quiet is False

    def test_every_command_has_the_logging_options(self):
        _, commands = build_parser()
        for parser in commands.values():
            dests = {a.dest for a in parser._actions}
            assert {"verbose", "quiet", "log_level"} <= dests

    def test_validation_rejects_bad_values(self):
        _, commands = build_parser()
        for argv in (["-i", "a", "-r", "r", "-o", "o", "--threads", "0"],
                     ["-i", "a", "-r", "r", "-o", "o", "--min-gap-size", "0"],
                     ["-i", "a", "-r", "r", "-o", "o", "--plot-dpi", "-1"]):
            with pytest.raises(SystemExit):
                commands["evaluate"].parse_args(argv)
        with pytest.raises(SystemExit):
            commands["map"].parse_args(["-i", "a", "-r", "r", "--min-identity", "120"])


# ------------------------------------------------------------------ mappers

class TestMapperOptions:
    def _args(self, *extra):
        _, commands = build_parser()
        return commands["map"].parse_args(["-i", "a", "-r", "r", *extra])

    def test_pblat_options_warn_under_minimap2(self, caplog):
        """F-26"""
        args = self._args("--mapper", "minimap2", "--min-identity", "95", "--min-score", "20")
        with caplog.at_level(logging.WARNING):
            resolve_mapper_options(args)
        assert "--min-score, --min-identity apply to pblat only" in caplog.text
        assert args.minimap2_preset == "sr"

    def test_preset_warns_under_pblat_and_defaults_are_filled(self, caplog):
        args = self._args("--minimap2-preset", "map-ont")
        with caplog.at_level(logging.WARNING):
            resolve_mapper_options(args)
        assert "--minimap2-preset has no effect" in caplog.text
        assert (args.min_match, args.min_score, args.min_identity) == (2, 30, 90.0)

    def test_no_warning_for_matching_mapper(self, caplog):
        args = self._args("--mapper", "minimap2", "--minimap2-preset", "map-ont")
        with caplog.at_level(logging.WARNING):
            resolve_mapper_options(args)
        assert caplog.text == ""
        assert args.minimap2_preset == "map-ont"

    def test_filtered_file_keeps_the_mapper_format(self, dataset, tmp_path, fake_pblat, run_cli):
        """F-25: with pblat the filtered file is .psl (the fake mapper writes PSL)."""
        out = tmp_path / "map"
        run_cli(["map", "-i", dataset["baits"], "-r", dataset["reference"], "-o", out, "--prefix", "run"])
        assert (out / "run-mapping_filtered.psl").exists()


# ------------------------------------------------------- evaluate / compare

class TestEvaluateSwitches:
    def test_max_hits_excludes_multi_mapping_baits(self, dataset, tmp_path):
        """A second hit for A_00 makes it multi-mapping; --max-hits 1 drops it."""
        psl = tmp_path / "multi.psl"
        extra = "\t".join(str(c) for c in [120, 0, 0, 0, 0, 0, 0, 0, "+", "A_00", 120, 0, 120,
                                           "chrB", 2000, 1880, 2000, 1, "120,", "0,", "1880,"])
        psl.write_text(dataset["psl"].read_text() + extra + "\n")
        plain = CoverageAnalyzer(psl, dataset["reference"], oligos_file=dataset["baits"], min_length=0).analyze()
        strict = CoverageAnalyzer(psl, dataset["reference"], oligos_file=dataset["baits"], min_length=0,
                                  max_hits=1)
        strict_stats = strict.analyze()
        assert strict_stats["mapped_oligos"] == plain["mapped_oligos"] - 1
        assert strict.hit_table.set_index("oligo_id").loc["A_00", "multi_mapping_excluded"]
        assert strict_stats["per_reference"]["chrA"]["covered_bases"] == plain["per_reference"]["chrA"]["covered_bases"] - 120

    def test_evaluate_writes_hits_and_exports_and_scores_without_optional_steps(
            self, dataset, tmp_path, fake_pblat, run_cli):
        """F-17 and F-18."""
        out = tmp_path / "eval"
        run_cli(["evaluate", "-i", dataset["baits"], "-r", dataset["reference"], "-o", out,
                 "--no-html-report", "--no-interactive-plots", "--no-benchmarking",
                 "--temp-dir", tmp_path / "tmp", "--keep-intermediates", "--max-hits", "1"])
        hits = pd.read_csv(out / "hits.tsv", sep="\t")
        assert len(hits) == dataset["expected"]["n_mapped"]
        assert (out / "data" / "coverage.bedgraph").exists()
        # every run below the minimum depth: the 300 bp hole and the two uncovered tails
        assert (out / "data" / "gap_regions.bed").read_text().count("\n") == 3
        assert not (out / "data" / "coverage_per_position.csv").exists()
        assert (out / "intermediates" / "mapping.psl").exists()
        data = json.loads((out / "evaluation.json").read_text())
        assert data["quality_score"] is not None and data["reference_summary"] is not None
        assert data["arguments"]["max_hits"] == 1 and data["arguments"]["temp_dir"] == str(tmp_path / "tmp")

    def test_compare_passes_analysis_parameters(self, dataset, tmp_path, fake_pblat, run_cli):
        out = tmp_path / "cmp"
        run_cli(["compare", "-r", dataset["reference"], "-o", out,
                 "--sets", f"A:{dataset['baits']}", f"B:{dataset['baits']}",
                 "--min-gap-size", "400", "--reference-analysis-window", "500", "--max-hits", "2",
                 "--temp-dir", tmp_path / "tmp"])
        data = json.loads((out / "comparison.json").read_text())
        assert data["parameters"]["min_gap_size"] == 400
        assert data["parameters"]["reference_analysis_window"] == 500
        # The 300 bp hole is below the 400 bp floor, so no gap is reported
        assert all(s["gap_analysis"]["total_gaps"] == 0 for s in data["sets"])
        assert data["reference_analysis"]["coverage_correlations"]["window_size"] == 500
