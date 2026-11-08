"""
Regression tests for phase 1 of the second audit (docs/code_audit_2026-11.md):
hangs, crashes and silent failures on realistic input. Each test failed
before its fix.
"""

import json
import time

import numpy as np
import pytest

from baitUtils.interactive_plots import InteractivePlotter


# ------------------------------------------------------------------ N-51

class TestCoverageProfileShapes:
    def test_thousands_of_gaps_build_quickly(self, tmp_path):
        """N-51: one add_vrect per gap was quadratic (about an hour for 2800 gaps)."""
        length, n_gaps = 1_000_000, 3000
        gaps = [{'chromosome': 'big', 'start': i * 300, 'end': i * 300 + 150, 'length': 150}
                for i in range(n_gaps)]
        plotter = InteractivePlotter({'per_reference': {'big': {'length': length}}, 'parameters': {}},
                                     {'gaps': gaps}, {}, {}, tmp_path,
                                     coverage_arrays={'big': np.ones(length)})
        start = time.time()
        fig = plotter.coverage_profile_figure()
        assert time.time() - start < 20
        rects = [s for s in fig.layout.shapes if s.type == 'rect']
        assert len(rects) == n_gaps
        assert (rects[0].x0, rects[0].x1) == (0, 150)

    def test_shapes_sit_on_their_own_subplot(self, tmp_path):
        arrays = {'long': np.ones(5000), 'short': np.ones(4000)}
        gaps = [{'chromosome': 'short', 'start': 10, 'end': 500, 'length': 490}]
        plotter = InteractivePlotter(
            {'per_reference': {'long': {'length': 5000}, 'short': {'length': 4000}}, 'parameters': {}},
            {'gaps': gaps}, {}, {}, tmp_path, coverage_arrays=arrays)
        fig = plotter.coverage_profile_figure()
        rect = [s for s in fig.layout.shapes if s.type == 'rect'][0]
        assert (rect.xref, rect.yref) == ('x2', 'y2 domain')     # second subplot: the shorter reference
        lines = [s for s in fig.layout.shapes if s.type == 'line']
        assert {(s.xref, s.yref) for s in lines} == {('x domain', 'y'), ('x2 domain', 'y2')}


# ------------------------------------------------------------------ N-52

def _fasta(path, records):
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
    return path


class TestFigureRobustness:
    def test_dollar_in_reference_name_does_not_crash(self, tmp_path):
        """N-52: '$...$' in a reference name was parsed as mathtext and aborted evaluate."""
        from baitUtils.coverage_viz import CoverageVisualizer
        names = ["ctg$^$x", "seq$1$_a"]
        stats = {'per_reference': {n: {'length': 1000, 'coverage_breadth': 80.0, 'mean_depth': 0.8,
                                       'max_depth': 1, 'gaps': 1} for n in names},
                 'coverage_breadth': 80.0, 'mean_depth': 0.8, 'parameters': {}}
        gaps = {'gaps': [{'chromosome': n, 'start': 100, 'end': 300, 'length': 200} for n in names]}
        arrays = {n: np.r_[np.ones(100), np.zeros(200), np.ones(700)] for n in names}
        written = CoverageVisualizer(stats, gaps, tmp_path, coverage_arrays=arrays).generate_all_plots()
        assert {'coverage_profile', 'per_reference_coverage', 'coverage_heatmap'} <= set(written)

    def test_one_failing_figure_does_not_abort_the_others(self, tmp_path, monkeypatch):
        from baitUtils.coverage_viz import CoverageVisualizer
        stats = {'per_reference': {'chrA': {'length': 100, 'coverage_breadth': 50.0, 'mean_depth': 1.0,
                                            'max_depth': 2, 'gaps': 1}},
                 'coverage_breadth': 50.0, 'mean_depth': 1.0, 'parameters': {}}
        viz = CoverageVisualizer(stats, {'gaps': []}, tmp_path, coverage_arrays={'chrA': np.ones(100)})

        def boom():
            raise ValueError("broken figure")
        monkeypatch.setattr(viz, 'plot_coverage_profile', boom)
        written = viz.generate_all_plots()
        assert 'coverage_profile' not in written
        assert 'coverage_overview' in written

    def test_style_shows_dollar_literally(self):
        import matplotlib
        from baitUtils import design
        assert design.matplotlib_rc()["text.parse_math"] is False
        assert matplotlib.rcParams["text.parse_math"] is True      # global style untouched


# ------------------------------------------------------------------- N-1

PSLAYOUT_HEADER = (
    "psLayout version 3\n"
    "\n"
    "match\tmis- \trep. \tN's\tQ gap\tQ gap\tT gap\tT gap\tstrand\tQ        \tQ   \tQ    \tQ  \tT        "
    "\tT   \tT    \tT  \tblock\tblockSizes \tqStarts\t tStarts\n"
    "     \tmatch\tmatch\t   \tcount\tbases\tcount\tbases\t      \tname     \tsize\tstart\tend\tname     "
    "\tsize\tstart\tend\tcount\n"
    "---------------------------------------------------------------------------------------------------"
    "------------------------------------------------------------\n"
)


class TestPslLayoutHeader:
    def _with_header(self, dataset, tmp_path):
        path = tmp_path / "pblat.psl"
        path.write_text(PSLAYOUT_HEADER + dataset["psl"].read_text())
        return path

    def test_pblat_header_is_not_reported_as_malformed(self, dataset, tmp_path, caplog):
        """N-1: the header line that starts with spaces was counted as a malformed record."""
        from baitUtils.mapping_utils import parse_psl
        path = self._with_header(dataset, tmp_path)
        with caplog.at_level("WARNING"):
            hits = list(parse_psl(path))
        assert len(hits) == len(list(parse_psl(dataset["psl"])))
        assert "malformed" not in caplog.text

    def test_filtered_psl_keeps_the_whole_header(self, dataset, tmp_path):
        from baitUtils.mapping_utils import PSLParser
        path = self._with_header(dataset, tmp_path)
        out = tmp_path / "filtered.psl"
        PSLParser.parse_psl_file(str(path), min_identity=90.0, filtered_output=str(out))
        lines = out.read_text().splitlines()
        assert "\n".join(lines[:5]) + "\n" == PSLAYOUT_HEADER
        assert lines[5].split("\t")[0].isdigit()

    def test_malformed_record_after_header_is_still_counted(self, dataset, tmp_path, caplog):
        from baitUtils.mapping_utils import parse_psl
        path = tmp_path / "bad.psl"
        path.write_text(PSLAYOUT_HEADER + dataset["psl"].read_text() + "12\tnot\ta\trecord\n")
        with caplog.at_level("WARNING"):
            list(parse_psl(path))
        assert "Skipped 1 malformed line" in caplog.text


# ------------------------------------------------------------------ N-32

class TestEmptyMapperOutput:
    def test_minimap2_with_no_alignments_is_not_a_failure(self, dataset, tmp_path, run_cli, fake_minimap2):
        """N-32: minimap2 writes an empty PAF when nothing maps; evaluate exited 1."""
        baits = _fasta(tmp_path / "unknown.fa", [("NOT_IN_FIXTURE", "ACGT" * 30)])
        out = tmp_path / "eval"
        run_cli(["evaluate", "-i", baits, "-r", dataset["reference"], "-o", out, "--mapper", "minimap2",
                 "--no-interactive-plots", "--no-html-report"])
        data = json.loads((out / "evaluation.json").read_text())
        assert data["coverage_stats"]["mapped_oligos"] == 0
        assert data["coverage_stats"]["coverage_breadth"] == 0.0


# ---------------------------------------------------------- N-14 to N-16

def _has_bedtools() -> bool:
    import shutil
    try:
        import pybedtools  # noqa: F401
    except ImportError:
        return False
    return shutil.which("bedtools") is not None


class TestFillStopping:
    def _selector_run(self, totals, regions, stall_rounds=3, max_passes=20):
        """Run the multi-pass loop with a scripted coverage calculator."""
        from baitUtils.gap_filling_algorithm import MultiPassSelector, OligoMapping
        mappings = {'chrA': [OligoMapping(f'b{i}', 'chrA', i * 200, i * 200 + 120) for i in range(30)]}
        calls = []

        def coverage(selected, all_mappings, min_coverage):
            i = min(len(calls), len(totals) - 1)
            calls.append(i)
            # one 50 bp region at the start of each of the first regions[i] baits,
            # so unselected candidates never run out
            return {'chrA': [(200 * j, 200 * j + 50) for j in range(regions[i])]}, totals[i]

        selector = MultiPassSelector()
        selector.multi_pass_selection(mappings, set(), coverage, 1.0, 0, 1, max_passes,
                                      max_oligos_per_pass=1, stall_rounds=stall_rounds)
        return len(calls)

    def test_splitting_gaps_is_progress(self):
        """N-14: more uncovered regions with fewer uncovered bases counted as a stall."""
        totals = [5000, 4400, 3800, 3200, 2600, 2000, 1400]
        regions = [2, 6, 10, 10, 12, 12, 12]
        # seven passes with fewer bases each, then three without: ten coverage calls.
        # The region-count rule stopped at the fourth call (10 regions, as before).
        assert self._selector_run(totals, regions) == 10

    def test_stops_after_stall_rounds_without_fewer_bases(self):
        totals = [5000, 4000, 4000, 4000, 4000, 4000]
        regions = [10, 10, 10, 10, 10, 10]
        # 5000, 4000, then three passes without fewer bases: stop on the third
        assert self._selector_run(totals, regions, stall_rounds=3) == 5

    def test_zero_overlap_scores_zero_with_min_contribution_zero(self):
        """N-16: log2(0) raised 'math domain error'."""
        from baitUtils.gap_filling_algorithm import OligoScorer, OligoMapping
        assert OligoScorer.score_oligo(OligoMapping('S', 'chrB', 5, 125), {'chrB': [(300, 310)]}, 0) == 0.0
        assert OligoScorer.score_oligo(OligoMapping('S', 'chrB', 5, 125), {'chrB': [(100, 110)]}, 0) > 0


@pytest.mark.skipif(not _has_bedtools(), reason="requires pybedtools and bedtools")
class TestFillForcedOligos:
    def test_forced_oligos_without_hit_are_kept_with_a_warning(self, dataset, tmp_path, run_cli, capsys):
        """N-15: forced baits without a passing hit vanished from the output."""
        forced = tmp_path / "forced.txt"
        forced.write_text("A_00\nNO_SUCH_BAIT\n")
        out = tmp_path / "sel.txt"
        run_cli(["fill", "--alignments", dataset["psl"], "--reference", dataset["reference"],
                 "--output", out, "--forced-oligos", forced])
        selected = out.read_text().split()
        assert {"A_00", "NO_SUCH_BAIT"} <= set(selected)
        assert "NO_SUCH_BAIT" in capsys.readouterr().err

    def test_fill_reaches_the_greedy_optimum_with_small_passes(self, dataset, tmp_path, run_cli, capsys):
        """N-14 end to end: 5 baits per pass used to stop at 3200 uncovered bases."""
        out = tmp_path / "sel.txt"
        run_cli(["fill", "--alignments", dataset["psl"], "--reference", dataset["reference"],
                 "--output", out, "--max-oligos-per-pass", "5", "--max-passes", "20"])
        assert "rounds without improvement" not in capsys.readouterr().err
        assert len(out.read_text().split()) > 15


# ---------------------------------------------------------- N-24, N-25

class TestReferenceRecords:
    def test_duplicate_ids_are_an_error(self, tmp_path):
        """N-24: two records named 'dup' were merged into one 1000 bp reference."""
        from baitUtils.cli_utils import InputError
        from baitUtils.mapping_utils import read_reference_fasta
        path = _fasta(tmp_path / "dup.fa", [("dup first", "A" * 3000), ("dup second", "C" * 1000)])
        with pytest.raises(InputError, match="Duplicate sequence id"):
            read_reference_fasta(path)

    def test_duplicate_ids_stop_evaluate_with_one_message(self, dataset, tmp_path, run_cli, capsys, fake_pblat):
        reference = _fasta(tmp_path / "dup.fa", [("chrA", "ACGT" * 100), ("chrA", "ACGT" * 50)])
        with pytest.raises(SystemExit) as excinfo:
            run_cli(["evaluate", "-i", dataset["baits"], "-r", reference, "-o", tmp_path / "out"])
        assert excinfo.value.code == 1
        err = capsys.readouterr().err
        assert err.count("Duplicate sequence id") == 1
        assert "Traceback" not in err

    def test_empty_record_is_skipped(self, dataset, tmp_path):
        """N-25: an empty record crashed CoverageAnalyzer (np.max of an empty array)."""
        from baitUtils.coverage_stats import CoverageAnalyzer
        reference = tmp_path / "ref.fa"
        reference.write_text(dataset["reference"].read_text().rstrip("\n") + "\n>chrEmpty\n")
        stats = CoverageAnalyzer(dataset["psl"], reference).analyze()
        assert set(stats["per_reference"]) == {"chrA", "chrB"}
        assert stats["reference_length"] == 5000

    @pytest.mark.skipif(not _has_bedtools(), reason="requires pybedtools and bedtools")
    def test_empty_record_does_not_stop_check(self, dataset, tmp_path, run_cli, capsys):
        reference = tmp_path / "ref.fa"
        reference.write_text(">chrEmpty\n" + dataset["reference"].read_text())
        run_cli(["check", "--alignments", dataset["psl"], "--reference", reference])
        assert "Skipped 1 empty record" in capsys.readouterr().err


# ------------------------------------------------------------ N-6 to N-8

class TestStatsFilters:
    def _stats(self, **overrides):
        stats = {'n_count': 0, 'gc_content': 50.0, 'melting_temperature': 70.0, 'masked_percentage': 0.0,
                 'hairpin_dg': -1.0, 'self_dimer_dg': -2.0}
        stats.update(overrides)
        return stats

    def test_na_fails_a_filter_the_user_set(self):
        """N-6: a sequence whose Tm or dG was NA passed the Tm and dG filters."""
        from baitUtils.sequence_statistics import SequenceFilter
        seq = "ACGT" * 30
        assert not SequenceFilter(min_tm=60).passes_filters(seq, self._stats(melting_temperature='NA'))
        assert not SequenceFilter(max_tm=90).passes_filters(seq, self._stats(melting_temperature='NA'))
        assert not SequenceFilter(min_hairpin_dg=-3).passes_filters(seq, self._stats(hairpin_dg='NA'))
        assert not SequenceFilter(min_dimer_dg=-6).passes_filters(seq, self._stats(self_dimer_dg='NA'))

    def test_na_passes_when_the_filter_is_not_set(self):
        from baitUtils.sequence_statistics import SequenceFilter
        stats = self._stats(melting_temperature='NA', hairpin_dg='NA', self_dimer_dg='NA')
        assert SequenceFilter(min_gc=40).passes_filters("ACGT" * 30, stats)

    def test_energy_filters_without_viennarna_are_an_error(self, dataset, tmp_path, run_cli, capsys, monkeypatch):
        """N-7: --min-hairpin-dg was silently ignored without ViennaRNA."""
        from baitUtils import sequence_analysis
        monkeypatch.setattr(sequence_analysis, "HAS_VIENNA_RNA", False)
        with pytest.raises(SystemExit) as excinfo:
            run_cli(["stats", "-i", dataset["baits"], "-o", tmp_path / "st", "-p", "1",
                     "--min-hairpin-dg", "-3", "-f"])
        assert excinfo.value.code == 1
        assert "need ViennaRNA" in capsys.readouterr().err


class TestPlotWithoutEnergies:
    def _table(self, tmp_path):
        """A stats table as written without ViennaRNA: the energy columns are empty."""
        import pandas as pd
        rng = np.random.default_rng(1)
        n = 40
        table = pd.DataFrame({
            'sequence_id': [f"b{i}" for i in range(n)],
            'length': 120,
            'gc_content': rng.uniform(35, 65, n),
            'melting_temperature': rng.uniform(70, 85, n),
            'entropy': rng.uniform(1.8, 2.0, n),
            'n_count': 0,
            'hairpin_dg': np.nan,
            'self_dimer_dg': np.nan,
            'kept': [i % 3 != 0 for i in range(n)],
        })
        path = tmp_path / "sequence_statistics.tsv"
        table.to_csv(path, sep="\t", index=False, na_rep="NA")
        return path

    def test_pca_and_coloured_histograms_skip_empty_columns(self, tmp_path, run_cli):
        """N-8: --color kept crashed and PCA dropped every row but reported success."""
        out = tmp_path / "plots"
        run_cli(["plot", "-i", self._table(tmp_path), "-o", out, "--plot-type", "pca", "histogram",
                 "--color", "kept"])
        names = {p.name for p in out.iterdir()}
        assert "histogram_gc_content.png" in names
        assert "histogram_hairpin_dg.png" not in names
        assert any(name.startswith("pca") for name in names)

    def test_pca_that_cannot_run_is_a_failure(self, tmp_path, run_cli):
        import pandas as pd
        path = tmp_path / "t.tsv"
        pd.DataFrame({'a': [1.0, np.nan, np.nan], 'b': [np.nan, 2.0, np.nan]}).to_csv(path, sep="\t", index=False)
        with pytest.raises(SystemExit) as excinfo:
            run_cli(["plot", "-i", path, "-o", tmp_path / "p", "--plot-type", "pca"])
        assert excinfo.value.code == 1


# ---------------------------------------------------------- N-38, N-39

class TestConfigValidation:
    def _run(self, run_cli, tmp_path, config, *argv):
        path = tmp_path / "cfg.json"
        path.write_text(config if isinstance(config, str) else json.dumps(config))
        with pytest.raises(SystemExit) as excinfo:
            run_cli(["--config", path, *argv])
        return excinfo.value.code

    @pytest.mark.parametrize("config, message", [
        ({"evaluate": {"strand": "reverse"}}, "invalid choice 'reverse'"),
        ({"evaluate": {"enable_html_report": "false"}}, "expected true or false"),
        ({"evaluate": {"plot_dpi": 0}}, "'evaluate.plot_dpi'"),
        ({"evaluate": {"gap_extend": -200}}, "'evaluate.gap_extend'"),
        ({"compare": {"sets": "A:x.fa"}}, "expected a list"),
        ({"common": {"threads": 2.5}}, "'common.threads'"),
        ({"evaluate": [1, 2]}, "must be a JSON object"),
        ("{bad json", "not valid JSON"),
    ])
    def test_invalid_values_are_usage_errors(self, run_cli, tmp_path, capsys, config, message):
        """N-38: config values skipped type and choices checks."""
        code = self._run(run_cli, tmp_path, config, "evaluate", "-i", "x.fa", "-r", "y.fa", "-o", tmp_path / "o")
        assert code == 2
        err = capsys.readouterr().err
        assert message in err
        assert err.count("'common.threads'") <= 1

    def test_missing_file_is_a_usage_error(self, run_cli, tmp_path, capsys):
        with pytest.raises(SystemExit) as excinfo:
            run_cli(["--config", tmp_path / "missing.json", "evaluate", "-i", "x", "-r", "y", "-o", "z"])
        assert excinfo.value.code == 2
        assert "cannot read configuration file" in capsys.readouterr().err

    def test_valid_values_are_converted(self):
        from baitUtils.__main__ import build_parser
        from baitUtils.config import apply_config
        parser, commands = build_parser()
        apply_config(commands, {"evaluate": {"strand": "minus", "min_identity": 95, "enable_html_report": False},
                                "compare": {"sets": ["A:a.fa", "B:b.fa"]}})
        args = parser.parse_args(["evaluate", "-i", "x", "-r", "y", "-o", "z"])
        assert (args.strand, args.min_identity, args.enable_html_report) == ("minus", 95.0, False)
        assert isinstance(args.min_identity, float)
        args = parser.parse_args(["compare", "-r", "y", "-o", "z"])
        assert args.sets == ["A:a.fa", "B:b.fa"]

    @pytest.mark.parametrize("flag", ["--conf", "--conf=", "--confi"])
    def test_abbreviated_config_option_is_refused(self, run_cli, tmp_path, flag):
        """N-39: --conf was accepted by the main parser and the file silently ignored."""
        path = tmp_path / "cfg.json"
        path.write_text('{"evaluate": {"strand": "minus"}}')
        argv = [flag + str(path)] if flag.endswith("=") else [flag, path]
        with pytest.raises(SystemExit) as excinfo:
            run_cli([*argv, "evaluate", "-i", "x", "-r", "y", "-o", "z"])
        assert "write --config in full" in str(excinfo.value.code)
