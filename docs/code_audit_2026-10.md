# baitUtils 0.5.1: code audit, gap analysis and switch analysis

Date: 2026-10-26. Audited revision: `22ab15a` (main, after PR #15).

## 1. Scope and method

- Every module under `baitUtils/` (31 files, 10 623 lines) was read in full,
  together with `tests/`, `README.md`, `CLAUDE.md`, `CONTRIBUTING.md`,
  `CHANGELOG.md`, `setup.py`, `requirements*.txt`, `pytest.ini`,
  `MANIFEST.in` and the two GitHub workflows.
- Every claim marked **verified** below was reproduced by running code
  (the CLI, a short script, or the test suite) on Python 3.13.17 with
  numpy 2.5.3, pandas 3.0.5, scipy 1.18.1, matplotlib 3.11.2,
  seaborn 0.13.2, plotly 7.1.0, biopython 1.88, ViennaRNA 2.7.2 and
  pybedtools 0.12. pblat, minimap2 and the bedtools binary were not
  available, so the `check`/`fill` end-to-end tests and the real mapper
  runs were covered by code reading only; the fake-pblat CLI tests ran.
  Claims marked **likely** rest on code reading alone.
- Tool runs: `pytest tests/` (245 collected: 240 passed, 5 skipped for
  missing bedtools/minimap2; one additional failure appears with
  plotly >= 6, see F-37), `ruff`, `flake8 --max-complexity`, `bandit`,
  `vulture`, plus a script that dumps every argparse option of every
  command and cross-checks it against the documentation (appendix A).

The audit has three parts: findings (section 2), the gap analysis
between documentation, intent and implementation (section 3), and the
switch analysis of the command-line interface (section 4). Section 5 is
a prioritised remediation plan. No code was changed in this audit.

## 2. Findings

Severity: **High** = wrong results or a crash on realistic input;
**Medium** = misleading output, silently ignored options, or robustness
holes; **Low** = cleanliness, dead code, small inconsistencies.

### 2.1 High

| # | Location | Finding | Status |
|---|---|---|---|
| F-1 | `quality_scorer.py:294` | `ZeroDivisionError` in `_score_reference_difficulty` when `expected_difficulty == 2.0` (one extreme-GC, challenging reference with no high-repeat fraction, e.g. a single AT-rich contig). `evaluate` and `compare` die after mapping with no report. For `expected_difficulty > 2` the denominator goes negative and the "better than expected" bonus is never applied. | verified (reproduced end to end) |
| F-2 | `coverage_viz.py:298-303` | `ZeroDivisionError` in `plot_coverage_heatmap` when every reference has zero depth, zero max depth and zero gaps, i.e. whenever no hit survives `--min-identity`/`--min-length`/`--strand`. The empty-PSL guard in `evaluate.py:370` only catches a zero-byte file. `evaluate` crashes instead of reporting "nothing mapped". | verified |
| F-3 | `sequence_analysis.py:120` | `calculate_self_dimer_dg("")` calls `RNA.fold_compound("&")`, which segfaults the interpreter; `except Exception` cannot catch it. An empty FASTA record (header followed by a blank line) kills `stats` when ViennaRNA is installed. | verified (faulthandler trace) |
| F-4 | `comparative_analyzer.py:297-302`, `comparative_report_generator.py:145,155,222,480,500`, `comparative_visualizations.py:459` | The quality score is 0 to 1 (`quality_scorer.py`, README, CHANGELOG 0.3.0), but the whole `compare` side still treats it as 0 to 10: the ranking adds `overall_score * 0.4` (range 0 to 0.4) to terms in the range 0 to 2.5, so the documented 40 % weight is about 6 % and the ranking is driven by breadth and mapping efficiency; the HTML report prints "0.8/10"; the radar plot divides by 10. `tests/test_comparative_analysis.py:30-62,581` assert the wrong scale. | verified (fixture: composite 5.92, quality contributes 0.34) |
| F-5 | `evaluate.py:395` vs `reference_analyzer.py:98` | `evaluate` writes `reference_analysis.get('summary')` but the analyzer returns `analysis_summary`, so `evaluation.json` always has `"reference_summary": null`. Nothing else consumes `coverage_correlations`, `breadth_by_gc`, `challenging_windows` or `window_features`: the Spearman correlations, GC-decile breadth and challenging windows that README (lines 200 to 206) and CHANGELOG 0.4.0 promise are computed and discarded. | verified |
| F-6 | `compare.py:201-207` vs `comparative_report_generator.py:51,59` | `--plot-format` and `--plot-dpi` are honoured by the first `ComparativeVisualizer`, then `ComparativeReportGenerator` builds a second visualizer with defaults and regenerates every plot into the same directory. With `--plot-format pdf` both PDF and PNG files appear and the report embeds the PNG; with the default format the 300 dpi pass overwrites the user's dpi. Every plot and every statistical test is computed three times. | verified |
| F-7 | `interactive_plots.py:428-445` | The "Quality Trends" panel of `quality_assessment.html` plots a fabricated series (`0.4` as "initial", `overall + 0.2` as "potential"); the comment says "would be actual historical data in a real implementation". The neighbouring "Theoretical vs Actual" panel and the threshold lines read keys (`theoretical_comparison`, `quality_thresholds`) that `QualityScore.to_dict()` never emits, so they are empty. The 0.3.0 changelog says simulated values were removed; this one remains. | verified in generated HTML |
| F-8 | `__main__.py:88-91` | `--config` is detected by scanning `sys.argv` for the literal `--config`: `--config=file.json` is silently ignored (argparse still accepts it, so no error), and a trailing `--config` raises an uncaught `IndexError`. | verified |
| F-9 | `config.py:40-53` and README config example | A `common` block is applied to every command, so the README's own example (`"common": {"threads": 4, "log_level": "WARNING"}`) emits seven warnings and `log_level` only takes effect for `check` and `fill`; `stats`, `plot`, `map`, `evaluate` and `compare` have no `log_level` destination (`stats` has `processes`, not `threads`). A misspelled command key (`"evaluat"`) is silently ignored. | verified |

### 2.2 Medium

| # | Location | Finding | Status |
|---|---|---|---|
| F-10 | `report_generator.py:119-131` | The executive summary recomputes the quality category with thresholds 0.8/0.6/0.4 instead of reading `quality_scores['category']` (scorer: 0.85/0.70/0.50). On the test fixture `evaluation.json` says 0.843 "Good" while the HTML report says "Excellent". | verified |
| F-11 | `comparative_report_generator.py:140-514`, `report_generator.py:378-386` | Set names from `--sets "Name:path"` and reference ids from FASTA headers are interpolated into HTML without `html.escape`; `--sets "A<b>x</b>:..."` produces raw markup eight times. Reports are meant to be shared. | verified |
| F-12 | `differential_analysis.py:126-147,199-213,257-260,294-301` | `interpretation` is computed on the raw p-value before correction, while stars and CSS class use the adjusted one; with Bonferroni over three metrics and raw p = 0.031 the block reads "Difference ... at alpha 0.05" with stars "ns". | verified |
| F-13 | `differential_analysis.py:91,304-310`, `comparative_report_generator.py:632` | Star labels use fixed 0.001/0.01/0.05 cut-offs regardless of `--significance-level`; the legend always says "p<=0.05". | verified |
| F-14 | `differential_analysis.py:178` vs `coverage_stats.py:355` vs `gap_analysis.py:106-130` | The paired "gap count" test uses `per_reference[ref]['gaps']` (every run below `min_coverage`, no size floor) while `total_gaps`, the gap-size test and the comparison matrix use `GapAnalyzer` gaps at or above `min_gap_size` (100). Fixture: per-reference 2 + 1 versus total 1. | verified |
| F-15 | `comparative_analyzer.py:356,388` | Gap overlap uses closed-interval logic (`end1 < start2`) on half-open gaps, so touching gaps count as overlapping with length 0 and are excluded from "unique gaps". | verified |
| F-16 | `comparative_report_generator.py:582` | Interactive plot links use the basename only; the files live in `interactive_comparative_plots/`, so every link is broken. | verified |
| F-17 | `report_generator.py:538-562` | The "Data Exports" section links to `data/coverage_per_position.csv`, `data/gap_regions.bed` and `data/statistics_summary.json`, but `evaluate.py` never calls `CoverageAnalyzer.export_coverage_data` or `export_data_files`; the `coverage.bedgraph` export added in 0.5.1 is unreachable from the CLI. | verified |
| F-18 | `evaluate.py:251-306` | With `--no-html-report --no-interactive-plots --no-benchmarking` the reference analysis and the quality score are skipped as well (`quality_score: null`, nothing printed), although README describes them as part of `evaluate` and the flags as skipping only those outputs. | verified by reading |
| F-19 | `quality_scorer.py:63-70,190-213` | The documented `depth_uniformity 0.8 (1 - Gini)` target is only used by `BenchmarkAnalyzer`; the scorer scores uniformity from `coverage_cv` with hard-coded 0.5/2.0 cut-offs. The docstring claim that targets can be overridden is false for this component, and `tests/test_quality_scorer.py:155-176` (varies Gini) is vacuous. | verified |
| F-20 | `quality_scorer.py:55-61` vs `279-299` | Docstring says "credit for hard references"; the implementation subtracts up to 0.8 for hard references and only adds a bonus when breadth exceeds an expectation, so at equal breadth a harder reference always scores lower. | verified |
| F-21 | `quality_scorer.py:120-155` | Weight overrides are not normalised: `weights={'coverage_breadth': 1.3}` gives an overall score of 1.56, above every category threshold. | verified |
| F-22 | `report_generator.py:484`, `interactive_plots.py:251-314` | The report's reference heatmap reads `complexity` (analyzer emits `shannon_entropy`), so that row is always 0; the gap-feature panels read `gap_analysis['gap_features']`, which `GapAnalyzer` keeps in an attribute but never puts in its result, so both panels are always empty. | verified |
| F-23 | `gap_analysis.py:486-656` | `_correlate_with_reference_features`, `_calculate_gap_feature_correlation`, `_identify_feature_driven_gaps` and `_calculate_problematic_overlap` (about 170 lines) only run when `analyze(reference_analysis=...)` is given; no caller passes it, so `reference_correlations` is always `{}`. The "correlation strength" inside is the fixed-formula heuristic the 0.4.0 changelog says was removed. | verified by grep |
| F-24 | `coverage_analysis.py:466-472`, `gap_filling.py:323-335` | `--uncovered_fasta` without `--n_split_fasta` skips the uncovered FASTA too; the warning says only the N-split output is skipped. | verified by reading |
| F-25 | `sequence_mapping.py:106-107` | The filtered alignment file is always named `<prefix>-mapping_filtered.psl`, also for `--mapper minimap2` where it contains PAF rows; README says `run-mapping_filtered.*`. `evaluate --keep-intermediates` likewise copies a PAF to `intermediates/mapping.psl` (`evaluate.py:318`). | verified by reading |
| F-26 | `sequence_mapping.py`, `mapping_utils.py:282-295` | pblat-only switches (`--minMatch`, `--minScore`, `--minIdentity`) are silently ignored under `--mapper minimap2`, and `--minimap2-preset` under pblat; `validate_identity_parameters` still enforces `--filterIdentity >= --minIdentity` even though `--minIdentity` has no effect with minimap2. minimap2 is run with a hard-coded `-N10` secondary limit. | verified by reading |
| F-27 | `sequence_statistics.py:378-380` | A missing `stats` input logs an error and returns, exit code 0. Every other command exits 1. | verified |
| F-28 | `coverage_stats.py:313-317`, `differential_analysis.py:198`, `json_export.py:38-41` | `coverage_cv` is `inf` when nothing is covered and the paired effect size is `inf` when all differences are equal; `to_jsonable` maps NaN to `null` but not `inf`, so `evaluation.json` / `comparison.json` contain the non-standard token `Infinity`, which strict JSON parsers reject. | verified |
| F-29 | `coverage_stats.py:349-350` | `plus_hits`/`minus_hits` per reference are computed by scanning all mappings once per reference: O(references x hits). A fragmented assembly (10^4 contigs) with 10^5 hits costs 10^9 comparisons. | verified by reading |
| F-30 | `coverage_stats.py`, `gap_analysis.py:96`, `reference_analyzer.py:103` | `evaluate` parses the reference FASTA three times and holds it as a Python string plus an int32 depth array plus a boolean N mask (about 6 bytes per base, so a 1 Gb reference needs about 6 GB plus pandas frames). Acceptable for targets, not for whole genomes; README says only "linear in reference length". | verified by reading |
| F-31 | `gap_filling_algorithm.py:152-174` | The spacing constraint rejects any candidate whose start is within `--spacing_distance` of an already selected start on the same reference, so `fill --min_coverage N` with N > 1 can never stack baits on one locus; depth above 1 is reachable only from distinct baits offset by at least the spacing. Not documented. | verified by reading |
| F-32 | `gap_filling_algorithm.py:305-317` | Because the base-count check breaks the loop whenever bases do not decrease, `--stall_rounds` can only count passes where the region count stalls while the base count still drops; a pass whose selections did not help is kept in the final selection. | verified by reading |
| F-33 | `coverage_evaluation.py:92` vs `gap_filling.py:353` vs `evaluate.py:98` | `--min_coverage` defaults to 10 in `check` but 1 in `fill`, `evaluate` and `compare`; `check --min_similarity` defaults to 95 while `map --filterIdentity` and `evaluate --min-identity` default to 90; `check`/`fill --min_length` default 100 while `map --min-length` defaults 0. The README workflow (`map` then `fill`) therefore re-filters the alignments with different thresholds without saying so. | verified |
| F-34 | `pytest.ini:1` | The section header is `[tool:pytest]` (the `setup.cfg` form); in `pytest.ini` only `[pytest]` is read, so `addopts`, `markers` and `filterwarnings` are all ignored. The `--strict-warnings` in `addopts` is not a pytest option and would abort every run if the file were read. | verified |
| F-35 | `CLAUDE.md:78-89`, `CONTRIBUTING.md:36` | Both recommend `python -m unittest discover tests/`, which runs 134 of the 245 tests (every pytest-fixture test, including all end-to-end CLI tests, is skipped silently). | verified |
| F-36 | `CLAUDE.md:212,233` | Usage examples use `--enable-html-report` and `--enable-statistical-analysis`, both removed in 0.3.0; the commands exit with "unrecognized arguments". | verified |
| F-37 | `tests/test_phase2_integration.py:220-239` | `test_report_generation` patches `builtins.open` globally; with plotly >= 6 (validators loaded lazily from JSON) it fails. `setup.py` allows `plotly>=5.0.0`, so CI will fail on a fresh resolver. | verified (fails with plotly 7.1.0) |
| F-38 | `coverage_viz.py:56`, `setup.py` | `plt.style.use('seaborn-v0_8')` needs matplotlib >= 3.6; the pin is `>=3.4.0`. `plotting_utils.py:58,76,111,129` pass `palette` without `hue`, which seaborn 0.13 deprecates and 0.14 removes; the pin is `>=0.11`. | likely / verified warning |
| F-39 | `plotting_utils.py:207-212` | `plot --plot_type pca` passes the frame straight to `StandardScaler`/`PCA`; a single NaN aborts the command after earlier plots were written. | verified |
| F-40 | `comparative_report_generator.py:600-631`, `report_generator.py:564-587` | Methodology sections describe "pblat" only, "t-test/ANOVA", "Cohen's d, eta-squared" and A to F grades; the code uses Wilcoxon/Friedman/Mann-Whitney with rank-biserial r and Kendall's W, Excellent/Good/Fair/Poor categories, and supports minimap2. The CSS has `.grade-a..f` while the HTML emits `grade-excellent/good/fair/poor`, so grade colouring never applies. | verified |

### 2.3 Low

| # | Location | Finding |
|---|---|---|
| F-41 | `mapping_utils.py:461` | `PSLParser.last_hit_table` is a class attribute: global mutable state written by every `parse_psl_file` call. |
| F-42 | `mapping_utils.py:377,460` vs `coverage_analysis.py:69,317` | Two unrelated classes named `SequenceLoader` and two named `PSLParser` in different modules. |
| F-43 | `evaluate.py:476`, `quality_scorer.py:312`, `gap_analysis.py:415`, `comparative_report_generator.py:462` | Four independent recommendation generators with different thresholds (breadth < 80 in one, gap count > 100 vs > 50 in two others). |
| F-44 | `evaluate.py:469` | `benchmark_analysis.txt` is rendered through `BenchmarkAnalyzer({}, {}, {}, None)`, so the "Overall quality score" header is always omitted. |
| F-45 | `quality_scorer.py:263-266` | Without reference analysis the difficulty component silently returns 0.7, or 1.0 when `analysis_summary` is missing; undocumented. |
| F-46 | `quality_scorer.py:390-408`, `interactive_plots.py:388-409` | Benchmark ratios are unclipped (gap ratio 20 with zero gaps) and plotted against a "target = 1" ring. |
| F-47 | `benchmark.py:86` | Missing `coverage_gini` defaults to 0.5, producing a fabricated uniformity line and recommendation. |
| F-48 | `sequence_features.py:188-192` | Window `repeat_density` numerator can exceed the denominator (101.1 % for 1000 bp windows, k = 12) before clipping. |
| F-49 | `reference_analyzer.py:50,68` | `coverage_data` parameter is stored and never read. |
| F-50 | `reference_analyzer.py:291` vs `193` | Window homopolymer flag uses a hard-coded 10 % while `LONG_HOMOPOLYMER = 10` is a run length; both labelled "homopolymer" in reasons. |
| F-51 | `reference_analyzer.py:341` | `coverage_correlations['applicable']` is True even when every rho is NaN. |
| F-52 | `sequence_analysis.py:385-409` vs `sequence_features.py:53-60` | `stats` entropy counts N as a symbol; the reference analyzer ignores N. The two "entropy" columns are not comparable and the docs do not say so. |
| F-53 | `sequence_analysis.py:253-260` | Only `ImportError` is caught around the ViennaRNA DNA parameter loader; an old ViennaRNA raises `AttributeError` at import and breaks every command. |
| F-54 | `gap_analysis.py:294-299` | `gap_percentage` falls back to `reference_length = 1`, giving percentages in the thousands when the coverage dict is incomplete. |
| F-55 | `gap_analysis.py:306-312` | Size bins start at `min_gap_size`; with `--min-gap-size 1000` the first bin is "1000-500". |
| F-56 | `gap_analysis.py:243-263` | `_calculate_low_complexity` recomputes a `Counter` for every 20 bp window in pure Python: about one second per 100 kb of gap sequence. |
| F-57 | `coverage_viz.py:233-234` | Gap labels longer than 20 characters drop the coordinates entirely. |
| F-58 | `comparative_visualizations.py:86,162` | Ranking bar colours are applied positionally after re-sorting, so a set can get another set's colour; `Set3` has 12 colours, so more than 12 sets fall back to defaults. |
| F-59 | `comparative_visualizations.py:339-341,457-460` | 0/0 normalisation gives NaN heatmap cells when all sets have zero gaps; `Mean_Depth_x / 20` is an arbitrary constant. |
| F-60 | `comparative_visualizations.py:61-62` | Sets global `plt.rcParams` as a side effect, twice. |
| F-61 | `compare.py:283,327` | Set specs split on the first `:` (Windows drive letters become the set name); `print_summary` takes an unused `comparison_matrix`. |
| F-62 | `differential_analysis.py:109` | `np.concatenate(means)` raises on an all-empty reference. |
| F-63 | `statistical_plots.py:135`, `plotting_utils.py:393-394` | Help text says ".txt or .txt.gz" (stats writes `.tsv`); the "remove kept from numeric columns" branch is unreachable because `kept` is bool. |
| F-64 | `coverage_evaluation.py:109` | `check --longest_uncovered_out` still defaults to `longest_uncovered.txt` in the working directory; 0.3.0 removed the same default from `fill` only. `fill --output` and `map -o` also default to the working directory. |
| F-65 | `gap_filling.py:215,231` | `args.temp_dir` is overwritten with the run's work directory; `if args.reference:` is dead because `--reference` is required. |
| F-66 | `coverage_analysis.py:100-122` | References present in the alignment file but absent from the reference FASTA are missing from the genome file, so bedtools fails with an opaque error instead of a clear message. |
| F-67 | `sequence_statistics.py:359-361` | `-p` defaults to `cpu_count() - 1`, so a one-sequence run forks a process pool; `--sample` takes the first N records, not a sample. |
| F-68 | `json_export.py:52`, `mapping_utils.py:267,366,370` | `subprocess.run` without `check=` (bandit B603, all low; arguments are lists, no shell). |
| F-69 | `setup.py`, `requirements.txt` | `tqdm` is in `install_requires` but never imported; `pybedtools` is in `requirements.txt` but only an extra in `setup.py`; `MANIFEST.in` includes `docs/` and `examples/`, which do not exist; no `pyproject.toml` (the `setup.py`-only editable install fails on Debian-patched setuptools with Python 3.13: `AttributeError: install_layout`); classifiers stop at 3.12 while the suite passes on 3.13. |
| F-70 | `.github/workflows/ci.yml` | `actions/setup-python@v4`, `actions/cache@v3`, `codecov/codecov-action@v3` are superseded (Node 16 runtime); the integration job only runs `stats` and `plot`; the `security-check` and lint-report steps never fail the build. |
| F-71 | `sequence_statistics.py:309`, others | Each subcommand's `--version` reports a different program name (`sequence-stats`, `coverage-evaluation`, `gap-filling`, `sequence-mapping`, `statistical-plots`) that does not match the command name. |
| F-72 | `mapping_utils.py:274,331`, `sequence_features.py:22`, `sequence_analysis.py:347,365`, `report_generator.py:627`, `mapping_utils.py:418` | Dead symbols: `MAPPERS`, `STRAND_CHOICES`, `BASE_CODES`, `clean_sequence`, `validate_dna_sequence`, `export_data_files`, `SequenceLoader.get_sequence_ids`, `PSLHit.target_span`, `HAS_PYBEDTOOLS` (set, never read). |
| F-73 | 32 `except Exception` blocks, 223 root-logger calls, 8 functions over the complexity 10 limit | Mostly log-and-reraise; the broad excepts in `comparative_visualizations.py:501-534` and `interactive_plots.py:652-667` swallow plotting failures so a report can be generated with figures missing and only a log line. |

### 2.4 Things checked and found correct

- PSL and PAF parsing, BLAT identity (`milli_bad`), aligned-block coverage, uncovered-run merging across reference boundaries, the fake-pblat end-to-end numbers in `tests/conftest.py`.
- Wilcoxon/Friedman/Mann-Whitney argument shapes, Kendall's W, rank-biserial r, Holm and Benjamini-Hochberg implementations, per-family correction as documented; "not applicable" handling for identical sets, a single reference, zero gaps and NaN p-values.
- All producer keys read by the comparative modules exist; `plot` supports exactly the 11 plot types README lists and honours `--format` and `--color` for each.
- Tm calculation, DNA parameters at `--hyb-temp`, strand filtering, `--orient-to-reference`, config-file precedence (command line wins, required options become optional), the global `--log-level` versus per-command `-v/-l/--quiet` precedence.
- No `np.random`, no leftover TODO/FIXME markers, ruff's CI selection and the undefined-name checks pass, no shell=True subprocess calls, pandas 3 / numpy 2.5 compatible (no `applymap`, `append`, `is_categorical_dtype`).

## 3. Gap analysis

### 3.1 Documented but not delivered

| Promise (where) | Reality | Finding |
|---|---|---|
| Spearman feature/coverage correlations, breadth by GC decile, challenging windows "in the report" (README 200-206, CHANGELOG 0.4.0) | Computed, never written: `reference_summary` is always null, nothing consumes the correlation tables | F-5 |
| Quality score "weighted sum ... weights and targets documented in QualityScorer and can be overridden" (README 206-210) | Uniformity target unused; overrides unnormalised; difficulty component penalises rather than credits | F-19, F-20, F-21 |
| Quality categories Excellent 0.85 / Good 0.70 / Fair 0.50 (README) | HTML report uses 0.8/0.6/0.4; compare report prints "x/10" | F-10, F-4 |
| `compare` ranking "composite scoring" with a 40 % quality weight (CLAUDE.md) | Effective weight about 6 % | F-4 |
| `--plot-format`, `--plot-dpi` on compare | Overwritten by a second default-settings pass | F-6 |
| `--no-html-report --no-interactive-plots --no-benchmarking` "skip those steps" | Also skips reference analysis and the quality score | F-18 |
| `coverage.bedgraph`, per-position CSV, gap BED exports (CHANGELOG 0.5.1; report "Data Exports" links) | Unreachable from the CLI; report links are dead | F-17 |
| `run-mapping_filtered.*` for minimap2 (README 107) | Always `.psl` | F-25 |
| `--config` "supplies option defaults ... `common` block" with the README example | Example emits seven warnings; `log_level` ignored by five commands; `--config=` form ignored | F-8, F-9 |
| "simulated values removed" (CHANGELOG 0.3.0) | Fabricated quality-trend series remains in the interactive report; heuristic correlation remains in dead gap-analysis code | F-7, F-23 |
| `python -m unittest discover tests/` (CLAUDE.md, CONTRIBUTING.md) | Runs 55 % of the suite | F-35 |
| `--enable-html-report`, `--enable-statistical-analysis` (CLAUDE.md) | Removed in 0.3.0 | F-36 |
| pytest configuration (`pytest.ini`) | Never read | F-34 |
| Methodology sections of both HTML reports | Describe tests, grades and a mapper the code does not use | F-40 |
| Help text `plot -i` ".txt or .txt.gz" | Input is the `.tsv` that `stats` writes | F-63 |

### 3.2 Capability gaps (not promised, but expected from the design)

| Gap | Why it matters |
|---|---|
| No off-target control in `evaluate`/`compare`: every passing hit of a multi-mapping bait adds depth at every locus, and there is no `--max-hits` as in `map`. The per-bait hit table is built (`coverage_stats.py:149`) but never written. | Depth, breadth and the quality score are inflated by repetitive baits; the number is not inspectable. |
| `compare` cannot set `--min-gap-size`, `--gap-extend`, `--reference-analysis-window`, `--offline-plots`, nor skip the HTML/benchmark steps; the values are hard-coded in `comparative_analyzer.py:142-172`. | Gap counts in `comparison.json` are not comparable with an `evaluate` run that used other values. |
| `evaluate`/`compare` have no `--temp-dir` (check/fill do) and always map into the system temp directory. | Large mappings on a small `/tmp` fail without recourse. |
| `fill` cannot build depth above 1 at a locus because of the spacing rule (F-31); `--min_coverage` above 1 is accepted without a warning. | Users asking for 2x fill get a selection that cannot satisfy it. |
| Only `stats` and `plot` read gzipped input; `map`, `evaluate`, `compare`, `check`, `fill` do not. | Inconsistent for the same FASTA. |
| No empty-input handling: an empty bait record segfaults `stats` (F-3), zero surviving hits crash `evaluate` (F-2). | Pipelines on filtered subsets hit these routinely. |
| Memory: the whole reference is held three ways (F-30); no streaming or per-reference processing. | Whole-genome references are out of reach on ordinary machines. |
| No machine-readable output for `check` and `fill` beyond TSV; no JSON with versions and arguments as `evaluate` has. | Provenance differs between the two halves of the toolkit. |
| No `--seed`/determinism concerns: fine. No `--dry-run`/`--overwrite` protection: every command overwrites existing output silently. | Reruns clobber results. |

### 3.3 Test coverage gaps

- Modules with no test importing them: `sequence_mapping`, `coverage_evaluation`, `gap_filling`, `gap_filling_algorithm`, `statistical_plots`, `bedtools_support` (the `map`, `check`, `fill` and `plot` commands are exercised only through the end-to-end fixtures, and `check`/`fill` only when bedtools is installed, which the default CI matrix lacks).
- The only `compare` end-to-end runs compare a set with itself, so every statistical test is "not applicable" and the p-value, star, effect-size and Friedman paths of the report are never rendered.
- `tests/test_comparative_analysis.py:428-483` and `tests/test_phase2_integration.py:298-368` wrap plotting in `try/except` or mock `plotly.graph_objects.Figure`/`seaborn.heatmap`, so a failing plot still passes and key mismatches (F-22) are invisible. Fixture dicts there have shapes the producers never emit.
- `tests/test_quality_scorer.py` asserts `score < 5.0` and `0 <= overall <= 10` on 0-to-1 values; the Gini test cannot fail (F-19); custom weights are assigned after construction rather than through the constructor.
- Not covered anywhere: HTML escaping, `--plot-format`/`--plot-dpi`, gap overlap/unique-gap logic, zero-hit mappings through the visualizer, empty sequences through `SequenceAnalyzer`, `--config=` form, `fill --min_coverage > 1`, `map --mapper minimap2` without the binary (a fake minimap2 like the fake pblat would do), formats other than PNG.
- The CI `integration-test` job runs only `stats` and `plot`; `external-tools` is the only job that runs `check`/`fill`.

## 4. Switch analysis

Appendix A lists all 118 switches with type, default and documentation
status (generated from the parsers, not from the docs).

### 4.1 Naming conventions

Four styles coexist, often inside one command:

| Command | word (`--strand`) | kebab (`--min-length`) | snake (`--min_length`) | camel (`--minIdentity`) |
|---|---|---|---|---|
| stats | 16 | 3 | 0 | 1 (`--noNs`) |
| plot | 6 | 0 | 1 (`--plot_type`) | 0 |
| map | 7 | 5 | 0 | 5 |
| check | 3 | 1 (`--log-level`) | 12 | 0 |
| fill | 4 | 1 | 17 | 0 |
| evaluate | 8 | 15 | 0 | 0 |
| compare | 8 | 10 | 0 | 0 |

`check`/`fill` are snake_case except for `--log-level`; `map` mixes
pblat's camelCase (`--minIdentity`, `--filterIdentity`, `--minMatchCount`)
with kebab (`--min-length`, `--max-hits`); `stats` is single-word with
one camelCase (`--noNs`) and three kebab additions. The config file
uses destination names, so the user must know that `--min-length` is
`min_length` but `--minIdentity` is `minIdentity`.

### 4.2 Same concept, different switch, different default

| Concept | stats | plot | map | check | fill | evaluate | compare |
|---|---|---|---|---|---|---|---|
| bait FASTA | `-i/--input` | (`-i` = stats TSV) | `-i/--input` | (via alignments) | (via alignments) | `-i/--input` | `--sets Name:path` |
| reference FASTA | | | `-q/--query` | `--reference` | `--reference` | `-r/--reference` | `-r/--reference` |
| output location | `-o/--outdir` (required) | `-o/--outdir` (`plots`) | `-o/--outdir` (`.`) + `--prefix` | per-file `--*_out` | `--output` (`selected_oligos.txt` in cwd) | `-o/--output` (required) | `-o/--output` (required) |
| threads | `-p/--processes` (cpu-1) | | `-X/--threads` (1) | | | `--threads` (1) | `--threads` (1) |
| identity floor | | | `--minIdentity` 90 (pblat) + `--filterIdentity` 90 (int) | `--min_similarity` 95.0 | `--min_similarity` 95.0 | `--min-identity` 90.0 | `--min-identity` 90.0 |
| aligned-length floor | | | `--min-length` 0 | `--min_length` 100 | `--min_length` 100 | `--min-length` 100 | `--min-length` 100 |
| depth threshold | | | | `--min_coverage` 10.0 | `--min_coverage` 1.0 | `--min-coverage` 1.0 | `--min-coverage` 1.0 |
| alignment input | | | | `--alignments`/`--psl` | `--alignments`/`--psl` | (internal) | (internal) |
| temp directory | | | | `--temp_dir` | `--temp_dir` | none | none |
| verbosity | `-v/--verbose` | `-l/--log` | `-l/--log` | `--log-level` | `--log-level` | `-l/--log`, `--quiet` | `-l/--log`, `--quiet` |
| strand | | | none (hits.tsv records it) | `--strand` | `--strand` | `--strand` | `--strand` |
| mapper | | | `--mapper`, `--minimap2-preset` | | | same | same |

Short flags are reused with unrelated meanings: `-q` is the minimum Tm
in `stats` and the reference FASTA in `map`; `-l` is `--log` in four
commands but `check`/`fill` use `--log-level`; `-X` for threads exists
only in `map`; `-n`/`-x`/`-z`/`-K`/`-L`/`-N`/`-c`/`-s` are `stats`-only
single letters with no mnemonic link to their long names.

### 4.3 Switches that are accepted but have no effect

| Switch | Condition | Finding |
|---|---|---|
| `map --minMatch`, `--minScore`, `--minIdentity` | `--mapper minimap2` | F-26 |
| `map --minimap2-preset`, `evaluate/compare --minimap2-preset` | `--mapper pblat` | F-26 |
| `compare --plot-format`, `--plot-dpi` | always (second pass overwrites) | F-6 |
| `compare --significance-level` | star labels and legend | F-13 |
| `--config=file.json` (equals form) | always | F-8 |
| `common.log_level` in a config file | `stats`, `plot`, `map`, `evaluate`, `compare` | F-9 |
| `fill --min_coverage N` with N > 1 | cannot be satisfied by stacking | F-31 |
| `fill --stall_rounds` | mostly pre-empted by the base-count break | F-32 |
| `check --uncovered_fasta` | without `--n_split_fasta` | F-24 |
| `evaluate --gap-extend` | only widens the sequence window used for gap features; gap coordinates and sizes are unchanged (undocumented) | |
| `evaluate --keep-intermediates` | keeps a mislabelled `.psl` under minimap2 | F-25 |
| `stats --length` | only with `--complete` (documented) | |

### 4.4 Switches missing from the documentation

Undocumented in README: `stats -p/--processes`, `-s/--sample`;
`map --minMatch`, `--minScore`, `--minMatchCount`; `check`/`fill
--uncovered_length_cutoff`, `--temp_dir`; `fill --max_oligos_per_pass`;
`evaluate --gap-extend`, `--plot-dpi`, `--keep-intermediates`;
`compare --plot-dpi`. Documented but removed: CLAUDE.md's
`--enable-html-report`, `--enable-statistical-analysis` (F-36).

### 4.5 Switches that should exist

- `evaluate`/`compare`: `--max-hits` (off-target filter, as in `map`),
  `--temp-dir`, and for `compare` the `evaluate` analysis switches
  (`--min-gap-size`, `--reference-analysis-window`, `--offline-plots`,
  `--no-*`).
- `map`: `--strand` (restrict the hit table and the mapped set),
  `--max-secondary` for minimap2 (currently hard-coded 10).
- `check`/`fill`: a `--json` provenance output; `--quiet`/`--verbose`
  parity with the other commands.
- `fill`: a warning or `--allow-stacking` when `--min_coverage > 1`.
- All commands: a consistent `--force`/overwrite policy for existing
  output directories.

### 4.6 Validation gaps

`evaluate`/`compare` check identity and coverage ranges but not
`--min-length >= 0`, `--threads >= 1`, `--plot-dpi > 0`,
`--reference-analysis-window > 0`, or `--min-gap-size >= 1`; `check`/`fill`
validate nothing (a negative `--extend_region` silently shrinks regions;
`--spacing_distance 0` disables the constraint). `map`'s
`percent` type rejects `90.5` although pblat accepts fractional identity.

## 5. Recommended remediation plan

Phase 1, correctness (one release, breaking for nobody):

1. Guard the two divisions (F-1, F-2) and the empty-sequence ViennaRNA
   calls (F-3); add the three regression tests.
2. Fix the 0-to-1 scale in the `compare` ranking, report and radar
   (F-4) and the HTML category thresholds (F-10); fix the test fixtures.
3. Write `analysis_summary`, `coverage_correlations` and
   `challenging_windows` to `evaluation.json` and add the correlation
   and GC-decile tables to the HTML report (F-5); emit `gap_features`
   and use `shannon_entropy` (F-22); drop the fabricated trend panel and
   the dead `theoretical_comparison` branches (F-7).
4. Pass the configured visualizer into `ComparativeReportGenerator`
   (F-6); compute interpretation after correction and derive stars from
   the chosen alpha (F-12, F-13); use one gap definition (F-14); `<=`
   for half-open overlap (F-15); relative links (F-16).
5. `html.escape` every name and id in both report generators (F-11).
6. Map `inf` to `null` in `to_jsonable` (F-28); exit 1 on a missing
   `stats` input (F-27); fix the `--uncovered_fasta` gating (F-24).

Phase 2, CLI and configuration (minor release, with aliases):

7. Parse `--config` through argparse (`parse_known_args` on a parent
   parser) so `--config=` works and a trailing `--config` errors
   cleanly; warn once per unknown command key; only warn about `common`
   keys that no command uses; add `log_level` (or map `--log-level`) to
   every command (F-8, F-9).
8. Adopt one naming style (kebab-case, matching `evaluate`), keep the
   old spellings as hidden aliases for one release, and align names and
   defaults for the shared concepts in section 4.2 (`--reference`,
   `--threads`, `--min-identity`, `--min-length`, `--min-coverage`,
   `--temp-dir`, `--verbose`/`--quiet`). Decide the `check`
   `--min-coverage` default (1, like everything else) and
   `--min-similarity` (90, like `map`).
9. Warn when a mapper-specific switch is given for the other mapper;
   drop the `filterIdentity >= minIdentity` check under minimap2; name
   the filtered file by format (F-25, F-26).
10. Add `--max-hits` and `--temp-dir` to `evaluate`/`compare`, and the
    hard-coded analysis parameters to `compare`; write the hit table
    from `evaluate`; call `export_coverage_data` or remove the dead
    links (F-17); compute the quality score unconditionally (F-18).

Phase 3, scoring and docs:

11. Use the uniformity target, normalise weight overrides, and make the
    difficulty component do what its docstring says or rename it
    (F-19 to F-21); consolidate the four recommendation generators.
12. Regenerate both methodology sections from the code's own constants
    (F-40); fix CLAUDE.md and CONTRIBUTING.md (F-35, F-36); fix
    `pytest.ini` to `[pytest]` and remove `--strict-warnings` (F-34).
13. Packaging: add `pyproject.toml`, drop `tqdm`, align
    `requirements.txt` with `setup.py`, bump the matplotlib/seaborn/plotly
    pins or fix the deprecations, update the GitHub actions, add 3.13
    to the matrix (F-37, F-38, F-69, F-70).

Phase 4, scale and tests:

14. Compute per-reference strand counts in one pass (F-29); consider a
    single shared reference loader and per-reference processing (F-30).
15. Add a fake minimap2 to the fixtures, a `compare` fixture with two
    different sets (so the statistical paths render), bedtools to the
    default CI job, and tests for the items in section 3.3.

## Appendix A: switch inventory

Generated from the argparse definitions. "In README" means the long
form appears in README.md. Global switches in front of the command:
`--version`, `--config FILE`, `--log-level LEVEL`.

| Command | Switch | Type | Default | Style | In README |
|---|---|---|---|---|---|
| stats | `-i --input` | str |  | word | no |
| stats | `-o --outdir` | str |  | word | no |
| stats | `-L --length` | int | 120 | word | yes |
| stats | `-c --complete` | flag | False | word | yes |
| stats | `-N --noNs` | flag | False | camel | yes |
| stats | `-n --mingc` | float |  | word | yes |
| stats | `-x --maxgc` | float |  | word | yes |
| stats | `-q --mint` | float |  | word | yes |
| stats | `-z --maxt` | float |  | word | yes |
| stats | `--na` | float | 50.0 | word | yes |
| stats | `--dnac1` | float | 250.0 | word | yes |
| stats | `--dnac2` | float | 250.0 | word | yes |
| stats | `--hyb-temp` | float | 65.0 | kebab | yes |
| stats | `--min-hairpin-dg` | float |  | kebab | yes |
| stats | `--min-dimer-dg` | float |  | kebab | yes |
| stats | `-K --maxmask` | float |  | word | yes |
| stats | `-p --processes` | int | 3 | word | no |
| stats | `-s --sample` | int |  | word | no |
| stats | `-f --filter` | flag | False | word | yes |
| stats | `-v --verbose` | flag | False | word | yes |
| plot | `-i --input` | str |  | word | no |
| plot | `-o --outdir` | str | plots | word | no |
| plot | `--columns` | str |  | word | yes |
| plot | `--plot_type` | str | ['histogram'] | snake | yes |
| plot | `--color` | str |  | word | yes |
| plot | `--format` | str | png | word | yes |
| plot | `-l --log` | flag | False | word | yes |
| map | `-i --input` | str |  | word | no |
| map | `-q --query` | str |  | word | no |
| map | `-o --outdir` | str | . | word | no |
| map | `--prefix` | str | out | word | yes |
| map | `--mapper` | str | pblat | word | yes |
| map | `--minimap2-preset` | str | sr | kebab | yes |
| map | `-X --threads` | int | 1 | word | yes |
| map | `--minMatch` | int | 2 | camel | no |
| map | `--minScore` | int | 30 | camel | no |
| map | `--minIdentity` | percent | 90 | camel | yes |
| map | `--filterIdentity` | percent | 90 | camel | yes |
| map | `--minMatchCount` | int | 0 | camel | no |
| map | `--min-length` | int | 0 | kebab | yes |
| map | `--max-hits` | int |  | kebab | yes |
| map | `--orient-to-reference` | flag | False | kebab | yes |
| map | `--fasta-output` | str | mapped | kebab | yes |
| map | `-l --log` | flag | False | word | yes |
| check | `--alignments --psl` | Path |  | word | yes |
| check | `--forced_oligos` | Path |  | snake | yes |
| check | `--reference` | Path |  | word | yes |
| check | `--min_coverage` | float | 10.0 | snake | yes |
| check | `--min_length` | int | 100 | snake | yes |
| check | `--strand` | str | both | word | yes |
| check | `--min_similarity` | float | 95.0 | snake | yes |
| check | `--uncovered_length_cutoff` | int | 0 | snake | no |
| check | `--coverage_out` | Path |  | snake | yes |
| check | `--longest_uncovered_out` | Path | longest_uncovered.txt | snake | yes |
| check | `--uncovered_fasta` | Path |  | snake | yes |
| check | `--n_split_fasta` | Path |  | snake | yes |
| check | `--extend_region` | int | 0 | snake | yes |
| check | `--min_oligo_length` | int | 120 | snake | yes |
| check | `--temp_dir` | Path |  | snake | no |
| check | `--log-level --log_level` | str | INFO | kebab | yes |
| fill | `--alignments --psl` | Path |  | word | yes |
| fill | `--forced_oligos` | Path |  | snake | yes |
| fill | `--reference` | Path |  | word | yes |
| fill | `--min_coverage` | float | 1.0 | snake | yes |
| fill | `--spacing_distance` | int | 30 | snake | yes |
| fill | `--min_contribution` | int | 5 | snake | yes |
| fill | `--max_passes` | int | 5 | snake | yes |
| fill | `--max_oligos_per_pass` | int |  | snake | no |
| fill | `--stall_rounds` | int | 3 | snake | yes |
| fill | `--min_length` | int | 100 | snake | yes |
| fill | `--strand` | str | both | word | yes |
| fill | `--min_similarity` | float | 95.0 | snake | yes |
| fill | `--uncovered_length_cutoff` | int | 0 | snake | no |
| fill | `--output` | Path | selected_oligos.txt | word | yes |
| fill | `--coverage_out` | Path |  | snake | yes |
| fill | `--longest_uncovered_out` | Path |  | snake | yes |
| fill | `--uncovered_fasta` | Path |  | snake | yes |
| fill | `--n_split_fasta` | Path |  | snake | yes |
| fill | `--extend_region` | int | 0 | snake | yes |
| fill | `--min_oligo_length` | int | 120 | snake | yes |
| fill | `--temp_dir` | Path |  | snake | no |
| fill | `--log-level --log_level` | str | INFO | kebab | yes |
| evaluate | `-i --input` | str |  | word | no |
| evaluate | `-r --reference` | str |  | word | yes |
| evaluate | `-o --output` | str |  | word | yes |
| evaluate | `--mapper` | str | pblat | word | yes |
| evaluate | `--strand` | str | both | word | yes |
| evaluate | `--minimap2-preset` | str | sr | kebab | yes |
| evaluate | `--min-identity` | float | 90.0 | kebab | yes |
| evaluate | `--min-length` | int | 100 | kebab | yes |
| evaluate | `--threads` | int | 1 | word | yes |
| evaluate | `--min-coverage` | float | 1.0 | kebab | yes |
| evaluate | `--target-coverage` | float | 10.0 | kebab | yes |
| evaluate | `--min-gap-size` | int | 100 | kebab | yes |
| evaluate | `--gap-extend` | int | 0 | kebab | no |
| evaluate | `--plot-format` | str | png | kebab | yes |
| evaluate | `--plot-dpi` | int | 300 | kebab | no |
| evaluate | `--no-html-report` | flag | True | kebab | yes |
| evaluate | `--no-interactive-plots` | flag | True | kebab | yes |
| evaluate | `--offline-plots` | flag | False | kebab | yes |
| evaluate | `--reference-analysis-window` | int | 1000 | kebab | yes |
| evaluate | `--no-benchmarking` | flag | True | kebab | yes |
| evaluate | `--keep-intermediates` | flag | False | kebab | no |
| evaluate | `--quiet` | flag | False | word | yes |
| evaluate | `-l --log` | flag | False | word | yes |
| compare | `-r --reference` | str |  | word | yes |
| compare | `-o --output` | str |  | word | yes |
| compare | `--sets` | str |  | word | yes |
| compare | `--mapper` | str | pblat | word | yes |
| compare | `--strand` | str | both | word | yes |
| compare | `--minimap2-preset` | str | sr | kebab | yes |
| compare | `--min-identity` | float | 90.0 | kebab | yes |
| compare | `--min-length` | int | 100 | kebab | yes |
| compare | `--threads` | int | 1 | word | yes |
| compare | `--min-coverage` | float | 1.0 | kebab | yes |
| compare | `--target-coverage` | float | 10.0 | kebab | yes |
| compare | `--no-statistical-analysis` | flag | True | kebab | yes |
| compare | `--significance-level` | float | 0.05 | kebab | yes |
| compare | `--multiple-comparison-correction` | str | fdr | kebab | yes |
| compare | `--plot-format` | str | png | kebab | yes |
| compare | `--plot-dpi` | int | 300 | kebab | no |
| compare | `--quiet` | flag | False | word | yes |
| compare | `-l --log` | flag | False | word | yes |
