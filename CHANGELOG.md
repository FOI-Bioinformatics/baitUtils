# Changelog

## Unreleased

### Changed (reports and figures)

The report and figure redesign in `docs/report_design_audit_2026-10.md`
(D-numbers refer to that report).

- One design system (`baitUtils/design.py`) for every figure and both
  reports: a fixed-order categorical palette that passes the colour-vision
  checks, one sequential blue ramp for magnitudes, status colours reserved
  for the quality category and shown with their label, one reserved red for
  gaps, shared fonts, number formats, axis styling and a Plotly template
  (D-10 to D-13).
- Chart forms follow the data: the donut and the radar charts are gone
  (D-1, D-2); the compare matrix shows raw values shaded within each column
  (D-3); p-values are lollipops against alpha, effect sizes are signed bars,
  and tests on fewer than five references carry a caveat (D-4); the depth
  thresholds are the user's `--min-coverage` and `--target-coverage` (D-5);
  count axes have integer ticks and histograms round bin widths (D-6);
  empty panels say why they are empty (D-7); one colour scale per figure
  (D-8).
- New figures: `plots/coverage_profile` (depth along the references with
  gaps shaded), `plots/quality_score`, and a coverage map in
  `plots/coverage_heatmap`. `coverage_dashboard.html` and
  `comparison_dashboard.html` are coverage profiles;
  `metric_comparison.html` replaces `performance_radar.html` (D-9).
- The evaluate report has a per-reference results table (length, GC, N,
  breadth, breadth over non-N bases, depth, gaps, hits per strand,
  difficulty) (D-9).
- Both reports share one stylesheet and script, with a table of contents,
  numbered sections, numbered figures and tables with captions, summary
  tiles, a metadata line and a footer with tool versions and the command
  line; `templates/comparative_report.css` and `.js` are removed (D-14 to
  D-18).
- Plot ids are prefixed `fig-` so that they cannot collide with section
  ids, which displaced three figures above their headings (D-19); legends
  sit above the subplot titles and are omitted for single series (D-20).
- `CoverageVisualizer`, `InteractivePlotter` and `InteractiveReportGenerator`
  take `coverage_arrays`; the visualizer takes `quality_scores`; both report
  generators take `metadata` (inputs, command, versions).

### Code audit remediation

The remediation plan in `docs/code_audit_2026-10.md` (finding numbers below
refer to that report).

### Changed (performance and tests)

- Per-reference strand counts are gathered in one pass over the hits
  instead of one pass per reference (F-29).
- `evaluate` and `compare` read the reference FASTA once:
  `CoverageAnalyzer.sequences()` hands the loaded sequences to `GapAnalyzer`
  and `ReferenceAnalyzer` (both accept `reference_sequences`), and the gap
  analysis upper-cases only the gap slices it characterises instead of a
  copy of the whole reference (F-30).
- `plot --plot-type pca` drops rows with missing values (with a warning)
  instead of aborting (F-39).
- Test fixtures: the fake pblat and minimap2 only emit hits for the baits
  in the query FASTA, and the dataset has a second bait set
  (`baits_subset`), so `compare` is tested with sets that differ and the
  statistical tests, the Friedman path and the report's p-values are
  exercised end to end. The gap filling algorithm, the fake minimap2 `map`
  run and the plot command's PCA and SVG paths have tests.

### Changed (quality score and packaging)

- The uniformity half of the depth component is `1 - Gini` relative to the
  `depth_uniformity` target (0.8), as documented; the coefficient of
  variation is used only when no Gini coefficient is available (F-19).
- The reference difficulty component is the breadth relative to the breadth
  expected for the reference's difficulty (the breadth target, reduced by up
  to half for the hardest references), so the same breadth scores higher on
  a harder reference, as the docstring always said; without a reference
  analysis the reference counts as easy instead of scoring a fixed 0.7
  (F-20, F-45).
- Weight overrides are normalised to sum to 1, so the overall score stays
  in 0-1 (F-21).
- `recommendations.txt` lists the gap analysis suggestions and the quality
  scorer's recommendations, the same lists as the HTML report, instead of a
  third set of thresholds (F-43); `benchmark_analysis.txt` carries the
  overall score header (F-44); the benchmark reports uniformity as not
  applicable instead of a fabricated 0.5 when no Gini coefficient exists
  (F-47).
- The evaluate report's methodology section is generated from the analysis
  parameters, the mapper used and the scorer's constants (F-40).
- Packaging moved to `pyproject.toml` (`setup.py` is a shim); `tqdm` is no
  longer a dependency; `requirements.txt` matches the project dependencies
  and `pybedtools` is the `bedtools` extra; matplotlib >= 3.6 and
  seaborn >= 0.12 are required, and the category plots pass `hue` so that
  seaborn 0.14 keeps working; Python 3.13 is tested and declared (F-37,
  F-38, F-69).
- CI uses current action versions, tests Python 3.13 and installs bedtools
  on Ubuntu so that the `check` and `fill` end-to-end tests run (F-70).

### Changed (command line)

- Every option is spelled in kebab-case. The previous spellings still work
  as hidden aliases, so existing scripts keep running, but the help text,
  the documentation and the configuration-file keys use the new names:
  - `stats`: `--no-ns`, `--min-gc`, `--max-gc`, `--min-tm`, `--max-tm`,
    `--max-mask`, `--threads` (was `-p/--processes`; `-p` remains),
    `-o/--output` (was `--outdir`; the configuration key stays `outdir`).
    Configuration keys: `no_ns`, `min_gc`, `max_gc`, `min_tm`, `max_tm`,
    `max_mask`, `threads`.
  - `plot`: `--plot-type`, `-o/--output`.
  - `map`: `-r/--reference` (was `-q/--query`), `--threads` (was `-X`),
    `--min-match`, `--min-score`, `--min-identity`, `--filter-identity`,
    `--min-match-count` (were camelCase), `-o/--output`. Configuration keys
    are the snake_case forms (`min_identity`, `filter_identity`, ...).
  - `check` and `fill`: `--alignments` (`--psl` remains an alias),
    `--forced-oligos`, `--min-coverage`, `--min-length`, `--min-identity`
    (was `--min_similarity`), `--uncovered-length-cutoff`, `--coverage-out`,
    `--longest-uncovered-out`, `--uncovered-fasta`, `--n-split-fasta`,
    `--extend-region`, `--min-oligo-length`, `--temp-dir`; `fill` also
    `--spacing-distance`, `--min-contribution`, `--max-passes`,
    `--max-oligos-per-pass`, `--stall-rounds`. Configuration key for the
    alignment file: `alignments` (was `psl`); for identity: `min_identity`.
  - Every command accepts `-v/--verbose`, `--quiet` and `--log-level`
    (`-l/--log` remain aliases of `--verbose` where they existed), so
    `common.log_level` in a configuration file applies everywhere (F-9).
- Defaults aligned across commands (F-33): `check --min-coverage` is 1
  (was 10); `check` and `fill --min-identity` are 90 (was 95), the values
  `evaluate` and `map` use. `check --longest-uncovered-out` no longer
  defaults to a file in the working directory (F-64).
- `--config` is parsed with argparse, so `--config=file.json` works, the
  option may follow the command, and a missing file name is a usage error
  (F-8). Sections that are not commands and `common` keys that no command
  has are reported; a `common` key is applied silently to the commands
  that have it (F-9).
- Numeric options are validated by the parser (`--threads` at least 1,
  `--min-length` not negative, `--min-gap-size`, `--plot-dpi` and the
  analysis window at least 1).

### Added

- `evaluate` and `compare`: `--max-hits` (off-target filter as in `map`),
  `--temp-dir`; `compare`: `--min-gap-size`, `--gap-extend` and
  `--reference-analysis-window`, which were fixed at the `evaluate`
  defaults. `evaluate` writes the per-bait hit table (`hits.tsv`), a
  run-length `data/coverage.bedgraph` and `data/gap_regions.bed`, and with
  `--per-position-coverage` the per-base CSV; the report's data links point
  to these files (F-17).
- `map` warns when a pblat-only option (`--min-identity`, `--min-score`,
  `--min-match`) is given with `--mapper minimap2`, and when
  `--minimap2-preset` is given with pblat (also in `evaluate` and
  `compare`); the `--filter-identity >= --min-identity` check applies to
  pblat only (F-26). The filtered alignment file keeps the mapper's format
  (`-mapping_filtered.paf` for minimap2), and `--keep-intermediates` keeps
  the mapper output under its own name (F-25).
- `fill` warns that `--min-coverage` above 1 cannot be reached by stacking
  baits because of the spacing rule (F-31).
- `evaluate` always computes the reference analysis and the quality score;
  `--no-html-report`, `--no-interactive-plots` and `--no-benchmarking` skip
  only those outputs (F-18).
- `pytest.ini` is read (its section header was the setup.cfg form) and
  CLAUDE.md and CONTRIBUTING.md describe the pytest suite (F-34, F-35,
  F-36).

### Fixed

- `evaluate` and `compare` no longer crash on a single extreme-GC reference
  (division by zero in the reference difficulty score, F-1) or when no hit
  passes the identity, length and strand filters (coverage heatmap, F-2).
- `stats` no longer segfaults on an empty FASTA record when ViennaRNA is
  installed; empty sequences get `NA` energies (F-3).
- `compare` treated the 0-1 quality score as 0-10: the ranking now combines
  terms that are all in 0-1 (quality 40 %, breadth 25 %, mapping efficiency
  20 %, gap term 15 %), so the composite score is itself in 0-1; the HTML
  report prints the score on its own scale instead of "x/10" and the radar
  plot no longer divides it by 10 (F-4).
- The evaluate HTML report takes the quality category from the scorer
  (Excellent >= 0.85, Good >= 0.70, Fair >= 0.50) instead of its own
  0.8/0.6/0.4 thresholds, so report and `evaluation.json` agree (F-10).
- `evaluation.json` now contains `reference_summary` (it was always null
  because of a key mismatch) and a `reference_analysis` block with the
  per-reference features, challenging regions, the feature-coverage
  Spearman correlations with breadth by GC decile, and the challenging
  windows; `comparison.json` carries the same block for the shared
  reference. The evaluate HTML report gains the correlation, GC-decile and
  challenging-window tables (F-5).
- The gap analysis result includes `gap_features`, so the gap-feature panels
  of `detailed_gap_analysis.html` are populated; the reference feature
  heatmap of the HTML report reads `shannon_entropy` and scales each row to
  0-1 (F-22).
- The fabricated "Quality Trends" panel and the always-empty "Theoretical
  vs Actual" panel of `quality_assessment.html` are replaced by the overall
  score against the category thresholds (F-7).
- `compare --plot-format` and `--plot-dpi` are kept: the report embeds the
  plots generated by the command instead of regenerating them with default
  settings; SVG is embedded with its own media type and PDF is linked (F-6).
- Statistical test interpretations are written from the adjusted p-value,
  so text, stars and colour agree; star labels and the legend follow
  `--significance-level` (F-12, F-13).
- The paired per-reference gap-count test counts the gaps of the gap
  analysis (at or above `--min-gap-size`), the same definition as
  `total_gaps` and the gap-size test (F-14).
- Gap overlap between sets uses half-open intervals, so gaps that only touch
  are not reported as overlapping (F-15); links to the interactive plots in
  the comparative report point to their directory (F-16).
- Set names, reference ids and test text are HTML-escaped in both reports
  (F-11).
- Infinite values are written as `null` in the JSON outputs instead of the
  non-standard `Infinity` token (F-28).
- `stats` exits with status 1 when the input file is missing (F-27).
- `check --uncovered_fasta` and `fill --uncovered_fasta` write the uncovered
  FASTA without `--n_split_fasta`; the N-split file is optional (F-24).
- The comparative report's methodology section is generated from the
  scorer's weights and thresholds and names the tests that are run; its
  grade colours apply to the Excellent/Good/Fair/Poor categories (F-40).

### Added

- `docs/code_audit_2026-10.md`: code audit with gap and switch analysis.
- `tests/test_phase1_regressions.py`.

## 0.5.1 (2026-09-25)

### Added

- Coverage breadth over assessable (non-N) bases is reported alongside the
  breadth over the whole reference (`coverage_breadth_non_n`, `n_bases`,
  `assessable_bases`, overall and per reference) and in
  `coverage_statistics.txt`.
- `CoverageAnalyzer.export_coverage_data` writes a run-length
  `coverage.bedgraph` in addition to the per-position CSV and gap BED.

### Changed

- The per-position export is built with numpy instead of one dictionary
  per base.
- Reports, plots and console output use ASCII only (no emoji, ">=" instead
  of the symbol, "C" for degrees), which also removes matplotlib font
  warnings.

## 0.5.0 (2026-09-25)

### Added

- `--strand {both,plus,minus}` on `evaluate`, `compare`, `check` and `fill`
  restricts hits to one target strand; coverage statistics report hit
  counts per strand overall and per reference.
- `map --orient-to-reference` writes mapped baits in reference orientation
  by reverse-complementing baits whose best hit is on the minus strand.
- `--config file.json` in front of the command name supplies option
  defaults per command with a `common` block; command-line values win.
- `--log-level` in front of the command name; one logging format for all
  commands (`baitUtils.logging_utils`). `check` and `fill` accept
  `--log-level` (the old `--log_level` spelling still works).
- `evaluation.json` and `comparison.json` record the versions of baitUtils,
  Python, numpy, pandas, scipy, biopython, ViennaRNA, pblat, minimap2 and
  bedtools, and the full parsed argument set.
- Content tests for the evaluate and compare HTML reports and interactive
  plots, and tests for the configuration file, logging, strand handling
  and orientation.

## 0.4.0 (2026-09-25)

### Added

- minimap2 as an alternative mapper for `map`, `evaluate` and `compare`
  (`--mapper minimap2`, `--minimap2-preset`). PAF files from minimap2 -c are
  parsed into the same hit representation as PSL, including aligned blocks
  from the cigar and mismatches from the NM tag.
- `check` and `fill` accept PAF as well as PSL through `--alignments`
  (`--psl` remains as an alias); the format is chosen by file extension.
- Reference analysis reports real Spearman correlations between window
  features and observed coverage (with p-values and window counts), breadth
  by GC decile, and challenging windows (low breadth with an extreme
  feature). Window-level coverage is taken from the per-base arrays.

### Changed

- Reference analysis is vectorized (new `sequence_features` module) and
  runs in linear time; a 5 Mb reference takes about two seconds. Repeat
  content is the fraction of duplicated 12-mers rather than a heuristic
  scan. The fixed-formula "correlations" are removed. The
  `window_features` and `coverage_correlations` entries in
  `evaluation.json` have a new structure.
- Gap analysis repeat content uses the same k-mer routine instead of a
  quadratic loop.

## 0.3.0 (2026-09-25)

This release makes the documented functionality work as described. An audit
found that `evaluate` and `compare` crashed on every run, that several
reported values were simulated or wrong, and that the test suite mocked the
broken code paths. The changes below are grouped by their effect on users.

### Breaking

- `stats`: the `mfe` column is replaced by `hairpin_dg` and `self_dimer_dg`
  (kcal/mol, ViennaRNA DNA parameters at `--hyb-temp`, default 65 C).
- `map`: `-o` is now the output directory and `--prefix` the file name
  prefix (`-Z` is removed). A per-bait hit table `<prefix>-hits.tsv` is
  written.
- `check` and `fill`: `--reference` is required and replaces
  `--fasta_reference` and `--reference_sequence`. Reference sizes come from
  the FASTA, so unmapped tails are reported. `--max_coverage` (check, fill)
  and `--force` (fill) are removed; neither had an effect.
  `fill --longest_uncovered_out` no longer defaults to a file in the
  working directory; `fill` also writes `<output>_mappings.tsv`.
- `evaluate`: `--enable-html-report`, `--enable-interactive-plots` and
  `--enable-benchmarking` were always on and could not be disabled; they
  are replaced by `--no-html-report`, `--no-interactive-plots` and
  `--no-benchmarking`. `compare`: `--enable-statistical-analysis` is
  replaced by `--no-statistical-analysis`; `--keep-intermediates` is removed.
- Quality scores are reported on their real 0 to 1 scale with the
  categories Excellent, Good, Fair and Poor. Earlier output printed them as
  a 0 to 10 score with letter grades.
- The benchmark compares observed metrics with the design targets in
  `QualityScorer.DEFAULT_BENCHMARKS`; the "theoretical optimum" is removed.
- `plot --color` takes `kept` (the column the stats command writes).
- Forced-oligo IDs are matched case-sensitively.
- biopython 1.80 or later is required.

### Fixed

- `evaluate` and `compare` crashed at the reference analysis step
  (constructor and method name mismatches).
- `stats` reported the melting temperature as `NA` for every sequence; the
  salt concentration was passed as a string. `--mint` and `--maxt` now
  filter. `--na`, `--dnac1` and `--dnac2` expose the Tm parameters.
- Secondary structure energies used RNA parameters at 37 C; they now use
  DNA parameters at the hybridization temperature.
- Gap coordinates in `evaluate` were simulated; they are now taken from the
  per-base coverage arrays.
- Mapping efficiency was always 100 %; it now counts oligos in the input
  FASTA.
- `compare` ran its distribution tests on random samples and its metric
  tests on one value per group (always NaN). Tests now use observed
  window-level depth and paired per-reference metrics; cases without valid
  input are reported as not applicable. The multiple-comparison correction
  is applied and the configured significance level is used in reports.
- `check` and `fill` failed with `'NoneType' object is not callable` when
  pybedtools was missing; they now exit with an installation message.
- Four PSL parsers with two identity formulas are replaced by one parser
  using BLAT's identity calculation, so filters match pblat's own
  `-minIdentity`. Coverage uses aligned blocks rather than the target span.
- Uncovered-region merging dropped a run at the end of one reference when
  the next reference started covered.
- `fill` selects individual mappings, so selecting an oligo no longer adds
  its other loci silently.
- `check` and `fill` left temporary files in the working directory.
- Composition metrics (complexity, homopolymers, dinucleotide bias) ignore
  case; soft-masked input no longer skews them.

### Added

- `evaluate` writes `evaluation.json`; `compare` writes `comparison.json`.
- `map` writes a per-bait hit table and gains `--min-length` and
  `--max-hits` for off-target filtering.
- `stats` gains `--min-hairpin-dg` and `--min-dimer-dg` filters.
- `evaluate --offline-plots` embeds plotly.js in the HTML report.
- End-to-end CLI tests with a fake pblat and hand-derived fixtures; a CI job
  with pblat, bedtools and ViennaRNA.

### Removed

- `coverage_checking.py` (merged into `coverage_analysis.py`), the
  theoretical-optimum benchmark, unused scorer methods and unused imports.
