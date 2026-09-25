# Changelog

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
