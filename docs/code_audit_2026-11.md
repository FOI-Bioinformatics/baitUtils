# baitUtils 0.6.0: second code audit

Date: 2026-11-08. Audited revision: `f17b9bb` (main, after PR #20).

## 1. Scope and method

- Every module under `baitUtils/` (36 files, about 11 000 lines) was read
  by one of five reviewers, each owning one area: mapping, `stats` and
  `plot`; coverage, `check` and `fill`; the `evaluate` pipeline with
  scoring, JSON export, config and CLI; `compare` and its statistics; and
  the HTML reports, figures and `design.py`.
- Each reviewer was given the first audit (`docs/code_audit_2026-10.md`)
  so that fixed items were not reported again, and was asked to prove
  every finding by running code.
- The audit ran with the real external tools installed, unlike the first
  one: bedtools 2.31.1, pybedtools 0.12.1, minimap2 2.26 and ViennaRNA
  2.7.2. pblat was replaced by a script returning precomputed hits. The
  "no ViennaRNA" paths were tested with a module that fails to import.
- With the tools installed the suite passes (315 passed, 1 skipped), so
  none of the findings below is caught by a test.
- Every finding was then checked a second time, independently of the
  reviewer. The **Check** column says how:
  - **R**: the reproduction was re-run and gave the reported result.
  - **C**: the code was read and does what the finding says.
  - **M**: measured by the reviewer and not re-run (timings, browser
    layout, BLAT source behaviour).
- Two reviewers sometimes found the same bug; each is listed once. One
  reported crash needed a different input to trigger and is listed with
  that input. One apparent bug (a gap percentage of 50 000 %) came from a
  test harness that omitted `reference_length` and is not a bug.

Severity: **H** wrong or missing output on ordinary input, or a hang;
**M** wrong results, crashes or silently ignored options on realistic
input; **L** edge cases, misleading text, latent problems.

## 2. Summary

| Area | H | M | L | Total |
|------|---|---|---|-------|
| Mapping and alignment parsing | 0 | 1 | 4 | 5 |
| `stats` and `plot` | 0 | 3 | 5 | 8 |
| `check` and `fill` | 1 | 4 | 5 | 10 |
| Reference input and coverage | 0 | 1 | 3 | 4 |
| `evaluate`: N bases, gaps, reference analysis, score | 0 | 5 | 5 | 10 |
| CLI and config | 0 | 2 | 1 | 3 |
| `compare` | 0 | 4 | 6 | 10 |
| Reports and figures | 1 | 2 | 8 | 11 |
| **Total** | **2** | **22** | **37** | **61** |

The most consequential findings:

1. **N-51** `evaluate` effectively hangs on references with thousands of
   gaps: the interactive coverage profile adds one Plotly shape per gap,
   which is quadratic, and is built twice.
2. **N-14** `fill` stops while coverage is still improving, because its
   stall counter watches the number of uncovered regions, which rises
   when a gap is split.
3. **N-1** every real pblat output triggers a "malformed line" warning,
   and the filtered PSL keeps only two of the five header lines.
4. **N-28 to N-30** N runs are treated as sequence: all-N windows count as
   0 % GC and zero entropy, which creates false GC-coverage correlations,
   wrong gap labels and suggestions, and a lower quality score for
   assemblies with scaffold gaps.
5. **N-32** `evaluate --mapper minimap2` fails when no bait maps.
6. **N-41 and N-42** two `compare` tests answer the wrong question: the
   identity test counts baits that `--max-hits` excluded, and the window
   tests ignore that windows are paired.
7. **N-38 and N-39** config values skip all validation, and `--conf`
   (an abbreviation) silently ignores the file.
8. **N-52** a `$` in a reference name crashes `evaluate` with no output.

## 3. Findings

### 3.1 Mapping and alignment parsing

| # | Sev | Location | Finding | Check |
|---|-----|----------|---------|-------|
| N-1 | M | `mapping_utils.py:24`, `:124-125`, `:493-500` | The fourth line of the psLayout header that pblat writes begins with spaces and a tab, so `is_psl_header` misses it. `parse_psl` warns "Skipped 1 malformed line" on every real pblat file (in `map`, `check`, `fill`, `evaluate` and `compare`), and the header copy stops there, so `*_filtered.psl` keeps 2 of 5 header lines. UCSC tools that skip a fixed header would then drop data lines. The test fixtures never write a header. Fix: treat a line whose first field is not an integer as header, and copy every line before the first data line. | R |
| N-2 | L | `sequence_mapping.py:53`, `:229-237`; `config.py:62-72` | `map --min-identity` sets pblat's seed `-minIdentity`, while every other command uses `min_identity` as the identity filter. A `common` config with `"min_identity": 95` makes `map` fail ("filter_identity (90.0) must be >= min_identity (95)"). Fix: give the pblat option its own destination, or make `map --min-identity` set the filter. | R |
| N-3 | L | `mapping_utils.py:137`, `:230`, `:243-256` | `alignment_format` accepts `.psl.gz` and `.paf.gz`, but the parsers use plain `open`, so `check --alignments hits.paf.gz` fails with a UTF-8 decode error. Fix: `gzip.open(..., 'rt')` for `.gz`. | R |
| N-4 | L | `mapping_utils.py:90-96` | BLAT truncates milliBad to an integer before computing identity; baitUtils keeps the float. A hit pblat reports at exactly 90 % (milliBad 100.92) is 89.91 % here and is dropped by `--filter-identity 90`. Fix: floor the milliBad. | C, M |
| N-5 | L | `mapping_utils.py:186-196` | PAF query block starts are wrong on the minus strand (for example `[0, 40, 79]` instead of `[79, 42, 0]`). Nothing reads them today. | M |

### 3.2 `stats` and `plot`

| # | Sev | Location | Finding | Check |
|---|-----|----------|---------|-------|
| N-6 | M | `sequence_statistics.py:84-99` | Tm, hairpin and dimer filters pass sequences whose value is NA. With `--min-tm 60 --max-tm 90 -f`, an empty record and an all-N record end up in `filtered_sequences.fasta`. Fix: a sequence whose value cannot be computed fails a filter the user set; log the count. | R |
| N-7 | M | `sequence_statistics.py:96-99`, `:350-355` | Without ViennaRNA, `--min-hairpin-dg` and `--min-dimer-dg` are silently ignored (7 of 7 pass; 3 of 7 with ViennaRNA). Fix: exit with an error, or warn, when these options are set and ViennaRNA is missing. | R |
| N-8 | M | `plotting_utils.py:40-43`, `:207-215`, `:404-407` | Without ViennaRNA the energy columns are all NaN but still selected: `plot --color kept` exits 1 ("'dataset' input should have multiple elements"), and `plot --plot-type pca` drops every row, logs an error, prints "All plots have been generated successfully" and exits 0 with no files. The first audit's fix (F-39) only handled scattered NaN. Fix: skip all-NaN columns with a warning, and return a failure when PCA cannot run. | R |
| N-9 | L | `sequence_statistics.py:224-229` | When no sequence passes, `filtered_sequences.fasta` is not rewritten, so a previous run's file (43 records) stays next to a log saying 0 passed. An empty input writes no table or summary and exits 0. Fix: always write the outputs, empty if need be, and warn on zero records. | R |
| N-10 | L | `sequence_analysis.py:141-146` | `self_alignment_score` is case-sensitive: a bait with a soft-masked half scores -12 where the upper- and lowercase versions both score 66. Fix: upper-case first, as the other metrics do. | R |
| N-11 | L | `sequence_analysis.py:76-83` | Tm ignores N: Biopython strips N, so a bait with 10 or 60 N in the middle gets the Tm of the joined flanks (76.75 in all three cases). Fix: report NA, or warn, for non-ACGT sequence. | R |
| N-12 | L | `sequence_analysis.py:207-241`, `:317-319` | N runs count as homopolymers (`max_homopolymer_length` 10 for a 10 N run), unlike `sequence_features.max_homopolymer`. Fix: count only A, C, G and T runs. | R |
| N-13 | L | `plotting_utils.py:404`, `statistical_plots.py:41-43` | Numeric sequence IDs are read as numbers and plotted as a statistic (`histogram_sequence_id.png`), and would enter PCA and the heatmap. Fix: read `sequence_id` as text. | R |

### 3.3 `check` and `fill`

| # | Sev | Location | Finding | Check |
|---|-----|----------|---------|-------|
| N-14 | H | `gap_filling_algorithm.py:305-317` | The stall counter rises whenever the number of uncovered regions does not fall. Filling the middle of a gap splits it, so good passes count as stalls: with `--max-oligos-per-pass 5` uncovered bases go 5000, 4400, 3800, 3200 and `fill` stops after 15 baits "after 3 rounds without improvement", with 64 % still uncovered. The first audit noted the loop logic (F-32) but it was not changed. Fix: count a stall only when uncovered bases do not fall. | R |
| N-15 | M | `gap_filling_algorithm.py:281-282`, `gap_filling.py:120`; `coverage_analysis.py:355` | Forced baits without a passing hit are left out of the `fill` output with no warning, although the README says forced baits are always included. `check` also ignores unmatched forced IDs silently. Fix: add all forced IDs to the output and warn about those without a hit. | R |
| N-16 | M | `gap_filling_algorithm.py:120` | `fill --min-contribution 0` is accepted and crashes with "math domain error" (`log2(0)` for a candidate with no overlap). Fix: score a zero-overlap candidate 0, or require a positive value. | R |
| N-17 | M | `gap_filling_algorithm.py:54`, `:103`, `:118` | Candidates are scored on their whole span, including target inserts, not on aligned blocks, so a gapped hit gets credit for bases it does not cover. In a crafted case a bait is selected that leaves uncovered bases unchanged. Fix: compute overlap from `m.blocks`. | R |
| N-18 | M | `coverage_analysis.py:467`, `gap_filling.py:177` | `--n-split-fasta` writes nothing unless `--uncovered-fasta` is also given; the first audit's fix (F-24) covered only the other direction. Fix: export when either path is set. | R |
| N-19 | L | `coverage_analysis.py:93`, `:170`, `:182` | bedtools treats lines starting with `track` or `browser` as headers, so a reference named `track_1` drops out of the coverage output of `check` and `fill`: with two 500 bp contigs and one 120 bp hit on each, `check` reports 380 bases under minimum coverage instead of 760, and `--coverage-out` has no `track_1` rows. `evaluate` is unaffected, so the commands disagree. Fix: reject such names with a clear error, or compute coverage with numpy. | R |
| N-20 | L | `coverage_analysis.py:109-121`; `coverage_stats.py` | A hit to a reference missing from the FASTA makes `check` fail with a bedtools error and `std::bad_alloc` (F-66, still open). `evaluate` drops such hits from coverage but still counts them in `mapped_oligos` and `strand_counts`. Fix: check names up front and warn and skip, or fail naming them. | R |
| N-21 | L | `coverage_analysis.py:167`, `gap_filling.py:140` | `fill` counts coverage only at the selected locus of each bait, but its output is a list of baits, which capture every locus. On a test, `fill` reports 10 bp uncovered and `check --forced-oligos` on the same list reports 0. Fix: count every mapping of a selected bait, at least in the final report. | R |
| N-22 | L | `gap_filling_algorithm.py:74` | The GC term of the score is case-sensitive, so soft-masked reference scores as 0 % GC (13 385 vs 12 213 for the same sequence). Fix: upper-case first. | R |
| N-23 | L | `gap_filling.py:79-84` | "Total oligos for selection: 40" counts BED blocks; there are 39 hits. | R |

A performance note that is plausible but not a correctness bug: each
`fill` pass scans every uncovered region for every mapping in pure Python,
about 1.6 s per pass for 10 000 mappings and 1 000 regions (M).

### 3.4 Reference input and coverage (shared by `evaluate`, `compare`, `check`)

| # | Sev | Location | Finding | Check |
|---|-----|----------|---------|-------|
| N-24 | M | `coverage_stats.py:116-121` | Duplicate reference IDs are merged silently: two records named `dup` (3000 and 1000 bp) become one reference of 1000 bp. Fix: reject duplicate IDs when loading. | R |
| N-25 | L | `coverage_stats.py:381`; `coverage_analysis.py:118-121` | An empty FASTA record crashes `evaluate` (`zero-size array to reduction operation maximum`) and `check` (bedtools "length equal to 0"). Fix: skip empty records with a warning. | R |
| N-26 | L | `evaluate.py:374`; `coverage_stats.py:320-351` | `--min-coverage 0` is accepted: with no hits counted, breadth is 100 % and the Gini coefficient is NaN. Fix: require a value above 0. | R |
| N-27 | L | `coverage_stats.py:322-326`, `benchmark.py:88-91`, `quality_scorer.py:360-364` | When nothing is covered, Gini is set to 1.0 and CV to infinity, so the benchmark's "not applicable" branch (the F-47 fix) is never reached: the report shows uniformity "Poor" with a redistribution advice and "Coverage depth is highly variable (CV=inf)". Fix: set both to None when nothing is covered and skip non-finite values. | R |

### 3.5 `evaluate`: N bases, gaps, reference analysis and score

| # | Sev | Location | Finding | Check |
|---|-----|----------|---------|-------|
| N-28 | M | `sequence_features.py:41-60`, `:181-184`; `reference_analyzer.py:144-148`, `:166-197` | A window with no A, C, G or T gets GC 0 % and entropy 0 instead of "undefined", and its breadth is 0 because N cannot be covered. N windows therefore produce false correlations (GC vs breadth rho 0.76 on data where coverage does not depend on GC), fill a "0-10 %" GC decile, and are flagged "extreme GC" and "low complexity". Window breadth also counts N in its denominator. Fix: NaN for GC and entropy without ACGT, drop such windows from the correlations, deciles and flags, and compute window breadth over non-N bases. | R |
| N-29 | M | `gap_analysis.py:173`, `:200-260` | Gap features count N as a base: a gap over 500 N has GC 0, a 500 bp homopolymer and "extreme composition", and the suggestions say all gaps have extreme GC and low complexity. Fix: compute gap features over ACGT only, with NaN when none remain. | R |
| N-30 | M | `quality_scorer.py:194`, `:254-256`; `gap_analysis.py:123`, `:309`; `reference_analyzer.py:211-216` | The score penalises N runs, which cannot be covered: a 20 kb reference with a 5 kb N run and every other base tiled (non-N breadth 99.6 %) gets breadth component 0.79 and gap component 0.47 (the largest gap is 5040 bp, 5000 of them N). An all-N contig is "Very Difficult" because of "Extreme GC content (0.0 %)". Fix: score breadth on non-N bases, keep N out of gap length and percentage (or report N gaps separately), and treat GC as undefined without ACGT. | R |
| N-31 | M | `gap_analysis.py:224-226`, `:386`, `:476` | Gap "repeat content" is the share of duplicated 6-mers, which grows with length on random sequence (15 % at 500 bp, 20 % at 1 kb, 38 % at 2 kb, 70 % at 5 kb). Every random gap of about 1.5 kb or more crosses the 30 % "high repeat" line and triggers a masking suggestion. Fix: use the reference analyzer's 12-mer measure or compare with the value expected for random sequence of that length. | R |
| N-32 | M | `evaluate.py:397-399` | With `--mapper minimap2`, `evaluate` exits 1 ("minimap2 produced no output") when no bait maps, because minimap2 writes an empty PAF. pblat writes a header and the same input works. Fix: treat an empty file from a mapper that exited 0 as no hits. | R |
| N-33 | L | `gap_analysis.py:168-170` | An empty gap sequence returns a feature dict without `extreme_composition`, so `evaluate` fails with `KeyError` when `--gap-extend` is negative and larger than half a gap. Only reachable through a config file today (see N-38). Fix: return the full key set. | R |
| N-34 | L | `gap_analysis.py:319-336` vs `:356-360` | The two gap size classifications disagree at the edges: `[lo, hi)` in one, `pd.cut`'s `(lo, hi]` in the other, so 500 and 1000 bp gaps land in different classes. With `--min-gap-size 1000` the first class is still "1000-500" (F-55). Fix: one set of half-open edges starting above the minimum gap size. | C |
| N-35 | L | `evaluate.py:491` | `gap_analysis.txt` writes 0-based half-open gaps as `chrA:1080-1380`, a notation read as 1-based closed by samtools and IGV. Fix: print `chrA:1081-1380` or state the convention. | R |
| N-36 | L | `quality_scorer.py:134-141` | A misspelt weight key (code API only) is accepted and caps the overall score: perfect statistics score 0.667 "Fair". Fix: reject unknown keys. | R |
| N-37 | L | `json_export.py:35-36` | `to_jsonable` passes NaN and infinity through inside numpy arrays (`[NaN, Infinity]`), which is invalid JSON. No array reaches the JSON files today (latent; F-28 was fixed for scalars). Fix: convert `tolist()` results recursively. | R |

### 3.6 CLI and config

| # | Sev | Location | Finding | Check |
|---|-----|----------|---------|-------|
| N-38 | M | `config.py:172-174`, `:128-134`, `:165` | Config values are applied with `set_defaults` and skip `type` and `choices`: `"strand": "reverse"` runs on both strands and is recorded as "reverse"; `"enable_html_report": "false"` still writes the report; `"gap_extend": -200` can crash (N-33); `"plot_dpi": 0` fails every figure but exits 0; a string `sets` value is split into characters. A missing file, invalid JSON or a non-object section gives a raw traceback. Fix: validate config values with the option's `type` and `choices` after parsing, and report file errors cleanly. | R |
| N-39 | M | `__main__.py:37`, `:57-59`, `:74` | The pre-parser that reads `--config` has `allow_abbrev=False`, but the main parser accepts `--conf` as an abbreviation and then overwrites it, so `baitUtils --conf cfg.json evaluate ...` runs without the config and without a warning. Fix: `allow_abbrev=False` on the main parser. | R |
| N-40 | L | `cli_utils.py:199-201` | The hidden alias spellings make prefixes ambiguous: `map --filter 95` fails with "could match --filter-identity, --filterIdentity", which are the same option. Fix: resolve prefixes whose matches share one action, or disallow abbreviations. | R |

### 3.7 `compare`

| # | Sev | Location | Finding | Check |
|---|-----|----------|---------|-------|
| N-41 | M | `differential_analysis.py:245-249` | The bait-identity test uses every bait in the hit table, including those `--max-hits` turned into unmapped. With `--max-hits 0` both sets have 0 mapped baits, yet the test reports n = 78, p = 1e-8, "***". Fix: drop rows flagged `multi_mapping_excluded`. | R |
| N-42 | M | `differential_analysis.py:120-146` | The window-depth tests (KS, Mann-Whitney, Levene) treat the two sets' windows as independent, although they are the same windows. With one set 1x deeper in every window, all three report "No difference" (p = 1.0, 0.705, 1.0); a paired Wilcoxon on the same windows gives p = 8e-6. Fix: test location with a paired signed-rank test, and present KS and Levene as descriptive or note that they ignore pairing. | R |
| N-43 | M | `comparative_analyzer.py:167-177`, `compare.py:305` | `comparison.json` stores one `reference_analysis`, computed from the first set's coverage, as if it applied to all sets: reordering `--sets` changes `breadth_by_gc` (89 vs 71). Fix: keep the sequence-only part shared and store the coverage correlations per set. | R |
| N-44 | M | `comparative_analyzer.py:324` | The ranking's gap term is `(100 - total_gaps) / 100`, so a set where no bait maps (one gap per reference) gets nearly the full 15 %, more than a well-covered design with many small gaps. Fix: base the term on the share of the reference in gaps. | R |
| N-45 | L | `differential_analysis.py:134-138`, `:231` | A Mann-Whitney test with all values tied is reported as applicable with p null and adjusted p 1.0, and still counts in its correction family. Fix: "Not applicable" when all values are tied or p is not finite. | R |
| N-46 | L | `comparative_analyzer.py:258-259`, `:445` | `Uniformity_Diff` in `pairwise_comparisons.csv` is the difference of Gini coefficients, so its sign is the opposite of its name. Fix: compute the difference of `1 - Gini`, or rename the column. | C |
| N-47 | L | `comparative_report_generator.py:390-399` | Recommendations name the best set by quality score and then `ranking[1]`; when the two orders differ, the same set is named twice and the composite winner not at all. Fix: base the text on the ranking and add a separate line for the best quality score when it differs. | C |
| N-48 | L | `comparative_analyzer.py:364-368` | "Distinct gap regions" counts exact coordinate tuples, so the same hole with different edges in two sets counts twice. Fix: merge overlapping intervals across sets before counting. | C |
| N-49 | L | `comparative_visualizations.py:217-230` | The effect-size panel says "positive: first set larger" for measures that cannot be negative (Levene ratio, KS D, Kendall's W). Fix: plot only signed effects there. | C |
| N-50 | L | `compare.py:330-375` | Empty or whitespace set names (`:baits.fa`) are accepted and appear blank in the ranking and report. Fix: reject them and compare stripped names for duplicates. | R |

Two cosmetic points in the same module: the Mann-Whitney interpretation
rounds a half-integer U (`U=14` next to 14.500), and the Wilcoxon `n`
counts references whose difference is zero, which scipy drops.

### 3.8 Reports and figures

| # | Sev | Location | Finding | Check |
|---|-----|----------|---------|-------|
| N-51 | H | `interactive_plots.py:127-129`; `report_generator.py:220` | The interactive coverage profile calls `fig.add_vrect` once per gap, which is quadratic in Plotly (0.8 s for 50 gaps, 2.6 s for 100, 10.8 s for 200), and `evaluate` builds the figure twice. A 1 Mb reference with 2778 gaps had not finished after 110 s; the estimate is about an hour per figure. Fix: build all shapes as dicts and set them once (0.6 s for 2778), or draw gaps as one filled trace per reference. | R |
| N-52 | M | `coverage_viz.py:162`, `:165`, `:177`; `design.py:165` | Matplotlib parses `$...$` in reference names as mathtext: `ctg$^$x` crashes `evaluate` with no report and no `evaluation.json`, and `seq$1$_a` renders as math. One failing figure also aborts the whole run because `generate_all_plots` does not catch per-figure errors. Fix: `"text.parse_math": False` in the design style, and catch and log per-figure failures. | R |
| N-53 | M | `interactive_plots.py:224` | In "Largest gaps", labels are truncated names plus coordinates; when truncation makes them equal, Plotly merges them into one category and stacks the bars, showing a single 1620 bp bar for three 540 bp gaps. Fix: unique category values with readable tick text. | R |
| N-54 | L | `report_generator.py:207` | The Gaps tile says "300 bp uncovered" while `uncovered_bases` is 440, because it counts only gaps at or above `--min-gap-size`. Fix: "300 bp in gaps", or show `uncovered_bases`. | R |
| N-55 | L | `interactive_plots.py:223`, `:243` | A histogram in the same figure makes Plotly's default bar gap 0, so the "Largest gaps" bars touch and equal-length gaps look like one block. Fix: set `bargap` explicitly or draw the histogram as bars. | C, M |
| N-56 | L | `coverage_viz.py:264-268` | The static depth histogram stops at the 99.5th percentile without saying so; a 0.4 % mass at depth 50 disappears. Fix: an overflow bin or a note. | C |
| N-57 | L | `coverage_viz.py:284`, `:299` | "Depth per reference" shows the 25 longest references without the "(25 longest of N)" note the other figures carry. | C |
| N-58 | L | `interactive_plots.py:117`, `:205-206`, `:213`, `:224`, `:235`, `:347`, `:372` | Plotly interprets `<b>`, `<br>` and entities in names used as titles, legend entries, ticks and hover text, so `chr<b>A</b>` renders bold. Scripts do not execute and the HTML tables are escaped. Fix: escape names passed to Plotly text. | M |
| N-59 | L | `interactive_plots.py:159-162`, `coverage_viz.py:228-232` | The cumulative depth curve makes one pass over all bases per depth value, three times per run: 1.38 s per call for 10 Mb at 150x, against 0.04 s with `bincount` and a reverse cumulative sum. | M |
| N-60 | L | `interactive_plots.py:95` | `--offline-plots` embeds plotly.js in the report but not in `interactive_plots/*.html`, which always load it from the CDN and are blank offline. | C |
| N-61 | L | `design.py:108-129`; `report_generator.py:200`, `:202`, `:309`, `:342`; `interactive_plots.py:303`; `coverage_viz.py:483` | `fmt_bp` picks the unit before rounding (999.5 gives "1000 bp", 999 500 gives "1000 kb"), and NaN gives "nan bp". The report and figures crash on `None` values, which is what NaN becomes when the dictionaries are reloaded from `evaluation.json`. Fix: round, then choose the unit; treat `None` like NaN. | R |

## 4. Items from the first audit that are still open

These were reported in October, are not marked fixed in the CHANGELOG, and
were seen again: F-23 (dead correlation code), F-41 (`PSLParser.last_hit_table`
is a class attribute), F-46 (unclipped benchmark ratios), F-48 (window
repeat density above 100 %), F-50 (hard-coded 10 % homopolymer flag),
F-51 (`applicable` true when every rho is NaN), F-52 (`stats` entropy counts
N), F-53 (only `ImportError` caught around the ViennaRNA loader) and F-56
(pure-Python low-complexity loop). F-24, F-32, F-39, F-47, F-55 and F-66
are partly or wholly open again as N-18, N-14, N-8, N-27, N-34 and N-20.

The `plot` command does not use `design.py`; its colours are hard-coded,
contrary to CLAUDE.md. This is a consistency issue, not a wrong result.

## 5. Proposed remediation

Grouped so that each phase can ship on its own, most valuable first.

1. **Hangs, crashes and silent failures on real data.** N-51 (shapes in
   one call), N-1 (PSL header), N-14 (fill stall), N-32 (empty PAF),
   N-52 (mathtext and per-figure isolation), N-24 and N-25 (duplicate and
   empty records), N-16 (min-contribution 0), N-15 (forced baits), N-6,
   N-7 and N-8 (`stats` filters and `plot` without ViennaRNA), N-38 and
   N-39 (config validation). Each gets a regression test that fails today.
2. **N bases and evaluate semantics.** One shared rule, "statistics over
   A, C, G and T; undefined without them", applied in
   `sequence_features`, `gap_analysis`, `reference_analyzer` and the
   scorer: N-28, N-29, N-30, N-12, N-11, F-52. Then N-31 (repeat measure),
   N-27 (nothing covered), N-26, N-33, N-34, N-35.
3. **`compare` statistics.** N-41, N-42 (paired window test; a design
   change that alters reported p-values, so it needs a CHANGELOG note),
   N-43, N-44 (ranking change, also user-visible), N-45 to N-50.
4. **`check` and `fill` correctness.** N-17, N-18, N-19, N-20, N-21,
   N-22, N-23, N-3.
5. **Smaller items and performance.** N-2, N-4, N-9, N-10, N-13, N-36,
   N-37, N-40, N-53 to N-61, and the open F items in section 4.

Phases 2 and 3 change numbers that users compare across runs (quality
scores of assemblies with N runs, `compare` p-values and ranking), so they
belong in a minor release with the change spelled out in the CHANGELOG.
