# baitUtils

[![CI](https://github.com/FOI-Bioinformatics/baitUtils/workflows/CI/badge.svg)](https://github.com/FOI-Bioinformatics/baitUtils/actions)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![codecov](https://codecov.io/gh/FOI-Bioinformatics/baitUtils/branch/main/graph/badge.svg)](https://codecov.io/gh/FOI-Bioinformatics/baitUtils)

Command-line tools for designing and evaluating oligonucleotide baits for
in-solution hybridization capture. baitUtils computes sequence properties of
bait sets, maps baits to a reference, quantifies coverage and gaps, selects
additional baits to close gaps, and compares alternative designs.

See [CHANGELOG.md](CHANGELOG.md) for the changes in 0.3.0, several of which
change output file names and options.

## Installation

Requirements: Python 3.10 or later. Two external tools are needed for parts
of the workflow:

- pblat or minimap2 for `map`, `evaluate` and `compare` (`--mapper`)
- bedtools (with the pybedtools package) for `check` and `fill`
- ViennaRNA (optional) for the hairpin and dimer energies in `stats`

```bash
conda create -n baitutils -c conda-forge -c bioconda python=3.12 pblat minimap2 bedtools pybedtools viennarna
conda activate baitutils
pip install baitutils          # or: pip install -e .  from a clone
```

Python dependencies (installed by pip): numpy, pandas, matplotlib, seaborn,
scikit-learn, biopython (1.80 or later), plotly, scipy, tqdm.

## Commands

| Command    | Purpose                                                      | Needs           |
|------------|--------------------------------------------------------------|-----------------|
| `stats`    | Per-bait sequence statistics and filtering                   | ViennaRNA (opt) |
| `plot`     | Plots from a statistics table                                |                 |
| `map`      | Map baits to a reference; per-bait hit table                 | pblat/minimap2  |
| `check`    | Coverage and uncovered regions from a PSL or PAF file        | bedtools        |
| `fill`     | Greedy multi-pass selection of baits to close gaps           | bedtools        |
| `evaluate` | Mapping, coverage, gaps, quality score and HTML report       | pblat/minimap2  |
| `compare`  | Evaluate several bait sets and test differences between them | pblat/minimap2  |

`baitUtils <command> --help` lists all options.

### stats

```bash
baitUtils stats -i baits.fasta -o stats/ --filter \
  --mingc 40 --maxgc 60 --mint 60 --maxt 80 --min-dimer-dg -15
```

Writes `sequence_statistics.tsv` (one row per bait), `summary.txt` and, with
`--filter`, `filtered_sequences.fasta` containing the baits that pass. All
baits stay in the table with a `kept` column.

Columns: `length`, `gc_content`, `melting_temperature` (nearest-neighbour Tm;
salt and strand concentrations from `--na`, `--dnac1`, `--dnac2`),
`hairpin_dg` and `self_dimer_dg` (kcal/mol at `--hyb-temp`, default 65 C,
`NA` without ViennaRNA), `entropy`, `complexity_2mer`, `complexity_3mer`,
`self_alignment_score`, `dinucleotide_bias`, `n_count`, `masked_count`,
`masked_percentage`, `homopolymer_runs`, `max_homopolymer_length`.

Filters: `--length` with `--complete`, `--noNs`, `--mingc/--maxgc`,
`--mint/--maxt`, `--maxmask`, `--min-hairpin-dg`, `--min-dimer-dg`. Free
energies are negative; a filter rejects baits whose value is below the given
threshold.

### plot

```bash
baitUtils plot -i stats/sequence_statistics.tsv -o plots/ \
  --columns gc_content melting_temperature hairpin_dg \
  --plot_type histogram boxplot scatterplot --color kept
```

Plot types: histogram, boxplot, violinplot, densityplot, boxenplot, swarmplot,
scatterplot, jointplot, heatmap, pairplot, pca. Output format with `--format`.

### map

```bash
baitUtils map -i baits.fasta -q reference.fasta -o mapping/ --prefix run \
  --threads 4 --minIdentity 90 --filterIdentity 95 --min-length 100 --max-hits 1
```

Runs the mapper (`-i` baits, `-q` reference) and writes, under `-o` with the
`--prefix`:

- `run-mapping.psl` (pblat) or `run-mapping.paf` (minimap2) and `run-mapping_filtered.*`
- `run-hits.tsv`: per bait, number of hits and targets, best and second-best
  identity, best locus and strand
- `run-mapped-sequence-ids.txt`, `run-unmapped-sequence-ids.txt`
- FASTA of mapped and/or unmapped baits (`--fasta-output`)

Identity follows BLAT's calculation, so `--filterIdentity` is on the same
scale as pblat's `-minIdentity`. `--max-hits` treats baits with more passing
hits as unmapped, which removes multi-mapping baits.

`--mapper minimap2` runs minimap2 with base-level alignment (`-c`, so
cigars and edit distances are available), secondary hits kept, and the
preset from `--minimap2-preset` (default `sr`). PAF input is converted to
the same hit representation as PSL, so all downstream numbers are
comparable between the two mappers apart from differences in the
alignments themselves.

In our tests pblat occasionally returned one hit fewer when run with several
threads on a very small input, and once terminated with a signal. `map`,
`evaluate` and `compare` retry once with a single thread after a signal;
for small inputs `--threads 1` is the safer choice.

### check

```bash
baitUtils check --alignments mapping/run-mapping.psl --reference reference.fasta \
  --min_coverage 1 --min_similarity 95 --min_length 100 \
  --longest_uncovered_out uncovered.tsv --coverage_out coverage.tsv \
  --uncovered_fasta uncovered.fasta --n_split_fasta uncovered_split.fasta
```

Converts PSL or PAF hits (aligned blocks) to BED, computes depth with bedtools and
reports runs below `--min_coverage`. Reference sizes come from the FASTA, so
regions after the last mapped bait are included. `--forced_oligos` restricts
the check to a list of bait IDs. `--uncovered_fasta` exports the uncovered
sequence, optionally extended by `--extend_region` and split at N runs into
pieces of at least `--min_oligo_length`.

### fill

```bash
baitUtils fill --alignments mapping/run-mapping.psl --reference reference.fasta \
  --output selected.txt --min_coverage 1 --spacing_distance 30 \
  --min_contribution 5 --max_passes 5
```

Greedy multi-pass selection: each pass scores every unselected mapping by
its overlap with uncovered regions (weighted by region length, local
difficulty and reference GC), and accepts candidates subject to
`--spacing_distance` between start positions. Selection is per mapping, so a
bait with several loci contributes only the locus that was chosen. Stops at
`--max_passes`, after `--stall_rounds` passes without improvement, or when the
uncovered base count stops decreasing. Baits in `--forced_oligos` are always
included.

Writes the selected bait IDs to `--output`, the selected loci to
`<output>_mappings.tsv`, and the same optional coverage and uncovered-region
outputs as `check`.

### evaluate

```bash
baitUtils evaluate -i baits.fasta -r reference.fasta -o evaluation/ \
  --min-identity 95 --min-coverage 1 --target-coverage 10 --min-gap-size 100 \
  --threads 4 --plot-format pdf --offline-plots
```

Output:

```
evaluation/
  evaluation.json                  all results, machine readable
  coverage_statistics.txt          breadth, depth, mapping efficiency
  gap_analysis.txt                 gap summary and largest gaps
  recommendations.txt
  benchmark_analysis.txt           observed metrics against design targets
  coverage_evaluation_report.html  report (plotly.js from CDN unless --offline-plots)
  plots/                           static figures
  interactive_plots/               plotly figures
```

Mapping efficiency counts baits in the input FASTA. Gap coordinates come from
the per-base coverage arrays.

Reference analysis computes sequence features per reference and per window
(`--reference-analysis-window`, default 1000 bp): GC, N content, entropy,
homopolymer content and repetitive content (fraction of duplicated
12-mers). Each window also carries its observed breadth and mean depth, and
the report gives Spearman correlations between features and coverage across
windows with p-values, breadth by GC decile, and the windows with low
breadth and an extreme feature. All of this is linear in reference length. The quality score is a weighted sum of five
component scores in 0 to 1 (breadth, depth, mapping efficiency, gaps,
reference difficulty) with categories Excellent (0.85 or more), Good (0.70),
Fair (0.50) and Poor. Weights and targets are documented in
`QualityScorer` and can be overridden in code. `--no-html-report`,
`--no-interactive-plots` and `--no-benchmarking` skip those steps.

### compare

```bash
baitUtils compare -r reference.fasta -o comparison/ \
  --sets "Design1:baits1.fasta" "Design2:baits2.fasta" \
  --significance-level 0.05 --multiple-comparison-correction fdr
```

Evaluates each set as `evaluate` does and writes `comparison.json`,
`comparison_matrix.csv`, `pairwise_comparisons.csv`,
`gap_overlap_analysis.txt`, `comparative_analysis_report.html` and plots.

Statistical tests are computed on observed data:

- Coverage distributions (two sets): Kolmogorov-Smirnov, Mann-Whitney U and
  Levene on mean depth per 1000 bp window
- Per-reference metrics (breadth, mean depth, gap count): Wilcoxon
  signed-rank for two sets, Friedman for more, each reference being one
  paired observation
- Best-hit identity per bait (two sets): Mann-Whitney U
- Gap sizes (two sets): Mann-Whitney U

P-values are adjusted within each family of tests with
`--multiple-comparison-correction` (bonferroni, holm, fdr). Tests without a
valid input, for example a single reference sequence, are reported as not
applicable rather than with a p-value. `--no-statistical-analysis` skips
this section.

## Workflow

```bash
baitUtils stats    -i baits.fasta -o stats/ --filter --mingc 35 --maxgc 65
baitUtils evaluate -i stats/filtered_sequences.fasta -r target.fasta -o eval/
baitUtils map      -i candidates.fasta -q target.fasta -o map/ --prefix cand
baitUtils fill     --alignments map/cand-mapping.psl --reference target.fasta --output add.txt
baitUtils compare  -r target.fasta -o cmp/ --sets "v1:baits.fasta" "v2:baits_v2.fasta"
```

## Testing

```bash
python -m pytest tests/
```

The suite includes end-to-end runs of every command on a small dataset with
hand-derived expected values; pblat is replaced by a script on PATH. Tests
that need bedtools, pblat or ViennaRNA are skipped when the tools are absent
and run in CI in a conda job.

## Citation

Sjodin, A. (2024). baitUtils: tools for oligo bait design, evaluation and
comparison. https://github.com/FOI-Bioinformatics/baitUtils

## License

MIT, see [LICENSE](LICENSE). Contributions are welcome; see
[CONTRIBUTING.md](CONTRIBUTING.md).
