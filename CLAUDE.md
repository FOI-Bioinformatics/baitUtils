# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

baitUtils is a Python CLI toolkit for analyzing and visualizing bait sequences used in in-solution hybridization. The package provides comprehensive statistical analysis, plotting, mapping, and coverage evaluation capabilities for bait sequence quality assessment.

## Architecture

The project follows a modular, well-factored architecture with clear separation of concerns:

### Main Entry Point
- **`__main__.py`**: CLI dispatcher using argparse subcommands

### Commands and Their Modules
- **`stats`** → `sequence_statistics.py` + `sequence_analysis.py` - Sequence statistics calculation
- **`plot`** → `statistical_plots.py` + `plotting_utils.py` - Statistical visualization generation  
- **`map`** → `sequence_mapping.py` + `mapping_utils.py` - Sequence mapping against references
- **`check`** → `coverage_evaluation.py` + `coverage_analysis.py` - Coverage evaluation and gap analysis
- **`fill`** → `gap_filling.py` + `gap_filling_algorithm.py` + `coverage_analysis.py` - Multi-pass gap filling
- **`evaluate`** → `evaluate.py` + supporting modules - Comprehensive oligo set evaluation
- **`compare`** → `compare.py` + comparative analysis modules - Multi-set comparison

### Core Utilities
- **`sequence_analysis.py`**: Core sequence analysis functions (GC content, Tm, MFE, entropy, complexity)
- **`coverage_analysis.py`**: hits to BED blocks (via the parser in `mapping_utils.py`), bedtools coverage, uncovered regions, N-split FASTA export for `check` and `fill`
- **`plotting_utils.py`**: Comprehensive plotting utilities with multiple chart types
- **`mapping_utils.py`**: pblat and minimap2 runners, the single alignment parser (`PSLHit`, `parse_psl`, `parse_paf`, `parse_alignments`, BLAT identity) and the per-bait hit table
- **`sequence_features.py`**: vectorized sequence features (GC, entropy, homopolymers, duplicated k-mers, window features and window coverage)
- **`bedtools_support.py`**: dependency check for the commands that need pybedtools and bedtools
- **`json_export.py`**: JSON writer, tool version record and argument record used by evaluate and compare
- **`config.py`**: JSON configuration file applied as argparse defaults (`--config` anywhere on the command line)
- **`cli_utils.py`**: `add_argument` with hidden aliases for the pre-0.6 option spellings, the shared logging options, argparse value types
- **`logging_utils.py`**: single logging configuration; the dispatcher configures it, command modules call `ensure_logging`
- **`templates/`**: CSS and JavaScript for the HTML reports
- **`gap_filling_algorithm.py`**: Multi-pass greedy selection algorithms

### Advanced Analysis Modules
- **`quality_scorer.py`**: Quality score in 0-1 with documented weights, targets and thresholds
- **`benchmark.py`**: Comparison of observed metrics with the design targets in `QualityScorer.DEFAULT_BENCHMARKS`
- **`reference_analyzer.py`**: Per-reference and per-window features with Spearman correlations against observed coverage
- **`design.py`**: the design system shared by every figure and both reports (palettes, chrome, number formats, matplotlib style, Plotly template, CSS variables, `plotly_layout`, `subplot_note`, `histogram_bins`)
- **`coverage_viz.py`**: static evaluate figures (coverage profile, overview, depth, gaps, per-reference, coverage map, quality score)
- **`interactive_plots.py`**: Interactive Plotly visualizations
- **`report_generator.py`**: HTML report generation

### Comparative Analysis
- **`comparative_analyzer.py`**: Multi-set comparison framework
- **`differential_analysis.py`**: Window-level coverage tests, paired per-reference tests and per-bait identity tests with multiple-comparison correction
- **`comparative_visualizations.py`**: Comparative plotting suite
- **`comparative_report_generator.py`**: Comparative HTML reports

## Development Commands

### Installation and Setup
Python 3.12 or later (`requires-python = ">=3.12"`; CI runs 3.12 and 3.13).
```bash
# External tools: pblat or minimap2 (map, evaluate, compare), bedtools +
# pybedtools (check, fill), ViennaRNA (optional, hairpin/dimer energies in stats)
conda create -n baitutils -c conda-forge -c bioconda python=3.12 pblat minimap2 bedtools pybedtools viennarna
conda activate baitutils

# Package in development mode with test and lint tools
pip install -e ".[dev]"
```
Packaging lives in `pyproject.toml` (`setup.py` is a shim); the version is
`baitUtils/_version.py`.

### Testing
```bash
# Run all tests (pytest; fixtures in tests/conftest.py build a small dataset
# with hand-derived expected values and a fake pblat on PATH). The suite
# uses pytest fixtures, so unittest discovery runs only part of it.
python -m pytest tests/

# Tests that need pblat, bedtools and ViennaRNA skip unless the tools are on PATH
# conda create -n baitutils-test -c conda-forge -c bioconda python=3.12 pybedtools bedtools pblat minimap2 viennarna
# PATH=$CONDA_PREFIX/bin:$PATH python -m pytest tests/

# Run specific test files
python -m pytest tests/test_evaluate.py tests/test_comparative_analysis.py

# Run tests with coverage (if coverage.py installed)
coverage run -m pytest tests/
coverage report
```

### Running the Application
```bash
# Run via installed command (after pip install -e .)
baitUtils --help
baitUtils stats --help
baitUtils evaluate --help

# Run via Python module
python -m baitUtils --help
python -m baitUtils stats -i sequences.fasta -o results/
```

### Code Quality
The checks CI enforces (the rest of ruff and flake8 runs with `--exit-zero`):
```bash
ruff check baitUtils/ --select=E9,F63,F7,F82,F401,F811,F841
flake8 baitUtils/ --select=E9,F63,F7,F82
```

### Releasing
1. Set `__version__` in `baitUtils/_version.py` and date the CHANGELOG section.
2. Merge to `main`, then either push a `v<version>` tag or run the Release
   workflow manually on `main` with `version: v<version>`.
3. The workflow checks that the tag matches `_version.py`, runs CI, builds
   the sdist and wheel and creates a **draft** GitHub release with them; a
   maintainer publishes it (publishing creates the tag for a manual run).
   The package is not on PyPI.

## Key Components

### Command Structure
Each command follows a consistent pattern:
- `add_arguments(parser)` - defines CLI arguments for the subcommand
- `main(args)` - entry point that receives parsed arguments
- Processor classes handle the main logic (e.g., `CoverageEvaluationProcessor`, `GapFillingProcessor`, `SequenceMappingProcessor`)
- Utility classes handle specific functionality

### Dependencies
- **Scientific Python stack**: numpy, pandas, matplotlib, seaborn, scikit-learn
- **Bioinformatics**: biopython for sequence handling
- **Interactive plotting**: plotly for enhanced visualizations
- **Optional external tools**: 
  - ViennaRNA (for hairpin and self-dimer free energies)
  - pblat or minimap2 (for sequence mapping, `--mapper`)
  - bedtools with pybedtools (for `check` and `fill`)

### Input/Output Patterns
- **Primary input**: FASTA/FASTA.GZ files containing sequences
- **Output formats**: TSV statistics, various plot formats (PNG/PDF/SVG), filtered FASTA, HTML reports
- **Directory structure**: Results organized in structured output directories

## Module Architecture Details

### Sequence Analysis (`sequence_analysis.py`)
Core sequence analysis functionality:
- GC content calculation using BioPython
- Melting temperature via BioPython's MeltingTemp module
- Hairpin and self-dimer free energies with ViennaRNA DNA parameters at the hybridization temperature (optional)
- Shannon entropy and sequence complexity metrics
- Homopolymer run analysis
- Masked base counting

### Statistics Processing (`sequence_statistics.py`)
Main statistics calculation workflow:
- `SequenceStatsCalculator`: Orchestrates analysis with parallel processing support
- `SequenceFilter`: Configurable filtering based on multiple criteria
- Supports both sequential and parallel processing modes
- Handles compressed FASTA files

### Plotting System (`plotting_utils.py`, `statistical_plots.py`)
Comprehensive visualization system:
- Multiple plot types: histograms, boxplots, scatterplots, violin plots, PCA
- Color coding support for categorical data
- Batch generation of pairwise combination plots
- Statistical plotting with seaborn integration

### Mapping System (`mapping_utils.py`, `sequence_mapping.py`) 
Sequence mapping workflow:
- pblat or minimap2 (`--mapper`), with BLAT identity on both so filters are comparable
- One parser for PSL and PAF (`parse_alignments`) used by every command
- Results categorization into mapped/unmapped sequences
- FASTA output generation for different categories

### Coverage Analysis (`coverage_analysis.py`, `coverage_evaluation.py`)
Coverage evaluation system:
- PSL to BED conversion with filtering
- Coverage calculation using bedtools integration
- Uncovered region identification and analysis
- FASTA export of uncovered regions with N-splitting

### Gap Filling (`gap_filling_algorithm.py`, `gap_filling.py`)
Multi-pass optimization system:
- Greedy selection algorithm with multiple scoring criteria
- Coverage pattern analysis for difficult regions
- Spacing constraints for oligo selection
- Multi-pass iteration with convergence detection

## Testing Strategy

The suite runs under pytest (tests/conftest.py provides a small dataset with hand-derived expected values and a fake pblat):

### Test Coverage Areas
- **Core algorithms**: Sequence analysis functions, statistics calculations
- **File I/O operations**: FASTA reading, results writing with mocking
- **Command workflows**: End-to-end command testing
- **Edge cases**: Error handling, malformed inputs, empty datasets
- **Integration tests**: Multi-module workflows

### Fixtures and external tools
- `tests/conftest.py` builds the dataset and puts fake `pblat` and `minimap2`
  scripts on PATH that return precomputed hits for the baits in the query FASTA
- The `run_cli` fixture runs commands end to end through the dispatcher
- Tests that need bedtools, pblat or ViennaRNA skip when the tool is absent;
  CI's `external-tools` job runs them in a conda environment
- Some older unit tests mock matplotlib's `savefig`; keep `plt.savefig` or
  `Figure.savefig` calls where those tests expect them

## Enhanced Coverage Evaluation (evaluate command)

The `evaluate` command provides comprehensive oligo set coverage analysis:

### Usage Examples
```bash
# Basic usage
baitUtils evaluate -i oligos.fasta -r reference.fasta -o coverage_report/

# Advanced analysis with custom parameters
baitUtils evaluate -i oligos.fasta -r reference.fasta -o results/ \
  --min-identity 95 --target-coverage 10 --threads 4 \
  --plot-format pdf --offline-plots
```

### Architecture
- **`evaluate.py`**: Main orchestrator integrating all analysis components
- **`coverage_stats.py`**: Coverage statistics computation
- **`quality_scorer.py`**: Quality score in 0-1 with Excellent/Good/Fair/Poor categories
- **`reference_analyzer.py`**: Per-reference and per-window features with Spearman correlations against observed coverage
- **`gap_analysis.py`**: Gap characterization with sequence correlation

## Comparative Analysis (compare command)

### Multi-Set Comparison
```bash
# Compare multiple oligo sets
baitUtils compare -r reference.fasta -o comparison_report/ \
  --sets "Design1:oligos1.fasta" "Design2:oligos2.fasta" "Design3:oligos3.fasta"

# With statistical analysis
baitUtils compare -r reference.fasta -o results/ \
  --sets "Set1:design1.fasta" "Set2:design2.fasta" \
  --significance-level 0.01 --multiple-comparison-correction holm
```

### Comparative Framework
- Each set is evaluated as in `evaluate` (same options)
- Statistical tests on observed data: window depths (KS, Mann-Whitney, Levene),
  per-reference metrics (Wilcoxon or Friedman), bait identity and gap sizes,
  adjusted with `--multiple-comparison-correction`; tests without valid input
  are reported as not applicable
- Ranking by a composite score in 0-1 (quality 40 %, breadth 25 %,
  mapping efficiency 20 %, gap term 15 %)
- HTML report plus `comparison.json`, `comparison_matrix.csv`, `pairwise_comparisons.csv`

## Reports and figures

- All colours, fonts, number formats and chart chrome come from `design.py`;
  do not hard-code colours or palettes in plotting code.
- Categorical colours follow the entity in a fixed order (`design.series_color`,
  `comparative_visualizations.set_colors`); status colours are only for the
  quality category; `design.CHROME["gap_fill"]` only for gaps.
- Matplotlib figures are drawn inside `@design.styled`, which applies the
  style in an `rc_context`. Never set global rcParams: a leaked
  `patch.linewidth = 0` broke seaborn plots in `plot` under matplotlib 3.10.
- Plotly figures use `design.plotly_layout` (template, legend band) and
  `design.subplot_note` for empty panels; report div ids are prefixed `fig-`
  so they never collide with section ids.
- Both reports share `templates/report.css` and `templates/report.js`.
  See `docs/report_design_audit_2026-10.md`.

## Error Handling and Robustness

### Graceful Degradation
- Optional dependencies handled with try/except imports
- Fallback methods when external tools unavailable
- Clear error messages for missing dependencies
- Continuation with reduced functionality when possible

### Input Validation
- argparse value types in `cli_utils.py` (`percent`, `positive_int`, `non_negative_int`)
- `stats`, `plot`, `evaluate` and `compare` check that their input files exist
  before starting; `map`, `check` and `fill` fail when the mapper or bedtools reads them

## Performance Considerations

- `stats` uses a process pool (`--threads`, default all CPUs but one)
- Coverage is held as one numpy depth array per reference; sequence
  features are vectorized (`sequence_features.py`) and linear in reference length
- `evaluate` and `compare` read the reference FASTA once and share it
  (`CoverageAnalyzer.sequences()`)
- Coverage is exported run-length encoded (`coverage.bedgraph`); the
  per-base CSV only with `--per-position-coverage`

## Integration Workflow

The complete baitUtils workflow supports iterative oligo design:

1. **Initial Analysis** (`evaluate`): Assess current oligo set quality
2. **Gap Identification** (`check`): Identify coverage gaps and problematic regions  
3. **Gap Filling** (`fill`): Generate improved oligo selections
4. **Comparison** (`compare`): Evaluate multiple design iterations
5. **Visualization** (`plot`): Generate publication-quality figures
6. **Statistics** (`stats`): Detailed sequence property analysis

This modular architecture enables flexible workflows adapted to specific research needs while maintaining code quality and testability.