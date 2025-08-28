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
- **`check`** → `coverage_evaluation.py` + `coverage_checking.py` - Coverage evaluation and gap analysis
- **`fill`** → `gap_filling.py` + `gap_filling_algorithm.py` + `coverage_analysis.py` - Multi-pass gap filling
- **`evaluate`** → `evaluate.py` + supporting modules - Comprehensive oligo set evaluation
- **`compare`** → `compare.py` + comparative analysis modules - Multi-set comparison

### Core Utilities
- **`sequence_analysis.py`**: Core sequence analysis functions (GC content, Tm, MFE, entropy, complexity)
- **`coverage_analysis.py`**: Coverage calculation utilities, PSL parsing, BED operations
- **`plotting_utils.py`**: Comprehensive plotting utilities with multiple chart types
- **`mapping_utils.py`**: Sequence mapping utilities with pblat integration
- **`coverage_checking.py`**: Coverage checking and uncovered region analysis
- **`gap_filling_algorithm.py`**: Multi-pass greedy selection algorithms

### Advanced Analysis Modules
- **`quality_scorer.py`**: Quantitative quality assessment system
- **`benchmark.py`**: Performance benchmarking against theoretical optimal
- **`reference_analyzer.py`**: Reference sequence complexity analysis
- **`interactive_plots.py`**: Interactive Plotly visualizations
- **`report_generator.py`**: HTML report generation

### Comparative Analysis
- **`comparative_analyzer.py`**: Multi-set comparison framework
- **`differential_analysis.py`**: Statistical testing for set comparisons
- **`comparative_visualizations.py`**: Comparative plotting suite
- **`comparative_report_generator.py`**: Comparative HTML reports

## Development Commands

### Installation and Setup
```bash
# Create conda environment with dependencies
conda create -n baitutils_env python=3.12 numpy pandas matplotlib-base seaborn scikit-learn biopython
conda activate baitutils_env

# Install optional external tools
conda install -c bioconda pblat bedtools
# Note: ViennaRNA installation is optional for MFE calculations
conda install -c bioconda viennarna

# Install package in development mode
pip install -e .
```

### Testing
```bash
# Run all tests
python -m unittest discover tests/

# Run specific test files
python -m unittest tests.test_evaluate
python -m unittest tests.test_comparative_analysis

# Run tests with verbose output
python -m unittest discover tests/ -v

# Run tests with coverage (if coverage.py installed)
coverage run -m unittest discover tests/
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
```bash
# Run linting (if installed)
ruff check baitUtils/
flake8 baitUtils/

# Format code (if installed)  
black baitUtils/
```

## Key Components

### Command Structure
Each command follows a consistent pattern:
- `add_arguments(parser)` - defines CLI arguments for the subcommand
- `main(args)` - entry point that receives parsed arguments
- Processor classes handle the main logic (e.g., `SequenceStatsProcessor`)
- Utility classes handle specific functionality

### Dependencies
- **Scientific Python stack**: numpy, pandas, matplotlib, seaborn, scikit-learn
- **Bioinformatics**: biopython for sequence handling
- **Interactive plotting**: plotly for enhanced visualizations
- **Optional external tools**: 
  - ViennaRNA (for RNA folding/MFE calculations)
  - pblat (for sequence mapping)
  - bedtools (for genomic interval operations)

### Input/Output Patterns
- **Primary input**: FASTA/FASTA.GZ files containing sequences
- **Output formats**: TSV statistics, various plot formats (PNG/PDF/SVG), filtered FASTA, HTML reports
- **Directory structure**: Results organized in structured output directories

## Module Architecture Details

### Sequence Analysis (`sequence_analysis.py`)
Core sequence analysis functionality:
- GC content calculation using BioPython
- Melting temperature via BioPython's MeltingTemp module
- Minimum Free Energy (MFE) using ViennaRNA (optional)
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
- pblat integration with configurable parameters
- PSL file parsing with filtering capabilities
- Results categorization into mapped/unmapped sequences
- FASTA output generation for different categories

### Coverage Analysis (`coverage_checking.py`, `coverage_evaluation.py`)
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

The codebase uses unittest with comprehensive mocking for external dependencies:

### Test Coverage Areas
- **Core algorithms**: Sequence analysis functions, statistics calculations
- **File I/O operations**: FASTA reading, results writing with mocking
- **Command workflows**: End-to-end command testing
- **Edge cases**: Error handling, malformed inputs, empty datasets
- **Integration tests**: Multi-module workflows

### Mock Strategy
- File operations mocked to avoid filesystem dependencies
- External tool calls (pblat, bedtools) mocked in tests
- Optional dependencies handled gracefully with availability checks

## Enhanced Coverage Evaluation (evaluate command)

The `evaluate` command provides comprehensive oligo set coverage analysis:

### Usage Examples
```bash
# Basic usage
baitUtils evaluate -i oligos.fasta -r reference.fasta -o coverage_report/

# Advanced analysis with custom parameters
baitUtils evaluate -i oligos.fasta -r reference.fasta -o results/ \
  --min-identity 95 --target-coverage 10 --threads 4 \
  --plot-format pdf --enable-html-report
```

### Architecture
- **`evaluate.py`**: Main orchestrator integrating all analysis components
- **`coverage_stats.py`**: Coverage statistics computation
- **`quality_scorer.py`**: Quality assessment with A-F grading
- **`reference_analyzer.py`**: Reference sequence complexity analysis
- **`gap_analysis.py`**: Gap characterization with sequence correlation

## Phase 3: Comparative Analysis (compare command)

### Multi-Set Comparison
```bash
# Compare multiple oligo sets
baitUtils compare -r reference.fasta -o comparison_report/ \
  --sets "Design1:oligos1.fasta" "Design2:oligos2.fasta" "Design3:oligos3.fasta"

# With statistical analysis
baitUtils compare -r reference.fasta -o results/ \
  --sets "Set1:design1.fasta" "Set2:design2.fasta" \
  --enable-statistical-analysis --significance-level 0.01
```

### Comparative Framework
- **Multi-set analysis**: Simultaneous evaluation of multiple designs
- **Statistical testing**: Rigorous comparison with significance testing
- **Performance ranking**: Composite scoring and recommendation system
- **Interactive reporting**: Comprehensive HTML reports with embedded analysis

## Error Handling and Robustness

### Graceful Degradation
- Optional dependencies handled with try/except imports
- Fallback methods when external tools unavailable
- Clear error messages for missing dependencies
- Continuation with reduced functionality when possible

### Input Validation
- FASTA file format validation
- Parameter range checking
- File existence verification
- Memory usage considerations for large datasets

## Performance Considerations

### Parallel Processing
- Multi-threaded sequence analysis where beneficial
- Configurable process counts for compute-intensive operations
- Memory-efficient streaming for large files

### Caching and Optimization
- Intermediate file caching for complex workflows
- Optimized data structures for coverage calculations
- Efficient algorithms for gap filling optimization

## Integration Workflow

The complete baitUtils workflow supports iterative oligo design:

1. **Initial Analysis** (`evaluate`): Assess current oligo set quality
2. **Gap Identification** (`check`): Identify coverage gaps and problematic regions  
3. **Gap Filling** (`fill`): Generate improved oligo selections
4. **Comparison** (`compare`): Evaluate multiple design iterations
5. **Visualization** (`plot`): Generate publication-quality figures
6. **Statistics** (`stats`): Detailed sequence property analysis

This modular architecture enables flexible workflows adapted to specific research needs while maintaining code quality and testability.