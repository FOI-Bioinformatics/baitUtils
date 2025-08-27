# baitUtils

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

**Comprehensive toolkit for oligo/bait design, evaluation, and comparative analysis with advanced coverage assessment and interactive reporting.**

baitUtils is a powerful command-line toolkit designed for researchers working with oligonucleotide baits in targeted sequencing applications. It provides end-to-end functionality from basic sequence statistics to sophisticated comparative analysis with publication-quality visualizations and interactive reports.

## 🚀 Key Features

### 📊 **Comprehensive Analysis Pipeline**
- **Individual Evaluation**: Detailed coverage analysis with quality scoring (Phase 1 & 2)
- **Comparative Analysis**: Statistical comparison of multiple oligo sets (Phase 3)
- **Interactive Reporting**: Professional HTML reports with embedded visualizations
- **Publication-Quality Plots**: High-resolution static and interactive visualizations

### 🔬 **Advanced Metrics & Statistics**
- Coverage breadth and depth analysis with uniformity assessment
- Gap detection and characterization with sequence feature correlation
- Quality scoring system with benchmarking against theoretical optimal
- Statistical significance testing with multiple comparison corrections

### 🎯 **Intelligent Recommendations**
- Automated best performer identification
- Gap-filling optimization strategies
- Performance improvement suggestions
- Design validation and comparison

## 📦 Installation

### Prerequisites
- **Python 3.11+** (recommended: Python 3.12)
- **pblat**: Required for sequence mapping ([installation guide](https://github.com/icebert/pblat))
- **bedtools**: Required for coverage analysis ([installation guide](https://bedtools.readthedocs.io/))

### Conda Environment (Recommended)
```bash
# Create conda environment with dependencies
conda create -n baitutils python=3.12
conda activate baitutils

# Install external tools
conda install -c bioconda pblat bedtools

# Install baitUtils
pip install baitutils
```

### From Source
```bash
git clone https://github.com/FOI-Bioinformatics/baitUtils.git
cd baitUtils
pip install -e .
```

### Dependencies
- numpy>=1.21.0
- pandas>=1.3.0
- matplotlib>=3.4.0
- seaborn>=0.11.0
- scikit-learn>=0.24.0
- biopython>=1.78
- plotly>=5.0.0
- scipy>=1.7.0
- pybedtools>=0.9.0

## 🎮 Commands Overview

baitUtils provides a comprehensive suite of commands for different analysis needs:

```bash
baitUtils <command> [options]
```

| Command | Purpose | Key Features |
|---------|---------|-------------|
| `stats` | Basic sequence statistics | GC content, Tm, MFE, filtering |
| `plot` | Statistical visualizations | Histograms, boxplots, PCA, heatmaps |
| `map` | Sequence mapping | pblat-based mapping with filtering |
| `check` | Coverage evaluation | Basic coverage assessment |
| `fill` | Gap optimization | Iterative gap-filling strategies |
| `evaluate` | **Comprehensive analysis** | **Complete coverage evaluation pipeline** |
| `compare` | **Multi-set comparison** | **Statistical comparison of multiple designs** |

## 📋 Quick Start

### 1. Basic Sequence Statistics
```bash
# Calculate statistics and filter sequences
baitUtils stats -i oligos.fasta -o stats_output \
  --length 120 --mingc 40 --maxgc 60 --filter
```

### 2. Generate Visualizations
```bash
# Create plots from statistics
baitUtils plot -i stats_output/filtered-params.txt -o plots \
  --columns GC% Tm MFE --plot_type histogram boxplot scatterplot
```

### 3. Comprehensive Evaluation (Recommended)
```bash
# Complete coverage analysis with interactive reporting
baitUtils evaluate -i oligos.fasta -r reference.fasta -o evaluation_report \
  --min-identity 95 --target-coverage 10 --threads 4
```

### 4. Compare Multiple Designs
```bash
# Statistical comparison of multiple oligo sets
baitUtils compare -r reference.fasta -o comparison_report \
  --sets "Design1:oligos1.fasta" "Design2:oligos2.fasta" "Design3:oligos3.fasta" \
  --enable-statistical-analysis --significance-level 0.01
```

## 🔬 Detailed Command Guide

### `evaluate` - Comprehensive Coverage Analysis

**Purpose**: Complete evaluation of a single oligo set with quality assessment, gap analysis, and interactive reporting.

```bash
baitUtils evaluate -i oligos.fasta -r reference.fasta -o results/

# Advanced options
baitUtils evaluate -i oligos.fasta -r reference.fasta -o results/ \
  --min-identity 95.0 \
  --target-coverage 10 \
  --min-gap-size 100 \
  --threads 8 \
  --plot-format pdf \
  --enable-html-report \
  --enable-interactive-plots \
  --enable-benchmarking
```

**Output Structure**:
```
results/
├── interactive_report.html          # Main comprehensive report
├── coverage_statistics.txt          # Detailed numeric results  
├── gap_analysis.txt                 # Gap characteristics
├── recommendations.txt              # Improvement suggestions
├── benchmark_analysis.txt           # Performance benchmarking
├── plots/                           # Static visualizations
├── interactive_plots/               # Interactive visualizations
└── data/                           # Raw analysis data
```

**Key Metrics**:
- **Coverage Breadth**: % of reference sequence covered
- **Coverage Depth**: Distribution of oligo overlap
- **Gap Analysis**: Uncovered region characteristics
- **Quality Score**: 0-10 scale with letter grades (A-F)
- **Mapping Efficiency**: % of oligos successfully mapped
- **Uniformity**: Coverage distribution evenness (Gini coefficient)

### `compare` - Multi-Set Comparative Analysis

**Purpose**: Statistical comparison of multiple oligo set designs to identify the best performer.

```bash
baitUtils compare -r reference.fasta -o comparison_report/ \
  --sets "Original:design1.fasta" "Optimized:design2.fasta" "Alternative:design3.fasta"

# With statistical analysis
baitUtils compare -r reference.fasta -o comparison_report/ \
  --sets "Method_A:setA.fasta" "Method_B:setB.fasta" \
  --enable-statistical-analysis \
  --significance-level 0.01 \
  --multiple-comparison-correction fdr
```

**Statistical Tests**:
- **Kolmogorov-Smirnov**: Coverage distribution differences
- **Mann-Whitney U**: Median coverage comparisons
- **Levene's Test**: Coverage variability assessment
- **t-test/ANOVA**: Quality metric comparisons
- **Effect Sizes**: Cohen's d, eta-squared for practical significance

**Output Features**:
- Executive summary with best performer identification
- Detailed comparison tables and statistical results
- Gap overlap analysis (unique vs shared gaps)
- Performance ranking with composite scoring
- Interactive comparative visualizations

### Legacy Commands

#### `stats` - Basic Sequence Statistics
```bash
baitUtils stats -i sequences.fasta -o output_dir \
  --length 120 --mingc 40 --maxgc 60 --filter
```

#### `plot` - Statistical Visualizations  
```bash
baitUtils plot -i params.txt -o plots_dir \
  --columns GC% Tm --plot_type histogram scatterplot
```

#### `map` - Sequence Mapping
```bash
baitUtils map -i baits.fasta -q reference.fasta -o mapping_results \
  --threads 4 --minIdentity 90 --filterIdentity 95
```

## 📊 Analysis Workflow Examples

### Complete Single-Set Analysis
```bash
# Step 1: Comprehensive evaluation
baitUtils evaluate -i my_oligos.fasta -r target_genome.fasta -o evaluation/

# Step 2: Review interactive report
open evaluation/interactive_report.html

# Step 3: Optimize based on recommendations
# (Use gap_regions.bed for targeted improvement)
```

### Multi-Set Design Selection
```bash
# Compare multiple design approaches
baitUtils compare -r reference.fasta -o comparison/ \
  --sets "Conservative:conservative_design.fasta" \
         "Aggressive:aggressive_design.fasta" \
         "Balanced:balanced_design.fasta" \
  --enable-statistical-analysis

# Results will identify the best performer with statistical significance
```

### Publication-Quality Analysis
```bash
# Generate high-resolution plots for publications
baitUtils compare -r genome.fasta -o publication_analysis/ \
  --sets "Method1:approach1.fasta" "Method2:approach2.fasta" \
  --plot-format pdf --plot-dpi 600 \
  --significance-level 0.001 \
  --multiple-comparison-correction bonferroni
```

## 📈 Understanding Output

### Quality Scores & Grades
- **A (9.0-10.0)**: Excellent performance, ready for use
- **B (8.0-8.9)**: Good performance, minor optimizations possible  
- **C (7.0-7.9)**: Fair performance, consider improvements
- **D (6.0-6.9)**: Poor performance, significant optimization needed
- **F (<6.0)**: Unacceptable performance, major redesign required

### Key Performance Indicators
- **Coverage Breadth >90%**: Excellent target coverage
- **Gap Count <50**: Acceptable coverage continuity
- **Mapping Efficiency >85%**: Good oligo design specificity
- **Gini Coefficient <0.3**: Uniform coverage distribution

### Statistical Significance Levels
- **\*\*\* (p≤0.001)**: Highly significant difference
- **\*\* (p≤0.01)**: Very significant difference
- **\* (p≤0.05)**: Significant difference
- **ns (p>0.05)**: Not significant

## 🔧 Advanced Configuration

### Custom Quality Scoring Weights
```python
# Default weights can be customized in quality_scorer.py
weights = {
    'coverage': 0.30,      # Coverage breadth importance
    'depth': 0.25,         # Depth uniformity importance  
    'gaps': 0.25,          # Gap minimization importance
    'mapping': 0.20        # Mapping efficiency importance
}
```

### Analysis Parameters
- **--min-identity**: Minimum mapping identity (default: 90.0%)
- **--target-coverage**: Target depth for analysis (default: 10x)
- **--min-gap-size**: Minimum gap size to report (default: 100bp)
- **--reference-analysis-window**: Window for sequence analysis (default: 1000bp)

## 🐛 Troubleshooting

### Common Issues

**1. "pblat not found in PATH"**
```bash
# Install pblat via conda
conda install -c bioconda pblat
```

**2. "bedtools not found"**
```bash
# Install bedtools via conda  
conda install -c bioconda bedtools
```

**3. Memory issues with large datasets**
```bash
# Reduce thread count and increase system memory
baitUtils evaluate ... --threads 2
```

**4. Interactive plots not displaying**
- Ensure plotly>=5.0.0 is installed
- Use a modern web browser to view HTML reports

### Performance Optimization
- Use `--threads` parameter for parallel processing
- Set appropriate `--min-identity` thresholds to reduce computation
- Use `--quiet` flag to reduce logging overhead

## 📝 Citation

If you use baitUtils in your research, please cite:

```
Sjödin, A. (2024). baitUtils: Comprehensive toolkit for oligo/bait design, 
evaluation, and comparative analysis. GitHub repository: 
https://github.com/FOI-Bioinformatics/baitUtils
```

## 🤝 Contributing

Contributions are welcome! Please see our [contribution guidelines](CONTRIBUTING.md) for details.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🆘 Support

- **Issues**: [GitHub Issues](https://github.com/FOI-Bioinformatics/baitUtils/issues)
- **Discussions**: [GitHub Discussions](https://github.com/FOI-Bioinformatics/baitUtils/discussions)
- **Email**: andreas.sjodin@gmail.com

## 🙏 Acknowledgments

- Built with modern Python scientific stack (NumPy, Pandas, Matplotlib, Plotly)
- Bioinformatics tools integration (BioPython, pblat, bedtools)
- Interactive reporting powered by Plotly and HTML5
- Statistical analysis using SciPy and scikit-learn

---

**baitUtils** - *Empowering precise oligo design through comprehensive analysis* 🧬