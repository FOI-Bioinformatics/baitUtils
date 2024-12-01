
# baitUtils

## Overview

baitUtils is a comprehensive toolkit for the analysis and visualization of bait sequences used in in-solution hybridization. It provides tools for generating bait quality statistics and visualizations.

## Installation

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html): Please ensure you have conda installed to manage dependencies.
- Install the required packages and dependencies with:

```bash
conda create -n baitutils_env numpy pandas matplotlib-base seaborn scikit-learn biopython viennarna pblat
conda activate baitutils_env
```

### Clone the Repository

```bash
git clone https://github.com/FOI-Bioinformatics/baitUtils.git
cd baitUtils
pip install .
```

## Usage

### General Structure

baitUtils offers two main functionalities: `stats` and `plot`, both accessible as subcommands of the primary script.

```bash
baitUtils [command] [options]
```

### Commands

- `stats`: Calculate quality statistics of bait sequences.
- `plot`: Generate plots based on bait sequence statistics.

## Commands and Options

### baitUtils stats

Calculates statistics on bait sequences and filters them based on user-defined criteria.

#### Example Usage

```bash
baitUtils stats -i probes.fasta.gz -o results --length 120 --mingc 40 --maxgc 60 --filter
```

#### Options

- `-i, --input`: Path to the input FASTA or FASTA.GZ file.
- `-o, --outdir`: Output directory for results.
- `--length`: Requested bait length (default is 120).
- `--mingc`: Minimum GC content percentage.
- `--maxgc`: Maximum GC content percentage.
- `--filter`: Save filtered FASTA output.

### baitUtils plot

Generates plots based on the bait sequence statistics file.

#### Example Usage

```bash
baitUtils plot -i results/filtered-params.txt -o plots --columns GC% Tm MFE --plot_type histogram boxplot scatterplot
```

#### Options

- `-i, --input`: Path to the parameters file.
- `-o, --outdir`: Output directory for plots.
- `--columns`: List of columns to include in plots.
- `--plot_type`: Types of plots to generate.
- `--color`: Column to use for coloring plots.

### baitUtils map

Maps bait sequences against a reference genome using pblat and filters mappings based on identity percentage.

#### Example Usage

```bash
baitUtils map -i baits.fa -q genome.fa -o mapping_results --outdir mappings --threads 4 --minIdentity 90 --filterIdentity 95 --fasta-output both
```

#### Options

- `-i, --input`: Path to the input baits FASTA file.
- `-q, --query`: Path to the reference genome FASTA file to map against.
- `-o, --outprefix`: Prefix for the output files (default is out).
- `-Z, --outdir`: Output directory path (default is .).
- `--mapper`: Mapping tool to use (pblat is currently supported).
- `-X, --threads`: Number of threads to use for mapping (default is 1).
- `--minMatch`: Minimum number of tile matches (default is 2 for nucleotide sequences).
- `--minScore`: Minimum score for alignments (default is 30).
- `--minIdentity`: Minimum sequence identity percentage for mappings (default is 90).
- `--filterIdentity`: Filter mappings with identity percentage less than this value; must be ≥ minIdentity (default is 90).
- `--fasta-output`: Choose which probes to include in the FASTA output file: mapped, unmapped, both, or none (default is mapped).
- `-l, --log`: Enable detailed logging for execution insights.

## Examples

### Stats Example

To calculate and filter baits based on length, GC content, and other quality metrics:

```bash
baitUtils stats -i probes.fasta.gz -o stats_output --length 120 --mingc 40 --maxgc 60 --filter
```

### Plot Example

To generate histograms, boxplots, and scatterplots for GC content and melting temperature:

```bash
baitUtils plot -i stats_output/filtered-params.txt -o plots_output --columns GC% Tm --plot_type histogram scatterplot --color Kept
```

### Map Example

To map baits against a reference genome and output both mapped and unmapped probes:

```bash
baitUtils map -i baits.fa -q genome.fa -o mapping_results --outdir mappings --threads 4 --minIdentity 90 --filterIdentity 95 --fasta-output both
```


## License

MIT License. See `LICENSE` file for details.
