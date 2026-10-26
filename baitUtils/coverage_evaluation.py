#!/usr/bin/env python3

"""
coverage_evaluation.py

Main command for evaluating oligo coverage and reporting uncovered regions.
Refactored from check.py for better organization and maintainability.
"""

import argparse
import logging
import sys
import tempfile
from pathlib import Path

from baitUtils._version import __version__
from baitUtils.cli_utils import (add_argument, add_logging_arguments, logging_kwargs, percent,
                                 positive_int, non_negative_int)
from baitUtils.logging_utils import ensure_logging
from baitUtils.bedtools_support import require_bedtools
from baitUtils.coverage_analysis import CoverageChecker, PSLParser, ForcedOligoHandler


class CoverageEvaluationProcessor:
    """Processes coverage evaluation commands with all options."""
    
    def __init__(self):
        """Initialize coverage evaluation processor."""
        self.coverage_checker = CoverageChecker()
        self.psl_converter = PSLParser()
        self.forced_filter = ForcedOligoHandler()
    
    def process_command(self, args) -> None:
        """
        Process the coverage evaluation command with given arguments.
        
        Args:
            args: Parsed command-line arguments
        """
        # Set up logging
        ensure_logging(**logging_kwargs(args))
        require_bedtools()

        import pybedtools
        forced_oligos = self.forced_filter.read_forced_oligos(args.forced_oligos)
        
        if args.temp_dir:
            Path(args.temp_dir).mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="baitUtils_check_", dir=args.temp_dir) as tmp:
            work_dir = Path(tmp)
            previous_tempdir = tempfile.tempdir
            pybedtools.set_tempdir(str(work_dir))  # sets tempfile.tempdir globally
            try:
                bed = self.psl_converter.parse_psl_to_bed(
                    args.alignments, args.min_length, args.min_identity, work_dir, "temp_check.bed",
                    strand=args.strand
                )
                if bed.count() == 0:
                    logging.error("No intervals found in PSL after filtering. Exiting.")
                    sys.exit(1)
                
                self.coverage_checker.check_coverage(
                    bed,
                    forced_oligos,
                    args.min_coverage,
                    args.coverage_out,
                    args.longest_uncovered_out,
                    work_dir,
                    args.uncovered_length_cutoff,
                    args
                )
            finally:
                pybedtools.cleanup(remove_all=True)
                tempfile.tempdir = previous_tempdir
    


def add_coverage_arguments(parser: argparse.ArgumentParser, min_coverage_help: str) -> None:
    """Options shared by check and fill: alignment input, filters, uncovered-region outputs."""
    # Input files
    add_argument(parser, "--alignments", dest="alignments", type=Path, required=True, aliases=("--psl",),
                 help="Alignment file: PSL from pblat or PAF from minimap2 -c (by extension)")
    add_argument(parser, "--forced-oligos", dest="forced_oligos", type=Path, aliases=("--forced_oligos",),
                 help="File with oligo IDs that must be included")
    parser.add_argument("--reference", type=Path, required=True,
                        help="Reference FASTA; used for reference sizes, sequence-based scoring "
                             "and export of uncovered regions")
    
    # Coverage and filtering parameters
    add_argument(parser, "--min-coverage", dest="min_coverage", type=float, default=1.0,
                 aliases=("--min_coverage",), help=min_coverage_help)
    add_argument(parser, "--min-length", dest="min_length", type=non_negative_int, default=100,
                 aliases=("--min_length",),
                 help="Minimum aligned length in bases for a hit to count (default: 100)")
    add_argument(parser, "--min-identity", dest="min_identity", type=percent, default=90.0,
                 aliases=("--min_similarity", "--min-similarity"),
                 help="Minimum BLAT percent identity for a hit to count (default: 90)")
    parser.add_argument("--strand", choices=["both", "plus", "minus"], default="both",
                        help="Count hits on both strands, or only plus or minus strand hits (default: both)")
    add_argument(parser, "--uncovered-length-cutoff", dest="uncovered_length_cutoff", type=non_negative_int,
                 default=0, aliases=("--uncovered_length_cutoff",),
                 help="Report uncovered stretches of at least this length (default: 0)")
    
    # Output files
    add_argument(parser, "--coverage-out", dest="coverage_out", type=Path, aliases=("--coverage_out",),
                 help="Output file for coverage data")
    add_argument(parser, "--longest-uncovered-out", dest="longest_uncovered_out", type=Path,
                 aliases=("--longest_uncovered_out",),
                 help="Output file for uncovered regions, longest first")
    add_argument(parser, "--uncovered-fasta", dest="uncovered_fasta", type=Path, aliases=("--uncovered_fasta",),
                 help="Output FASTA file for uncovered regions")
    add_argument(parser, "--n-split-fasta", dest="n_split_fasta", type=Path, aliases=("--n_split_fasta",),
                 help="Output FASTA file for uncovered regions split at N runs")
    
    # FASTA export parameters
    add_argument(parser, "--extend-region", dest="extend_region", type=non_negative_int, default=0,
                 aliases=("--extend_region",),
                 help="Number of bases to extend uncovered regions by on each side (default: 0)")
    add_argument(parser, "--min-oligo-length", dest="min_oligo_length", type=positive_int, default=120,
                 aliases=("--min_oligo_length",),
                 help="Minimum piece length after N-splitting (default: 120)")
    
    # System parameters
    add_argument(parser, "--temp-dir", dest="temp_dir", type=Path, aliases=("--temp_dir",),
                 help="Parent directory for the run's temporary directory (default: system temp)")


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add command-line arguments for coverage evaluation."""
    parser.add_argument('--version', action='version', version=f"baitUtils check {__version__}")
    add_coverage_arguments(parser, "Minimum depth for a base to count as covered (default: 1)")
    add_logging_arguments(parser)

def main(args) -> None:
    """Main function for coverage evaluation command."""
    processor = CoverageEvaluationProcessor()
    try:
        processor.process_command(args)
    except Exception as e:
        logging.error(f"Coverage evaluation failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Evaluate coverage of oligo sets and report uncovered regions.'
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)