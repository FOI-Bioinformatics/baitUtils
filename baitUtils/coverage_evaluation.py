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
        self._setup_logging(args.log_level)
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
                    args.psl, args.min_length, args.min_similarity, work_dir, "temp_check.bed",
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
    
    def _setup_logging(self, log_level: str = "INFO") -> None:
        """Configure logging."""
        ensure_logging(level=log_level)

def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add command-line arguments for coverage evaluation."""
    parser.add_argument('--version', action='version', 
                       version=f"baitUtils coverage-evaluation {__version__}")
    
    # Input files
    parser.add_argument("--alignments", "--psl", dest="psl", type=Path, required=True,
                       help="Alignment file: PSL from pblat or PAF from minimap2 -c (by extension)")
    parser.add_argument("--forced_oligos", type=Path,
                       help="File with oligo IDs that must be included in coverage check")
    parser.add_argument("--reference", type=Path, required=True,
                       help="Reference FASTA file for exporting uncovered regions")
    
    # Coverage parameters
    parser.add_argument("--min_coverage", type=float, default=10.0,
                       help="Minimum coverage required per base (default=10.0)")
    
    # Filtering parameters
    parser.add_argument("--min_length", type=int, default=100,
                       help="Minimum mapping length to consider (default=100)")
    parser.add_argument("--strand", choices=["both", "plus", "minus"], default="both",
                       help="Count hits on both strands, or only plus or minus strand hits (default: both)")
    parser.add_argument("--min_similarity", type=float, default=95.0,
                       help="Minimum percent identity to consider (default=95.0)")
    parser.add_argument("--uncovered_length_cutoff", type=int, default=0,
                       help="Save uncovered stretches longer than this cutoff (default=0)")
    
    # Output files
    parser.add_argument("--coverage_out", type=Path,
                       help="Output file for coverage data")
    parser.add_argument("--longest_uncovered_out", type=Path,
                       default=Path("longest_uncovered.txt"),
                       help="Output file for uncovered regions")
    parser.add_argument("--uncovered_fasta", type=Path,
                       help="Output FASTA file for uncovered regions")
    parser.add_argument("--n_split_fasta", type=Path,
                       help="Output FASTA file for N-split regions")
    
    # FASTA export parameters
    parser.add_argument("--extend_region", type=int, default=0,
                       help="Number of bases to extend uncovered regions on each side")
    parser.add_argument("--min_oligo_length", type=int, default=120,
                       help="Minimum sub-region length after N-splitting (default=120)")
    
    # System parameters
    parser.add_argument("--temp_dir", type=Path,
                       help="Parent directory for the run's temporary directory (default: system temp)")
    parser.add_argument("--log-level", "--log_level", dest="log_level",
                       choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                       default="INFO", help="Set logging level")


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