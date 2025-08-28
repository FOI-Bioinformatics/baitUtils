#!/usr/bin/env python3

"""
gap_filling.py

Main command for multi-pass gap filling to maximize oligo coverage.
Refactored from fill.py for better organization and maintainability.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

from Bio import SeqIO

from baitUtils._version import __version__
from baitUtils.gap_filling_algorithm import MultiPassSelector
from baitUtils.coverage_analysis import (
    CoverageAnalysisOrchestrator, 
    ForcedOligoHandler,
    SequenceLoader,
    CoverageCalculator
)


class GapFillingProcessor:
    """Processes gap filling commands with all options."""
    
    def __init__(self):
        """Initialize gap filling processor."""
        self.multi_pass_selector = MultiPassSelector()
        self.coverage_orchestrator = CoverageAnalysisOrchestrator()
        self.forced_handler = ForcedOligoHandler()
        self.sequence_loader = SequenceLoader()
        self.coverage_calculator = CoverageCalculator()
    
    def process_command(self, args) -> None:
        """
        Process the gap filling command with given arguments.
        
        Args:
            args: Parsed command-line arguments
        """
        # Set up logging
        self._setup_logging(args.log_level)
        
        # Load forced oligos
        forced_oligos = self.forced_handler.read_forced_oligos(args.forced_oligos)
        
        # Set up analysis
        bed, mappings_dict, genome_file = self.coverage_orchestrator.setup_analysis(
            args.psl, args.min_length, args.min_similarity, args.temp_dir
        )
        
        count = bed.count()
        if count == 0:
            logging.error("No intervals found in PSL after filtering. Exiting.")
            sys.exit(1)
        
        logging.info(f"Total oligos for selection: {count}")
        
        # Load reference sequences if provided
        sequences = None
        if args.reference_sequence:
            sequences = self.sequence_loader.load_reference_sequences(
                args.reference_sequence
            )
        
        # Create coverage calculator for multi-pass selection
        coverage_calc_func = self.coverage_orchestrator.create_coverage_calculator(
            genome_file, args.temp_dir
        )
        
        # Multi-pass selection
        selected_oligos = self.multi_pass_selector.multi_pass_selection(
            mappings_dict,
            forced_oligos,
            coverage_calc_func,
            args.min_coverage,
            args.max_coverage,
            args.spacing_distance,
            args.min_contribution,
            args.max_passes,
            args.max_oligos_per_pass,
            sequences,
            args.force,
            args.uncovered_length_cutoff,
            args.stall_rounds
        )
        
        # Write results
        self._write_results(args, selected_oligos, mappings_dict, genome_file)
    
    def _setup_logging(self, log_level: str = "INFO") -> None:
        """Configure logging settings."""
        logging.basicConfig(
            level=getattr(logging, log_level),
            format="%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%H:%M:%S"
        )
    
    def _write_results(self, args, selected_oligos, mappings_dict, genome_file) -> None:
        """Write out all results and perform final coverage analysis."""
        # Write selected oligos
        final_count = len(selected_oligos)
        logging.info(f"Final selection contains {final_count} oligos.")
        
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            for oligo_id in sorted(selected_oligos):
                f.write(f"{oligo_id}\n")
        logging.info(f"Selected oligos written to {args.output}")
        
        # Final coverage check
        coverage_bed = self.coverage_calculator.coverage_from_selected(
            selected_oligos, mappings_dict, genome_file, args.temp_dir
        )
        
        uncovered, bases_under, bases_over, total_bases, coverage_data, _ = \
            self.coverage_calculator.calculate_coverage(
                coverage_bed,
                args.min_coverage,
                args.max_coverage,
                args.temp_dir,
                genome_file
            )
        
        if total_bases > 0:
            logging.info(f"Final coverage check: Bases under min coverage: {bases_under:,} "
                        f"({bases_under/total_bases*100:.1f}%)")
        
        # Write optional outputs
        self._write_optional_outputs(args, coverage_data, uncovered)
    
    def _write_optional_outputs(self, args, coverage_data, uncovered) -> None:
        """Write optional output files."""
        # Coverage data
        if args.coverage_out:
            with open(args.coverage_out, "w") as f:
                f.write("Reference\tStart\tEnd\tCoverage\n")
                for ref, start, end, cov in coverage_data:
                    f.write(f"{ref}\t{start}\t{end}\t{cov}\n")
            logging.info(f"Coverage data written to {args.coverage_out}")
        
        # Uncovered regions
        if args.longest_uncovered_out:
            from baitUtils.coverage_checking import UncoveredRegionAnalyzer
            uncovered_list = [
                (ref, s, e, 0.0) 
                for ref, intervals in uncovered.items() 
                for (s, e) in intervals
            ]
            count_uncovered = UncoveredRegionAnalyzer.write_uncovered_regions(
                uncovered_list, args.longest_uncovered_out, args.uncovered_length_cutoff
            )
            logging.info(f"Wrote {count_uncovered} uncovered stretches >= "
                        f"{args.uncovered_length_cutoff}bp to {args.longest_uncovered_out}")
        
        # Uncovered FASTA
        if args.uncovered_fasta and args.fasta_reference:
            if not args.n_split_fasta:
                logging.warning("--n_split_fasta not specified, skipping N-split output")
            else:
                from baitUtils.coverage_checking import SequenceProcessor
                SequenceProcessor.export_uncovered_fasta(
                    uncovered,
                    args.fasta_reference,
                    args.uncovered_fasta,
                    args.n_split_fasta,
                    args.extend_region,
                    args.uncovered_length_cutoff,
                    args.min_oligo_length
                )


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add command-line arguments for gap filling."""
    parser.add_argument('--version', action='version', 
                       version=f"baitUtils gap-filling {__version__}")
    
    # Input files
    parser.add_argument("--psl", type=Path, required=True,
                       help="Path to PSL-like file")
    parser.add_argument("--forced_oligos", type=Path,
                       help="File with oligo IDs that must be included")
    parser.add_argument("--reference_sequence", type=Path,
                       help="Reference sequence file for sequence-based scoring")
    parser.add_argument("--fasta_reference", type=Path,
                       help="Reference FASTA file for exporting uncovered regions")
    
    # Coverage parameters
    parser.add_argument("--min_coverage", type=int, default=1,
                       help="Minimum coverage required per base (default=1)")
    parser.add_argument("--max_coverage", type=float,
                       help="Maximum coverage allowed per base")
    
    # Selection parameters
    parser.add_argument("--spacing_distance", type=int, default=30,
                       help="Min distance between start positions of selected oligos (default=30)")
    parser.add_argument("--min_contribution", type=int, default=5,
                       help="Minimum coverage contribution for selection (default=5)")
    parser.add_argument("--max_passes", type=int, default=5,
                       help="Maximum number of selection passes (default=5)")
    parser.add_argument("--max_oligos_per_pass", type=int,
                       help="Maximum number of oligos to select in each pass")
    parser.add_argument("--stall_rounds", type=int, default=3,
                       help="Number of rounds without improvement before stopping (default=3)")
    parser.add_argument("--force", action="store_true",
                       help="Continue selection even if uncovered regions don't decrease")
    
    # Filtering parameters
    parser.add_argument("--min_length", type=int, default=100,
                       help="Minimum mapping length to consider (default=100)")
    parser.add_argument("--min_similarity", type=float, default=95.0,
                       help="Minimum percent identity to consider (default=95.0)")
    parser.add_argument("--uncovered_length_cutoff", type=int, default=0,
                       help="Save uncovered stretches longer than this cutoff (default=0)")
    
    # Output files
    parser.add_argument("--output", type=Path, default=Path("selected_oligos.txt"),
                       help="Output file path for selected oligos")
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
                       help="Number of bases to extend uncovered regions by on each side")
    parser.add_argument("--min_oligo_length", type=int, default=120,
                       help="Minimum oligo length for N-split regions (default=120)")
    
    # System parameters
    parser.add_argument("--temp_dir", type=Path,
                       help="Directory for temporary files")
    parser.add_argument("--log_level", 
                       choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                       default="INFO", help="Set logging level")


def main(args) -> None:
    """Main function for gap filling command."""
    processor = GapFillingProcessor()
    try:
        processor.process_command(args)
    except Exception as e:
        logging.error(f"Gap filling failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Multi-pass gap filling to maximize oligo coverage.'
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)