#!/usr/bin/env python3

"""
gap_filling.py

Main command for multi-pass gap filling to maximize oligo coverage.
Refactored from fill.py for better organization and maintainability.
"""

import argparse
import logging
import sys
import tempfile
from pathlib import Path


from baitUtils._version import __version__
from baitUtils.cli_utils import add_argument, add_logging_arguments, logging_kwargs, positive_int, non_negative_int
from baitUtils.coverage_evaluation import add_coverage_arguments
from baitUtils.logging_utils import ensure_logging
from baitUtils.bedtools_support import require_bedtools
from baitUtils.gap_filling_algorithm import MultiPassSelector
from baitUtils.coverage_analysis import (
    CoverageAnalysisOrchestrator,
    ForcedOligoHandler,
    SequenceLoader,
    CoverageCalculator,
    SequenceProcessor,
    write_uncovered_regions
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
        ensure_logging(**logging_kwargs(args))
        require_bedtools()

        import pybedtools
        forced_oligos = self.forced_handler.read_forced_oligos(args.forced_oligos)
        
        if args.temp_dir:
            Path(args.temp_dir).mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="baitUtils_fill_", dir=args.temp_dir) as tmp:
            work_dir = Path(tmp)
            previous_tempdir = tempfile.tempdir
            pybedtools.set_tempdir(str(work_dir))  # sets tempfile.tempdir globally
            try:
                self._run(args, forced_oligos, work_dir)
            finally:
                pybedtools.cleanup(remove_all=True)
                tempfile.tempdir = previous_tempdir
    
    def _run(self, args, forced_oligos, work_dir: Path) -> None:
        """Run selection with all temporary files under work_dir."""
        args.temp_dir = work_dir
        # Set up analysis
        bed, mappings_dict, genome_file = self.coverage_orchestrator.setup_analysis(
            args.alignments, args.min_length, args.min_identity, args.temp_dir, args.reference,
            strand=args.strand
        )
        
        count = bed.count()
        if count == 0:
            logging.error("No intervals found in PSL after filtering. Exiting.")
            sys.exit(1)
        
        logging.info(f"Total oligos for selection: {count}")
        if args.min_coverage > 1:
            logging.warning(
                f"--min-coverage {args.min_coverage:g}: the spacing constraint keeps selected oligos "
                f"at least {args.spacing_distance} bp apart, so depth above 1 at a position can only "
                "come from different oligos that start at least that far apart")
        
        # Reference sequences for the GC term of the oligo score
        sequences = self.sequence_loader.load_reference_sequences(args.reference)
        
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
            args.spacing_distance,
            args.min_contribution,
            args.max_passes,
            args.max_oligos_per_pass,
            sequences,
            args.uncovered_length_cutoff,
            args.stall_rounds
        )
        
        # Write results
        self._write_results(args, selected_oligos, mappings_dict, genome_file)
    
    def _write_results(self, args, selected_oligos, mappings_dict, genome_file) -> None:
        """Write out all results and perform final coverage analysis."""
        # Write selected oligos (IDs) and the selected mappings (loci)
        selected_ids = sorted({key[0] for key in selected_oligos})
        logging.info(f"Final selection contains {len(selected_ids)} oligos "
                     f"({len(selected_oligos)} mappings).")
        
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            for oligo_id in selected_ids:
                f.write(f"{oligo_id}\n")
        logging.info(f"Selected oligos written to {args.output}")
        
        mappings_out = args.output.with_name(args.output.stem + "_mappings.tsv")
        with open(mappings_out, "w") as f:
            f.write("oligo_id\treference\tstart\tend\n")
            for ref_id, mlist in sorted(mappings_dict.items()):
                for m in sorted(mlist, key=lambda x: x.start):
                    if m.key in selected_oligos:
                        f.write(f"{m.oligo_id}\t{ref_id}\t{m.start}\t{m.end}\n")
        logging.info(f"Selected mappings written to {mappings_out}")
        
        # Final coverage check
        coverage_bed = self.coverage_calculator.coverage_from_selected(
            selected_oligos, mappings_dict, genome_file, args.temp_dir
        )
        
        uncovered, bases_under, bases_over, total_bases, coverage_data, _ = \
            self.coverage_calculator.calculate_coverage(
                coverage_bed,
                args.min_coverage,
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
            count_uncovered = write_uncovered_regions(
                uncovered, args.longest_uncovered_out, args.uncovered_length_cutoff
            )
            logging.info(f"Wrote {count_uncovered} uncovered stretches >= "
                        f"{args.uncovered_length_cutoff}bp to {args.longest_uncovered_out}")
        
        # Uncovered FASTA
        if args.uncovered_fasta and args.reference:
            if not args.n_split_fasta:
                logging.info("--n_split_fasta not given; writing the uncovered FASTA only")
            SequenceProcessor.export_uncovered_fasta(
                uncovered,
                args.reference,
                args.uncovered_fasta,
                args.n_split_fasta,
                args.extend_region,
                args.uncovered_length_cutoff,
                args.min_oligo_length
            )



def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add command-line arguments for gap filling."""
    parser.add_argument('--version', action='version', version=f"baitUtils fill {__version__}")
    add_coverage_arguments(parser, "Minimum depth the selection aims for at every base (default: 1)")
    
    # Selection parameters
    add_argument(parser, "--spacing-distance", dest="spacing_distance", type=non_negative_int, default=30,
                 aliases=("--spacing_distance",),
                 help="Minimum distance between start positions of selected oligos (default: 30)")
    add_argument(parser, "--min-contribution", dest="min_contribution", type=non_negative_int, default=5,
                 aliases=("--min_contribution",),
                 help="Minimum number of uncovered bases an oligo must add to be selected (default: 5)")
    add_argument(parser, "--max-passes", dest="max_passes", type=positive_int, default=5,
                 aliases=("--max_passes",),
                 help="Maximum number of selection passes (default: 5)")
    add_argument(parser, "--max-oligos-per-pass", dest="max_oligos_per_pass", type=positive_int,
                 aliases=("--max_oligos_per_pass",),
                 help="Maximum number of oligos to select in each pass (default: no limit)")
    add_argument(parser, "--stall-rounds", dest="stall_rounds", type=positive_int, default=3,
                 aliases=("--stall_rounds",),
                 help="Passes without fewer uncovered regions before stopping (default: 3)")
    
    # Output
    parser.add_argument("--output", type=Path, default=Path("selected_oligos.txt"),
                        help="Output file for the selected oligo IDs (default: selected_oligos.txt); "
                             "the selected loci go to <output>_mappings.tsv")
    add_logging_arguments(parser)

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