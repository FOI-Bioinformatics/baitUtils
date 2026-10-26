#!/usr/bin/env python3

"""
sequence_mapping.py

Main command for mapping sequences against reference genomes using external tools.
Refactored from map.py for better organization and maintainability.
"""

import argparse
import logging
import os
import sys

from baitUtils._version import __version__
from baitUtils.cli_utils import (add_argument, add_logging_arguments, logging_kwargs, percent,
                                 positive_int, non_negative_int)
from baitUtils.logging_utils import ensure_logging
from baitUtils.mapping_utils import (
    SequenceMapper,
    SequenceLoader,
    PSLParser,
    MappingResultsWriter,
    alignment_format,
    create_output_directory,
    validate_identity_parameters
)


class SequenceMappingProcessor:
    """Processes sequence mapping commands with all options."""
    
    def __init__(self):
        """Initialize sequence mapping processor."""
        self.sequence_mapper = SequenceMapper()
        self.sequence_loader = SequenceLoader()
        self.psl_parser = PSLParser()
        self.results_writer = MappingResultsWriter()
    
    def process_command(self, args) -> None:
        """
        Process the sequence mapping command with given arguments.
        
        Args:
            args: Parsed command-line arguments
        """
        # Set up logging
        ensure_logging(**logging_kwargs(args))
        
        # Mapper-specific options and their defaults
        resolve_mapper_options(args)
        if args.mapper == 'pblat':
            validate_identity_parameters(args.min_identity, args.filter_identity)
        
        # Create output directory
        create_output_directory(args.outdir)
        
        # Load input sequences and count target sequences
        seq_records, all_sequence_ids = self._load_input_data(args)
        
        # Map sequences
        mapping_output = self._perform_mapping(args)
        
        # Parse mapping results
        mapped_sequences = self._parse_mapping_results(args, mapping_output)
        
        # Compute unmapped sequences
        unmapped_sequences = all_sequence_ids - mapped_sequences
        
        # Log results
        self._log_mapping_results(mapped_sequences, unmapped_sequences)
        
        # Write results
        self._write_results(args, mapped_sequences, unmapped_sequences, seq_records)
        
        logging.info("Mapping complete.")
    
    def _load_input_data(self, args) -> tuple:
        """Load input sequences and validate target file."""
        logging.info("Reading input FASTA file...")
        seq_records = self.sequence_loader.load_sequence_records(args.input)
        all_sequence_ids = set(seq_records.keys())
        
        return seq_records, all_sequence_ids
    
    def _perform_mapping(self, args) -> str:
        """Perform sequence mapping."""
        logging.info("Mapping sequences against target genome...")
        return self.sequence_mapper.map_sequences(
            args.input,
            args.reference,
            args.outdir,
            args.outprefix,
            args.mapper,
            args.threads,
            args.min_match,
            args.min_score,
            args.min_identity,
            args.minimap2_preset
        )
    
    def _parse_mapping_results(self, args, mapping_output: str) -> set:
        """Parse mapping output and apply filters."""
        logging.info("Parsing mapping output...")
        
        # Filtered output keeps the format of the mapper output
        mapping_filtered_output = os.path.join(
            args.outdir, f"{args.outprefix}-mapping_filtered.{alignment_format(mapping_output)}"
        )
        
        # Parse with filters
        mapped_sequences = self.psl_parser.parse_psl_file(
            mapping_output,
            args.filter_identity,
            args.min_match_count,
            mapping_filtered_output,
            args.min_length
        )
        
        # Per-bait hit table for off-target assessment
        self.hit_table = self.psl_parser.last_hit_table
        hits_file = os.path.join(args.outdir, f"{args.outprefix}-hits.tsv")
        self.hit_table.to_csv(hits_file, sep="\t", index=False)
        logging.info(f"Wrote per-bait hit table to {hits_file}")
        
        if args.max_hits is not None:
            multi = set(self.hit_table.loc[self.hit_table['n_hits'] > args.max_hits, 'oligo_id'])
            if multi:
                logging.info(f"Excluding {len(multi)} baits with more than {args.max_hits} hits")
            mapped_sequences = mapped_sequences - multi
        
        return mapped_sequences
    
    def _log_mapping_results(self, mapped_sequences: set, unmapped_sequences: set) -> None:
        """Log mapping results statistics."""
        logging.info(f"Number of mapped sequences: {len(mapped_sequences)}")
        logging.info(f"Number of unmapped sequences: {len(unmapped_sequences)}")
        
        total = len(mapped_sequences) + len(unmapped_sequences)
        if total > 0:
            mapped_pct = len(mapped_sequences) / total * 100
            logging.info(f"Mapping success rate: {mapped_pct:.1f}%")
    
    def _write_results(self, args, mapped_sequences: set, unmapped_sequences: set, 
                      seq_records: dict) -> None:
        """Write all result files."""
        # Write ID files
        mapped_ids_file = os.path.join(args.outdir, f"{args.outprefix}-mapped-sequence-ids.txt")
        unmapped_ids_file = os.path.join(args.outdir, f"{args.outprefix}-unmapped-sequence-ids.txt")
        
        self.results_writer.write_sequence_ids(mapped_sequences, mapped_ids_file)
        self.results_writer.write_sequence_ids(unmapped_sequences, unmapped_ids_file)
        
        # Write FASTA files based on option
        if args.fasta_output != 'none':
            self._write_fasta_outputs(args, mapped_sequences, unmapped_sequences, seq_records)
    
    def _write_fasta_outputs(self, args, mapped_sequences: set, unmapped_sequences: set,
                           seq_records: dict) -> None:
        """Write FASTA output files based on user selection."""
        if args.fasta_output in ['mapped', 'both']:
            fasta_file = os.path.join(args.outdir, f"{args.outprefix}-mapped-sequences.fa")
            logging.info(f"Writing mapped sequences to FASTA file: {fasta_file}")
            reverse = set()
            if args.orient_to_reference and self.hit_table is not None and len(self.hit_table):
                minus = self.hit_table.loc[self.hit_table['best_strand'] == '-', 'oligo_id']
                reverse = set(minus) & set(mapped_sequences)
                logging.info(f"Reverse-complementing {len(reverse)} baits whose best hit is on the minus strand")
            self.results_writer.write_sequences_fasta(mapped_sequences, seq_records, fasta_file, reverse)
        
        if args.fasta_output in ['unmapped', 'both']:
            fasta_file = os.path.join(args.outdir, f"{args.outprefix}-unmapped-sequences.fa")
            logging.info(f"Writing unmapped sequences to FASTA file: {fasta_file}")
            self.results_writer.write_sequences_fasta(unmapped_sequences, seq_records, fasta_file)



PBLAT_DEFAULTS = {'min_match': 2, 'min_score': 30, 'min_identity': 90.0}
MINIMAP2_DEFAULT_PRESET = 'sr'


def resolve_mapper_options(args) -> None:
    """
    Fill in the mapper-specific defaults and warn when an option of the
    other mapper was given, since it would be silently ignored.
    """
    if args.mapper == 'pblat':
        if args.minimap2_preset is not None:
            logging.warning("--minimap2-preset has no effect with --mapper pblat")
        args.minimap2_preset = MINIMAP2_DEFAULT_PRESET
    else:
        given = [f"--{key.replace('_', '-')}" for key in PBLAT_DEFAULTS if getattr(args, key) is not None]
        if given:
            logging.warning(f"{', '.join(given)} apply to pblat only and have no effect with --mapper minimap2")
        if args.minimap2_preset is None:
            args.minimap2_preset = MINIMAP2_DEFAULT_PRESET
    for key, default in PBLAT_DEFAULTS.items():
        if getattr(args, key) is None:
            setattr(args, key, default)


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add command-line arguments for sequence mapping."""
    parser.add_argument('--version', action='version', version=f"baitUtils map {__version__}")
    
    # Input/Output
    add_argument(parser, '-i', '--input', required=True,
                 help='Input sequences FASTA file')
    add_argument(parser, '-r', '--reference', dest='reference', required=True, aliases=('-q', '--query'),
                 help='Reference FASTA file to map against')
    add_argument(parser, '-o', '--output', dest='outdir', default='.', aliases=('--outdir',),
                 help='Output directory (default: current directory)')
    parser.add_argument('--prefix', dest='outprefix', default='out',
                        help='Prefix for output file names (default: out)')
    
    # Mapping parameters
    parser.add_argument('--mapper', choices=['pblat', 'minimap2'], default='pblat',
                        help='Mapping tool (default: pblat). minimap2 is run with -c and --minimap2-preset')
    parser.add_argument('--minimap2-preset', dest='minimap2_preset', default=None,
                        help=f'minimap2 -x preset when --mapper minimap2 (default: {MINIMAP2_DEFAULT_PRESET})')
    add_argument(parser, '--threads', type=positive_int, default=1, aliases=('-X',),
                 help='Number of threads (default: 1)')
    
    # pblat parameters
    add_argument(parser, '--min-match', dest='min_match', type=int, default=None, aliases=('--minMatch',),
                 help='pblat -minMatch, number of tile matches (default: 2)')
    add_argument(parser, '--min-score', dest='min_score', type=int, default=None, aliases=('--minScore',),
                 help='pblat -minScore (default: 30)')
    add_argument(parser, '--min-identity', dest='min_identity', type=percent, default=None,
                 aliases=('--minIdentity',),
                 help='pblat -minIdentity, percent 0-100 (default: 90)')
    
    # Filtering parameters
    add_argument(parser, '--filter-identity', dest='filter_identity', type=percent, default=90.0,
                 aliases=('--filterIdentity',),
                 help='Keep hits with BLAT identity at or above this percent, 0-100 '
                      '(with pblat, must be >= --min-identity; default: 90)')
    add_argument(parser, '--min-match-count', dest='min_match_count', type=non_negative_int, default=0,
                 aliases=('--minMatchCount',),
                 help='Minimum number of matching bases required (default: 0)')
    parser.add_argument('--min-length', dest='min_length', type=non_negative_int, default=0,
                        help='Minimum aligned length in bases for a hit to count (default: 0)')
    parser.add_argument('--max-hits', dest='max_hits', type=non_negative_int, default=None,
                        help='Treat baits with more passing hits than this as unmapped '
                             '(off-target filter; default: no limit)')
    
    # Output options
    parser.add_argument('--orient-to-reference', dest='orient_to_reference', action='store_true',
                        help='Write mapped baits in reference orientation: reverse-complement '
                             'baits whose best hit is on the minus strand')
    parser.add_argument('--fasta-output',
                        choices=['mapped', 'unmapped', 'both', 'none'],
                        default='mapped',
                        help='Which sequences to include in FASTA output '
                             '(default: mapped)')
    
    add_logging_arguments(parser, verbose_aliases=('-l', '--log'))

def main(args) -> None:
    """Main function for sequence mapping command."""
    processor = SequenceMappingProcessor()
    try:
        processor.process_command(args)
    except Exception as e:
        logging.error(f"Sequence mapping failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Map sequences against reference genomes using external tools.'
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)