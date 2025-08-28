#!/usr/bin/env python3

"""
sequence_statistics.py

Main command for calculating bait/oligo sequence statistics with filtering capabilities.
Refactored from stats.py for better organization and maintainability.
"""

import argparse
import logging
import multiprocessing
import os
import gzip
from pathlib import Path
from typing import List, Dict, Any, Optional, IO, TextIO
from Bio import SeqIO

from baitUtils._version import __version__
from baitUtils.sequence_analysis import SequenceAnalyzer, clean_sequence


class SequenceFilter:
    """Filter sequences based on various criteria."""
    
    def __init__(self, length: int = 120, complete: bool = False, no_ns: bool = False,
                 min_gc: Optional[float] = None, max_gc: Optional[float] = None,
                 min_tm: Optional[float] = None, max_tm: Optional[float] = None,
                 max_mask: Optional[float] = None):
        """
        Initialize sequence filter with criteria.
        
        Args:
            length: Target sequence length
            complete: Require sequences to be exact length
            no_ns: Exclude sequences with N bases
            min_gc: Minimum GC content percentage
            max_gc: Maximum GC content percentage
            min_tm: Minimum melting temperature
            max_tm: Maximum melting temperature
            max_mask: Maximum percentage of masked (lowercase) bases
        """
        self.length = length
        self.complete = complete
        self.no_ns = no_ns
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.max_mask = max_mask
        
        self.analyzer = SequenceAnalyzer()
    
    def passes_filters(self, sequence: str, stats: Dict[str, Any]) -> bool:
        """
        Check if sequence passes all filter criteria.
        
        Args:
            sequence: DNA sequence
            stats: Pre-calculated sequence statistics
            
        Returns:
            True if sequence passes all filters
        """
        # Length filter
        if self.complete and len(sequence) != self.length:
            return False
        
        # N bases filter
        if self.no_ns and stats['n_count'] > 0:
            return False
        
        # GC content filters
        if self.min_gc is not None and stats['gc_content'] < self.min_gc:
            return False
        if self.max_gc is not None and stats['gc_content'] > self.max_gc:
            return False
        
        # Melting temperature filters
        tm = stats['melting_temperature']
        if tm != 'NA':
            if self.min_tm is not None and tm < self.min_tm:
                return False
            if self.max_tm is not None and tm > self.max_tm:
                return False
        
        # Masked bases filter
        if self.max_mask is not None and stats['masked_percentage'] > self.max_mask:
            return False
        
        return True


class SequenceStatsCalculator:
    """Main class for calculating sequence statistics."""
    
    def __init__(self, num_processes: int = 1):
        """
        Initialize statistics calculator.
        
        Args:
            num_processes: Number of processes for parallel processing
        """
        self.num_processes = num_processes
        self.analyzer = SequenceAnalyzer()
    
    def process_sequences(self, input_file: str, seq_filter: Optional[SequenceFilter] = None,
                         sample_size: Optional[int] = None) -> List[Dict[str, Any]]:
        """
        Process sequences from FASTA file and calculate statistics.
        
        Args:
            input_file: Path to input FASTA file
            seq_filter: Optional sequence filter
            sample_size: Optional number of sequences to sample
            
        Returns:
            List of sequence statistics dictionaries
        """
        sequences = self._load_sequences(input_file, sample_size)
        
        if self.num_processes > 1:
            results = self._process_parallel(sequences)
        else:
            results = self._process_sequential(sequences)
        
        # Apply filters if provided
        if seq_filter:
            filtered_results = []
            for result in results:
                sequence = result.get('original_sequence', '')
                if seq_filter.passes_filters(sequence, result):
                    result['kept'] = True
                    filtered_results.append(result)
                else:
                    result['kept'] = False
                    filtered_results.append(result)
            return filtered_results
        else:
            # Mark all as kept if no filter
            for result in results:
                result['kept'] = True
            return results
    
    def _load_sequences(self, input_file: str, sample_size: Optional[int] = None) -> List[tuple]:
        """Load sequences from FASTA file."""
        sequences = []
        
        # Handle gzipped files
        if input_file.endswith('.gz'):
            file_handle = gzip.open(input_file, 'rt')
        else:
            file_handle = open(input_file, 'r')
        
        try:
            for i, record in enumerate(SeqIO.parse(file_handle, 'fasta')):
                if sample_size and i >= sample_size:
                    break
                sequences.append((record.id, str(record.seq)))
        finally:
            file_handle.close()
        
        return sequences
    
    def _process_sequential(self, sequences: List[tuple]) -> List[Dict[str, Any]]:
        """Process sequences sequentially."""
        results = []
        for seq_id, sequence in sequences:
            stats = self.analyzer.analyze_sequence(sequence, seq_id)
            stats['original_sequence'] = sequence
            results.append(stats)
        return results
    
    def _process_parallel(self, sequences: List[tuple]) -> List[Dict[str, Any]]:
        """Process sequences in parallel."""
        with multiprocessing.Pool(self.num_processes) as pool:
            results = pool.starmap(self._analyze_single_sequence, sequences)
        return results
    
    def _analyze_single_sequence(self, seq_id: str, sequence: str) -> Dict[str, Any]:
        """Analyze a single sequence (for parallel processing)."""
        stats = self.analyzer.analyze_sequence(sequence, seq_id)
        stats['original_sequence'] = sequence
        return stats
    
    def export_results(self, results: List[Dict[str, Any]], output_dir: Path,
                      save_filtered: bool = False) -> None:
        """
        Export results to files.
        
        Args:
            results: List of sequence statistics
            output_dir: Output directory
            save_filtered: Whether to save filtered FASTA
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Export statistics to TSV
        stats_file = output_dir / "sequence_statistics.tsv"
        self._export_tsv(results, stats_file)
        
        # Export filtered sequences if requested
        if save_filtered:
            kept_sequences = [r for r in results if r.get('kept', False)]
            if kept_sequences:
                filtered_file = output_dir / "filtered_sequences.fasta"
                self._export_fasta(kept_sequences, filtered_file)
                logging.info(f"Saved {len(kept_sequences)} filtered sequences to {filtered_file}")
        
        # Export summary
        summary_file = output_dir / "summary.txt"
        self._export_summary(results, summary_file)
        
        logging.info(f"Results exported to {output_dir}")
    
    def _export_tsv(self, results: List[Dict[str, Any]], output_file: Path) -> None:
        """Export results to TSV file."""
        if not results:
            return
        
        # Get all keys (excluding original_sequence and homopolymer_details)
        exclude_keys = {'original_sequence', 'homopolymer_details'}
        all_keys = set()
        for result in results:
            all_keys.update(k for k in result.keys() if k not in exclude_keys)
        
        header = sorted(all_keys)
        
        with open(output_file, 'w') as f:
            # Write header
            f.write('\t'.join(header) + '\n')
            
            # Write data
            for result in results:
                row = []
                for key in header:
                    value = result.get(key, 'NA')
                    if isinstance(value, float):
                        value = f"{value:.4f}"
                    row.append(str(value))
                f.write('\t'.join(row) + '\n')
    
    def _export_fasta(self, results: List[Dict[str, Any]], output_file: Path) -> None:
        """Export filtered sequences to FASTA."""
        with open(output_file, 'w') as f:
            for result in results:
                seq_id = result.get('sequence_id', 'unknown')
                sequence = result.get('original_sequence', '')
                f.write(f">{seq_id}\n{sequence}\n")
    
    def _export_summary(self, results: List[Dict[str, Any]], output_file: Path) -> None:
        """Export summary statistics."""
        if not results:
            return
        
        total_sequences = len(results)
        kept_sequences = len([r for r in results if r.get('kept', True)])
        filtered_out = total_sequences - kept_sequences
        
        with open(output_file, 'w') as f:
            f.write("SEQUENCE STATISTICS SUMMARY\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Total sequences analyzed: {total_sequences:,}\n")
            f.write(f"Sequences kept: {kept_sequences:,}\n")
            f.write(f"Sequences filtered out: {filtered_out:,}\n")
            if total_sequences > 0:
                f.write(f"Filter success rate: {kept_sequences/total_sequences*100:.1f}%\n\n")
            
            if kept_sequences > 0:
                f.write("STATISTICS FOR KEPT SEQUENCES:\n")
                f.write("-" * 30 + "\n")
                
                kept_results = [r for r in results if r.get('kept', True)]
                
                # Calculate summary stats
                gc_values = [r['gc_content'] for r in kept_results if isinstance(r['gc_content'], (int, float))]
                tm_values = [r['melting_temperature'] for r in kept_results if r['melting_temperature'] != 'NA']
                length_values = [r['length'] for r in kept_results]
                
                if gc_values:
                    f.write(f"GC content: {min(gc_values):.1f}% - {max(gc_values):.1f}% (mean: {sum(gc_values)/len(gc_values):.1f}%)\n")
                if tm_values:
                    f.write(f"Melting temperature: {min(tm_values):.1f}°C - {max(tm_values):.1f}°C (mean: {sum(tm_values)/len(tm_values):.1f}°C)\n")
                if length_values:
                    f.write(f"Sequence length: {min(length_values)} - {max(length_values)} bp (mean: {sum(length_values)/len(length_values):.0f} bp)\n")


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add command-line arguments for sequence statistics calculation."""
    parser.add_argument('--version', action='version', version=f"baitUtils sequence-stats {__version__}")
    
    # Input/output
    parser.add_argument('-i', '--input', required=True, 
                       help='Input FASTA or FASTA.GZ file')
    parser.add_argument('-o', '--outdir', required=True,
                       help='Output directory for results')
    
    # Sequence requirements
    parser.add_argument('-L', '--length', type=int, default=120,
                       help='Target bait length (default: 120)')
    parser.add_argument('-c', '--complete', action='store_true',
                       help='Require baits to be exact target length')
    parser.add_argument('-N', '--noNs', action='store_true',
                       help='Exclude sequences with N bases')
    
    # GC content filters
    parser.add_argument('-n', '--mingc', type=float,
                       help='Minimum GC content percentage')
    parser.add_argument('-x', '--maxgc', type=float,
                       help='Maximum GC content percentage')
    
    # Melting temperature filters
    parser.add_argument('-q', '--mint', type=float,
                       help='Minimum melting temperature (°C)')
    parser.add_argument('-z', '--maxt', type=float,
                       help='Maximum melting temperature (°C)')
    
    # Other filters
    parser.add_argument('-K', '--maxmask', type=float,
                       help='Maximum percentage of masked (lowercase) bases')
    
    # Processing options
    parser.add_argument('-p', '--processes', type=int,
                       default=max(1, multiprocessing.cpu_count() - 1),
                       help='Number of processes for parallel computation')
    parser.add_argument('-s', '--sample', type=int,
                       help='Sample only N sequences for analysis')
    
    # Output options
    parser.add_argument('-f', '--filter', action='store_true',
                       help='Save filtered sequences to FASTA')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose logging')


def main(args) -> None:
    """Main function for sequence statistics calculation."""
    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Validate input file
    if not os.path.exists(args.input):
        logging.error(f"Input file not found: {args.input}")
        return
    
    # Initialize filter
    seq_filter = SequenceFilter(
        length=args.length,
        complete=args.complete,
        no_ns=args.noNs,
        min_gc=args.mingc,
        max_gc=args.maxgc,
        min_tm=args.mint,
        max_tm=args.maxt,
        max_mask=args.maxmask
    )
    
    # Initialize calculator
    calculator = SequenceStatsCalculator(num_processes=args.processes)
    
    # Process sequences
    logging.info(f"Processing sequences from {args.input}")
    if args.sample:
        logging.info(f"Sampling {args.sample} sequences")
    
    results = calculator.process_sequences(
        args.input, 
        seq_filter=seq_filter,
        sample_size=args.sample
    )
    
    # Export results
    calculator.export_results(results, Path(args.outdir), save_filtered=args.filter)
    
    # Print summary
    total = len(results)
    kept = len([r for r in results if r.get('kept', False)])
    logging.info(f"Analysis complete: {kept}/{total} sequences passed filters")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate sequence statistics with filtering')
    add_arguments(parser)
    args = parser.parse_args()
    main(args)