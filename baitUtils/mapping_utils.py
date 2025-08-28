#!/usr/bin/env python3

"""
mapping_utils.py

Core mapping utilities for sequence mapping operations using external tools.
Provides functionality for running mappers, parsing results, and managing outputs.
"""

import logging
import os
import sys
import subprocess
from typing import Set, Dict, Optional, List
from pathlib import Path
from Bio import SeqIO


class SequenceLoader:
    """Handles loading and processing of FASTA sequences."""
    
    @staticmethod
    def load_sequence_records(fasta_path: str) -> Dict[str, SeqIO.SeqRecord]:
        """
        Load sequence records from FASTA file.
        
        Args:
            fasta_path: Path to FASTA file
            
        Returns:
            Dictionary mapping sequence IDs to SeqRecord objects
        """
        try:
            seq_records = SeqIO.to_dict(SeqIO.parse(fasta_path, 'fasta'))
            logging.info(f"Loaded {len(seq_records)} sequences from {fasta_path}")
            return seq_records
        except Exception as e:
            logging.error(f"Error reading FASTA file {fasta_path}: {e}")
            raise

    @staticmethod
    def count_sequences(fasta_path: str) -> int:
        """
        Count number of sequences in FASTA file without loading all into memory.
        
        Args:
            fasta_path: Path to FASTA file
            
        Returns:
            Number of sequences in file
        """
        try:
            count = sum(1 for _ in SeqIO.parse(fasta_path, 'fasta'))
            logging.info(f"Counted {count} sequences in {fasta_path}")
            return count
        except Exception as e:
            logging.error(f"Error counting sequences in {fasta_path}: {e}")
            raise

    @staticmethod
    def get_sequence_ids(fasta_path: str) -> Set[str]:
        """
        Get set of all sequence IDs from FASTA file.
        
        Args:
            fasta_path: Path to FASTA file
            
        Returns:
            Set of sequence IDs
        """
        try:
            ids = set()
            for record in SeqIO.parse(fasta_path, 'fasta'):
                ids.add(record.id)
            logging.info(f"Extracted {len(ids)} sequence IDs from {fasta_path}")
            return ids
        except Exception as e:
            logging.error(f"Error extracting sequence IDs from {fasta_path}: {e}")
            raise


class PblatRunner:
    """Handles running pblat mapping tool."""
    
    @staticmethod
    def run_pblat(
        baits_file: str,
        target_file: str,
        output_file: str,
        threads: int = 1,
        min_match: int = 2,
        min_score: int = 30,
        min_identity: int = 90
    ) -> None:
        """
        Run pblat to map sequences against target.
        
        Args:
            baits_file: Path to input FASTA file with sequences to map
            target_file: Path to target FASTA file to map against
            output_file: Path to output PSL file
            threads: Number of threads to use
            min_match: Minimum number of tile matches
            min_score: Minimum score
            min_identity: Minimum sequence identity percentage
        """
        cmd = [
            'pblat',
            f'-threads={threads}',
            f'-minMatch={min_match}',
            f'-minScore={min_score}',
            f'-minIdentity={min_identity}',
            target_file,
            baits_file,
            output_file
        ]

        logging.info(f"Running pblat: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, 
                          stderr=subprocess.DEVNULL)
            logging.info(f"pblat completed successfully. Output: {output_file}")
        except subprocess.CalledProcessError as e:
            logging.error(f"pblat failed with return code {e.returncode}")
            raise
        except FileNotFoundError:
            logging.error("pblat not found in PATH. Please install pblat.")
            raise
        except Exception as e:
            logging.error(f"Error running pblat: {e}")
            raise


class PSLParser:
    """Parses PSL files and extracts mapping information."""
    
    @staticmethod
    def parse_psl_file(
        psl_file: str,
        min_identity: float = 90.0,
        min_match_count: int = 0,
        filtered_output: Optional[str] = None
    ) -> Set[str]:
        """
        Parse PSL file and return set of sequence IDs that meet criteria.
        
        Args:
            psl_file: Path to PSL file
            min_identity: Minimum identity percentage for filtering
            min_match_count: Minimum number of matching bases required
            filtered_output: Optional path to save filtered PSL entries
            
        Returns:
            Set of sequence IDs that passed filters
        """
        mapped_sequences = set()
        header_lines = []
        filtered_hits = []
        
        try:
            with open(psl_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    
                    # Handle header lines
                    if PSLParser._is_header_line(line):
                        if filtered_output:
                            header_lines.append(line)
                        continue
                    
                    if not line:
                        continue
                    
                    # Parse PSL line
                    hit_info = PSLParser._parse_psl_line(line)
                    if hit_info is None:
                        continue
                    
                    matches, identity_pct, seq_id = hit_info
                    
                    # Apply filters
                    if matches >= min_match_count and identity_pct >= min_identity:
                        mapped_sequences.add(seq_id)
                        if filtered_output:
                            filtered_hits.append(line)

            # Write filtered output if requested
            if filtered_output and filtered_hits:
                PSLParser._write_filtered_psl(filtered_output, header_lines, filtered_hits)
                logging.info(f"Wrote {len(filtered_hits)} filtered hits to {filtered_output}")

        except Exception as e:
            logging.error(f"Error parsing PSL file {psl_file}: {e}")
            raise
            
        return mapped_sequences

    @staticmethod
    def _is_header_line(line: str) -> bool:
        """Check if line is a PSL header line."""
        return (line.startswith('psLayout') or 
                line.startswith('-') or 
                line.startswith('match') or 
                line.startswith('no matches'))

    @staticmethod
    def _parse_psl_line(line: str) -> Optional[tuple]:
        """
        Parse a single PSL data line.
        
        Returns:
            Tuple of (matches, identity_percentage, sequence_id) or None if invalid
        """
        cols = line.split()
        if len(cols) < 21:
            return None
        
        try:
            matches = int(cols[0])
            mismatches = int(cols[1])
            rep_matches = int(cols[2])
            n_count = int(cols[3])
            q_num_insert = int(cols[4])
            seq_id = cols[9]
            
            # Calculate identity percentage
            size = matches + mismatches + rep_matches + q_num_insert
            if size == 0:
                return None
            
            identity_pct = 100.0 * (matches + rep_matches) / size
            
            return matches, identity_pct, seq_id
            
        except (ValueError, IndexError) as e:
            logging.warning(f"Invalid PSL line format: {line[:50]}...")
            return None

    @staticmethod
    def _write_filtered_psl(output_path: str, headers: List[str], hits: List[str]) -> None:
        """Write filtered PSL hits to output file."""
        try:
            with open(output_path, 'w') as f:
                # Write headers
                for header in headers:
                    f.write(header + '\n')
                # Write hits
                for hit in hits:
                    f.write(hit + '\n')
        except Exception as e:
            logging.error(f"Error writing filtered PSL file {output_path}: {e}")
            raise


class MappingResultsWriter:
    """Writes mapping results to various output formats."""
    
    @staticmethod
    def write_sequence_ids(sequence_ids: Set[str], output_file: str) -> None:
        """
        Write sequence IDs to a text file, one per line.
        
        Args:
            sequence_ids: Set of sequence IDs to write
            output_file: Path to output text file
        """
        try:
            with open(output_file, 'w') as f:
                for seq_id in sorted(sequence_ids):
                    f.write(seq_id + '\n')
            logging.info(f"Wrote {len(sequence_ids)} sequence IDs to {output_file}")
        except Exception as e:
            logging.error(f"Error writing sequence IDs to {output_file}: {e}")
            raise

    @staticmethod
    def write_sequences_fasta(
        sequence_ids: Set[str],
        seq_records: Dict[str, SeqIO.SeqRecord],
        output_file: str
    ) -> None:
        """
        Write selected sequences to FASTA file.
        
        Args:
            sequence_ids: Set of sequence IDs to write
            seq_records: Dictionary of sequence records
            output_file: Path to output FASTA file
        """
        try:
            written_count = 0
            with open(output_file, 'w') as f:
                for seq_id in sequence_ids:
                    if seq_id in seq_records:
                        SeqIO.write(seq_records[seq_id], f, 'fasta')
                        written_count += 1
                    else:
                        logging.warning(f"Sequence ID {seq_id} not found in records")
            
            logging.info(f"Wrote {written_count} sequences to {output_file}")
        except Exception as e:
            logging.error(f"Error writing FASTA file {output_file}: {e}")
            raise


class SequenceMapper:
    """High-level interface for sequence mapping operations."""
    
    def __init__(self):
        """Initialize sequence mapper."""
        self.sequence_loader = SequenceLoader()
        self.pblat_runner = PblatRunner()
        self.psl_parser = PSLParser()
        self.results_writer = MappingResultsWriter()
    
    def map_sequences(
        self,
        input_file: str,
        target_file: str,
        output_dir: str,
        output_prefix: str,
        mapper: str = 'pblat',
        threads: int = 1,
        min_match: int = 2,
        min_score: int = 30,
        min_identity: int = 90
    ) -> str:
        """
        Map sequences against target using specified mapper.
        
        Args:
            input_file: Path to input FASTA file
            target_file: Path to target FASTA file
            output_dir: Output directory
            output_prefix: Output file prefix
            mapper: Mapping tool to use
            threads: Number of threads
            min_match: Minimum tile matches
            min_score: Minimum score
            min_identity: Minimum identity percentage
            
        Returns:
            Path to mapping output file
        """
        # Determine output file path
        mapping_output = os.path.join(output_dir, f"{output_prefix}-mapping.psl")
        
        if mapper == 'pblat':
            self.pblat_runner.run_pblat(
                input_file, target_file, mapping_output,
                threads, min_match, min_score, min_identity
            )
        else:
            raise ValueError(f"Unsupported mapper: {mapper}")
        
        return mapping_output


def create_output_directory(directory: str) -> None:
    """
    Create output directory if it doesn't exist.
    
    Args:
        directory: Path to directory to create
    """
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
            logging.info(f"Created output directory: {directory}")
        except Exception as e:
            logging.error(f"Failed to create output directory '{directory}': {e}")
            raise
    else:
        logging.info(f"Output directory already exists: {directory}")


def validate_identity_parameters(min_identity: int, filter_identity: int) -> None:
    """
    Validate that identity parameters are consistent.
    
    Args:
        min_identity: Minimum identity for mapping
        filter_identity: Minimum identity for filtering
    """
    if filter_identity < min_identity:
        raise ValueError(
            f"filter_identity ({filter_identity}) must be >= min_identity ({min_identity})"
        )