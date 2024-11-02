#!/usr/bin/env python3

"""
stats.py

Calculate bait statistics from a FASTA file.
"""

import argparse
from baitUtils._version import __version__  # Import version from _version.py
import logging
import multiprocessing
import os
import sys
import random
import math
import gzip
from typing import Tuple, Dict, Any, List, Optional, Union, IO
from Bio import SeqIO, Seq
from Bio.SeqUtils import gc_fraction, MeltingTemp as mt
import RNA  # ViennaRNA package for MFE calculation
from Bio.Align import PairwiseAligner
from collections import Counter
from math import log2

def add_arguments(parser: argparse.ArgumentParser):

    """
    Add command-line arguments specific to the 'stats' subcommand.

    Args:
        parser (argparse.ArgumentParser): The argument parser to which arguments will be added.
    """
    parser.add_argument('--version', action='version', version=f"baitUtils stats {__version__}")
    parser.add_argument('-i', '--input', required=True, help='Input FASTA or FASTA.GZ file')
    parser.add_argument('-L', '--length', type=int, default=120, help='Requested bait length (Default = 120)')
    parser.add_argument('-c', '--complete', action='store_true', help='Require baits be full length')
    parser.add_argument('-N', '--noNs', action='store_true', help='Exclude bait sequences with Ns')
    parser.add_argument('-n', '--mingc', type=float, help='Minimum GC content in %%')
    parser.add_argument('-x', '--maxgc', type=float, help='Maximum GC content in %%')
    parser.add_argument('-q', '--mint', type=float, help='Minimum melting temperature in °C')
    parser.add_argument('-z', '--maxt', type=float, help='Maximum melting temperature in °C')
    parser.add_argument('-K', '--maxmask', type=float, help='Maximum %% sequence consisting of masked elements')
    parser.add_argument('-J', '--maxhomopoly', type=int, help='Maximum homopolymer length')
    parser.add_argument('-y', '--minlc', type=float, help='Minimum sequence linguistic complexity')
    parser.add_argument('--minMFE', type=float, help='Minimum MFE (kcal/mol) to keep a bait')
    parser.add_argument('--maxMFE', type=float, help='Maximum MFE (kcal/mol) to keep a bait')
    parser.add_argument('--maxSelfScore', type=float, help='Maximum self-alignment score to keep a bait')
    parser.add_argument('--minEntropy', type=float, help='Minimum Shannon entropy to keep a bait')
    parser.add_argument('-o', '--outprefix', default='out', help='Output file prefix (Default = out)')
    parser.add_argument('-Z', '--outdir', default='.', help='Output directory path (Default is ./)')
    parser.add_argument('-l', '--log', action='store_true', help='Output detailed log')
    parser.add_argument('-C', '--collapse', action='store_true', help='Collapse ambiguities to a single nucleotide')
    parser.add_argument('-Y', '--rna', action='store_true', help='Output baits as RNA rather than DNA')
    parser.add_argument('-R', '--rc', action='store_true', help='Output reverse-complemented baits')
    parser.add_argument('-G', '--gaps', choices=['include', 'exclude', 'extend'], default='include',
                        help='Strategy to handle sequence gaps (-) (include, exclude, or extend) (Default = include)')
    parser.add_argument('-X', '--threads', type=int, default=1, help='Number of threads (Default = 1)')
    parser.add_argument('--rng', type=int, help='Random number seed (Default uses system entropy)')
    parser.add_argument('--gzip', action='store_true', help='Gzip output files')
    parser.add_argument('--filter', action='store_true', help='Save the filtered FASTA file')

def main(args):
    """
    Main function for the 'stats' subcommand.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    setup_logging(args)
    # Create output directory if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)
    # Set random seed if provided
    if args.rng:
        random.seed(args.rng)
    # Read sequences
    logging.info("Reading input sequences...")
    try:
        with open_fasta(args.input) as handle:
            seq_list: List[SeqIO.SeqRecord] = list(SeqIO.parse(handle, 'fasta'))
            if not seq_list:
                raise ValueError("No sequences found in the input FASTA file.")
    except Exception as e:
        logging.error(f"Error reading input FASTA file: {e}")
        sys.exit(1)
    logging.info(f"Total baits read: {len(seq_list)}")
    # Process sequences
    process_all_sequences(seq_list, args)
    logging.info("Processing complete.")
    # Additional verification
    # Parameters file is always generated
    final_params_file: str = os.path.join(args.outdir, f"{args.outprefix}-filtered-params.txt")
    try:
        with open(final_params_file, 'r') as pf:
            total_params_lines: int = sum(1 for line in pf) - 1  # Subtract header
        logging.info(f"Total entries in parameters file: {total_params_lines}")
        if total_params_lines != len(seq_list):
            logging.warning(f"Total entries in parameters file ({total_params_lines}) do not match total input sequences ({len(seq_list)}).")
    except Exception as e:
        logging.error(f"Error reading parameters file: {e}")

def setup_logging(args: argparse.Namespace) -> None:
    """
    Configure the logging settings based on the provided arguments.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    log_level = logging.DEBUG if args.log else logging.INFO
    logging.basicConfig(level=log_level, format='%(asctime)s %(levelname)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

def max_homopolymer_length(seq: str) -> int:
    """
    Calculate the maximum homopolymer length in a nucleotide sequence.

    Args:
        seq (str): DNA or RNA sequence.

    Returns:
        int: Length of the longest consecutive identical nucleotide.
    """
    if not seq:
        return 0
    max_len = 1
    current_len = 1
    prev_base = seq[0]
    for base in seq[1:]:
        if base == prev_base:
            current_len += 1
            max_len = max(max_len, current_len)
        else:
            current_len = 1
            prev_base = base
    return max_len

def linguistic_complexity(seq: str) -> float:
    """
    Calculate the linguistic complexity of a nucleotide sequence.

    Linguistic complexity is defined as the ratio of observed k-mers to possible k-mers.

    Args:
        seq (str): DNA or RNA sequence.

    Returns:
        float: Linguistic complexity value between 0 and 1.
    """
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return 1.0
    observed = 0
    possible = 0
    for k in range(1, n + 1):
        max_possible_kmers = min(4 ** k, n - k + 1)
        observed_kmers = set(seq[i:i + k] for i in range(n - k + 1))
        observed += len(observed_kmers)
        possible += max_possible_kmers
    return observed / possible if possible > 0 else 1.0

def collapse_ambiguities(seq: str) -> str:
    """
    Collapse ambiguity codes in a nucleotide sequence to a single nucleotide by randomly selecting one possibility.

    Args:
        seq (str): DNA or RNA sequence with potential ambiguity codes.

    Returns:
        str: Sequence with ambiguity codes resolved to single nucleotides.
    """
    ambig_codes: Dict[str, List[str]] = {
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T'],
    }
    new_seq = []
    for base in seq.upper():
        if base in ambig_codes:
            new_base = random.choice(ambig_codes[base])
            new_seq.append(new_base)
        else:
            new_seq.append(base)
    return ''.join(new_seq)

def calculate_mfe(sequence: str) -> Union[float, str]:
    """
    Calculate the Minimum Free Energy (MFE) of a nucleotide sequence using the ViennaRNA package.

    Args:
        sequence (str): DNA sequence.

    Returns:
        Union[float, str]: MFE value in kcal/mol or 'NA' if calculation fails.
    """
    try:
        md = RNA.md()
        md.material = 'DNA'  # Set the model to use DNA parameters
        fc = RNA.fold_compound(sequence, md)
        structure, mfe = fc.mfe()
        return mfe
    except Exception as e:
        logging.debug(f"MFE calculation error for sequence: {e}")
        return 'NA'

def calculate_self_alignment_score(seq: str) -> Union[float, str]:
    """
    Calculate the self-alignment score between a sequence and its reverse complement.

    Args:
        seq (str): DNA or RNA sequence.

    Returns:
        Union[float, str]: Alignment score or 'NA' if calculation fails.
    """
    try:
        aligner = PairwiseAligner()
        # Set alignment parameters
        aligner.mode = 'global'
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -1

        rc_seq = str(Seq.Seq(seq).reverse_complement())
        score = aligner.score(seq, rc_seq)
        return score
    except Exception as e:
        logging.debug(f"Self-alignment score calculation error for sequence: {e}")
        return 'NA'

def calculate_entropy(seq: str) -> float:
    """
    Calculate the Shannon entropy of a nucleotide sequence.

    Args:
        seq (str): DNA or RNA sequence.

    Returns:
        float: Shannon entropy value.
    """
    seq = seq.upper()
    length = len(seq)
    if length == 0:
        return 0.0
    # Count the occurrences of each nucleotide
    freq: Dict[str, float] = {base: 0.0 for base in ['A', 'C', 'G', 'T']}
    for base in ['A', 'C', 'G', 'T']:
        freq[base] = seq.count(base) / length
    # Calculate entropy
    entropy = 0.0
    for p in freq.values():
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy

def open_fasta(filename: str) -> IO[str]:
    """
    Open a FASTA or FASTA.GZ file and return a file handle.
    Automatically detects if the file is gzipped based on the file extension.

    Args:
        filename (str): Path to the input FASTA or FASTA.GZ file.

    Returns:
        IO[str]: File handle opened in text mode.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')  # Text mode
    else:
        return open(filename, 'r')

def filter_bait(seq: str, args: argparse.Namespace) -> Tuple[bool, Dict[str, Any]]:
    """
    Apply all filters to a bait sequence and collect its parameters.

    Args:
        seq (str): Nucleotide sequence of the bait.
        args (argparse.Namespace): Parsed command-line arguments.

    Returns:
        Tuple[bool, Dict[str, Any]]:
            - bool: Whether the sequence passes all filters.
            - Dict[str, Any]: Dictionary of calculated parameters.
    """
    passes_filter = True
    # Remove unexpected characters
    seq = seq.replace(' ', '').replace('\n', '').replace('\r', '')
    valid_bases = set('ACGTUacgtuRYSWKMBDHVN-')
    if not set(seq).issubset(valid_bases):
        logging.debug(f"Sequence contains invalid characters: {seq}")
        passes_filter = False
    # Handle gaps
    if args.gaps == 'exclude' and '-' in seq:
        passes_filter = False
    elif args.gaps == 'extend' and '-' in seq:
        seq = seq.replace('-', '')
        if len(seq) < args.length:
            # Extend by repeating the sequence itself
            repeats_needed = (args.length - len(seq)) // len(seq) + 1
            seq = (seq * repeats_needed)[:args.length]
    # Complete bait length check
    if args.complete and len(seq) < args.length:
        passes_filter = False
    # Exclude Ns
    if args.noNs and 'N' in seq.upper():
        passes_filter = False
    # Check for zero-length sequence after processing
    if len(seq) == 0:
        logging.debug("Sequence length is zero after processing.")
        passes_filter = False
    # Compute GC content
    try:
        gc_content = gc_fraction(seq)  # Fraction between 0 and 1
        gc_percent = gc_content * 100
    except Exception as e:
        logging.debug(f"GC calculation error for sequence: {e}")
        gc_content = 0.0
        gc_percent = 0.0
    # Melting temperature calculation
    try:
        tm_value = mt.Tm_NN(seq, nn_table=mt.DNA_NN1, dnac1=0.25, dnac2=0.25)  # Adjusted DNA concentration
    except Exception as e:
        logging.debug(f"Tm calculation error for sequence: {e}")
        tm_value = 'NA'
    # Masked sequence check (lowercase letters)
    num_masked = sum(1 for c in seq if c.islower())
    percent_masked = (num_masked / len(seq)) * 100 if len(seq) > 0 else 0
    # Homopolymer length
    max_homopoly = max_homopolymer_length(seq)
    # Linguistic complexity
    seq_for_lc = collapse_ambiguities(seq)
    seq_for_lc = seq_for_lc.replace('-', '').upper()
    seq_complexity = linguistic_complexity(seq_for_lc)
    # Minimum Free Energy (MFE) Calculation
    mfe_value = calculate_mfe(seq)
    # Self-alignment score
    self_score = calculate_self_alignment_score(seq)
    # Shannon entropy calculation
    entropy_value = calculate_entropy(seq)
    # Ns and Gaps
    ns_present = 'true' if 'N' in seq.upper() else 'false'
    gaps_present = 'true' if '-' in seq else 'false'
    # Collect parameters
    params: Dict[str, Any] = {
        'GC%': round(gc_percent, 2) if isinstance(gc_percent, float) else 'NA',
        'Tm': round(tm_value, 2) if isinstance(tm_value, float) else 'NA',
        'Masked%': round(percent_masked, 2),
        'MaxHomopolymer': max_homopoly,
        'SeqComplexity': round(seq_complexity, 12) if isinstance(seq_complexity, float) else 'NA',
        'MFE': round(mfe_value, 2) if isinstance(mfe_value, float) else 'NA',
        'SelfScore': round(self_score, 2) if isinstance(self_score, float) else 'NA',
        'Entropy': round(entropy_value, 4),
        'Ns': ns_present,
        'Gaps': gaps_present
    }
    # Apply filters based on computed parameters
    # GC content check
    if args.mingc is not None and gc_percent < args.mingc:
        passes_filter = False
    if args.maxgc is not None and gc_percent > args.maxgc:
        passes_filter = False
    # Melting temperature filters
    if isinstance(tm_value, float):
        if args.mint is not None and tm_value < args.mint:
            passes_filter = False
        if args.maxt is not None and tm_value > args.maxt:
            passes_filter = False
    # MFE filtering
    if isinstance(mfe_value, float):
        if args.minMFE is not None and mfe_value < args.minMFE:
            passes_filter = False
        if args.maxMFE is not None and mfe_value > args.maxMFE:
            passes_filter = False
    # Self-alignment score filtering
    if isinstance(self_score, float) and args.maxSelfScore is not None and self_score > args.maxSelfScore:
        passes_filter = False
    # Entropy filtering
    if args.minEntropy is not None and entropy_value < args.minEntropy:
        passes_filter = False
    # Masked sequence check
    if args.maxmask is not None and percent_masked > args.maxmask:
        passes_filter = False
    # Homopolymer length check
    if args.maxhomopoly is not None and max_homopoly > args.maxhomopoly:
        passes_filter = False
    # Linguistic complexity check
    if args.minlc is not None and seq_complexity < args.minlc:
        passes_filter = False
    return passes_filter, params

def process_sequences(seq_chunk: List[SeqIO.SeqRecord], thread_id: int, args: argparse.Namespace) -> None:
    """
    Process a chunk of sequences: apply filters and collect parameters.

    Args:
        seq_chunk (List[SeqIO.SeqRecord]): List of sequence records to process.
        thread_id (int): Identifier for the current thread.
        args (argparse.Namespace): Parsed command-line arguments.
    """
    filtered_sequences: List[SeqIO.SeqRecord] = []
    all_params: List[str] = []
    for record in seq_chunk:
        seq_id: str = record.id
        seq: str = str(record.seq)
        try:
            # Sanitize the sequence
            seq = seq.replace(' ', '').replace('\n', '').replace('\r', '')
            if args.rc:
                seq = str(Seq.Seq(seq).reverse_complement())
            if args.collapse:
                seq = collapse_ambiguities(seq)
            # Adjust length
            bait_length: int = args.length
            if len(seq) > bait_length:
                seq = seq[:bait_length]
            # Convert to RNA if needed
            if args.rna:
                seq = seq.replace('T', 'U').replace('t', 'u')
            # Apply filters and collect parameters
            passes_filter, params = filter_bait(seq, args)
            # Collect parameters for all baits
            params.update({
                'Bait': seq_id,
                'BaitLength': len(seq),
                'Kept': 'Yes' if passes_filter else 'No'
            })
            param_line: str = '\t'.join(str(params.get(key, 'NA')) for key in [
                'Bait', 'BaitLength', 'GC%', 'Tm', 'Masked%', 'MaxHomopolymer',
                'SeqComplexity', 'MFE', 'SelfScore', 'Entropy', 'Ns', 'Gaps', 'Kept'
            ])
            all_params.append(param_line)
            # Collect filtered sequences
            if passes_filter:
                # Write sequence
                filtered_sequences.append(SeqIO.SeqRecord(Seq.Seq(seq), id=seq_id, description=''))
        except Exception as e:
            logging.error(f"Thread {thread_id}: Error processing sequence {seq_id}: {e}")
            # Collect minimal parameters if possible
            params: Dict[str, Any] = {
                'Bait': seq_id,
                'BaitLength': len(seq) if seq else 'NA',
                'GC%': 'NA',
                'Tm': 'NA',
                'Masked%': 'NA',
                'MaxHomopolymer': 'NA',
                'SeqComplexity': 'NA',
                'MFE': 'NA',
                'SelfScore': 'NA',
                'Entropy': 'NA',
                'Ns': 'NA',
                'Gaps': 'NA',
                'Kept': 'No'
            }
            param_line = '\t'.join(str(params.get(key, 'NA')) for key in [
                'Bait', 'BaitLength', 'GC%', 'Tm', 'Masked%', 'MaxHomopolymer',
                'SeqComplexity', 'MFE', 'SelfScore', 'Entropy', 'Ns', 'Gaps', 'Kept'
            ])
            all_params.append(param_line)
            continue  # Skip to the next sequence
    # Write to temporary files within the output directory if '--filter' is set
    if args.filter:
        temp_seq_file = os.path.join(args.outdir, f"{args.outprefix}-filtered-baits.fa.thread{thread_id}")
        SeqIO.write(filtered_sequences, temp_seq_file, 'fasta')
    # Always write parameters
    temp_params_file = os.path.join(args.outdir, f"{args.outprefix}-filtered-params.txt.thread{thread_id}")
    with open(temp_params_file, 'w') as pf:
        header = '\t'.join([
            'Bait', 'BaitLength', 'GC%', 'Tm', 'Masked%', 'MaxHomopolymer',
            'SeqComplexity', 'MFE', 'SelfScore', 'Entropy', 'Ns', 'Gaps', 'Kept'
        ])
        pf.write(header + '\n')
        for line in all_params:
            pf.write(line + '\n')

def process_all_sequences(seq_list: List[SeqIO.SeqRecord], args: argparse.Namespace) -> None:
    """
    Distribute sequences across threads, process them, and merge results.

    Args:
        seq_list (List[SeqIO.SeqRecord]): List of all sequence records to process.
        args (argparse.Namespace): Parsed command-line arguments.
    """
    # Split sequences into chunks for each thread
    sequences_per_thread: List[List[SeqIO.SeqRecord]] = []
    k, m = divmod(len(seq_list), args.threads)
    for i in range(args.threads):
        start = i * k + min(i, m)
        end = (i + 1) * k + min(i + 1, m)
        sequences_per_thread.append(seq_list[start:end])
    # Log the number of sequences per thread
    for i, chunk in enumerate(sequences_per_thread):
        logging.debug(f"Thread {i} processing {len(chunk)} sequences")
    # Start multiprocessing
    processes: List[multiprocessing.Process] = []
    for i in range(args.threads):
        p = multiprocessing.Process(target=process_sequences, args=(sequences_per_thread[i], i, args))
        processes.append(p)
        p.start()
    for p in processes:
        p.join()
    # Concatenate output files
    if args.filter:
        final_seq_file: str = os.path.join(args.outdir, f"{args.outprefix}-filtered-baits.fa")
        with open(final_seq_file, 'w') as out_f:
            for i in range(args.threads):
                temp_file = os.path.join(args.outdir, f"{args.outprefix}-filtered-baits.fa.thread{i}")
                if os.path.exists(temp_file):
                    with open(temp_file, 'r') as temp_f:
                        out_f.write(temp_f.read())
                    os.remove(temp_file)
    # Always concatenate parameters files
    final_params_file: str = os.path.join(args.outdir, f"{args.outprefix}-filtered-params.txt")
    with open(final_params_file, 'w') as out_f:
        header = '\t'.join([
            'Bait', 'BaitLength', 'GC%', 'Tm', 'Masked%', 'MaxHomopolymer',
            'SeqComplexity', 'MFE', 'SelfScore', 'Entropy', 'Ns', 'Gaps', 'Kept'
        ])
        out_f.write(header + '\n')
        for i in range(args.threads):
            temp_file = os.path.join(args.outdir, f"{args.outprefix}-filtered-params.txt.thread{i}")
            if os.path.exists(temp_file):
                with open(temp_file, 'r') as temp_f:
                    lines = temp_f.readlines()
                    out_f.writelines(lines[1:])  # Skip header
                os.remove(temp_file)
    # Now, count kept and filtered out baits by reading the params file
    kept = 0
    filtered = 0
    try:
        with open(final_params_file, 'r') as pf:
            header_line = pf.readline().strip().split('\t')
            try:
                kept_index = header_line.index('Kept')
            except ValueError:
                logging.error("'Kept' column not found in parameters file.")
                kept_index = -1
            if kept_index == -1:
                logging.error("Cannot determine 'Kept' status without the 'Kept' column.")
            else:
                for line in pf:
                    fields = line.strip().split('\t')
                    if len(fields) <= kept_index:
                        logging.warning(f"Malformed line skipped: {line.strip()}")
                        continue  # Skip malformed lines
                    kept_status = fields[kept_index]
                    if kept_status == 'Yes':
                        kept += 1
                    elif kept_status == 'No':
                        filtered += 1
    except FileNotFoundError:
        logging.error(f"Parameters file {final_params_file} not found.")
    # Log total baits processed
    total_baits_processed = kept + filtered
    logging.info(f"Total baits processed: {total_baits_processed}")
    logging.info(f"Baits kept: {kept}")
    logging.info(f"Baits filtered out: {filtered}")
    # Verify total baits match input sequences
    if total_baits_processed != len(seq_list):
        logging.warning(f"Total baits processed ({total_baits_processed}) does not match total input sequences ({len(seq_list)}).")
    # Gzip output files if requested
    if args.gzip:
        files_to_gzip: List[str] = []
        if args.filter:
            files_to_gzip.append(final_seq_file)
        files_to_gzip.append(final_params_file)
        for filename in files_to_gzip:
            if filename and os.path.exists(filename):
                with open(filename, 'rb') as f_in, gzip.open(filename + '.gz', 'wb') as f_out:
                    f_out.writelines(f_in)
                os.remove(filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate bait statistics from a FASTA file.')
    add_arguments(parser)
    args = parser.parse_args()
    main(args)