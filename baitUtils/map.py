#!/usr/bin/env python3

"""
map.py

Map baits against a genome sequence using pblat.

Features:
- Map baits using pblat
- Output mapped/unmapped baits and their IDs
- Allow customization of pblat parameters
- Filter mappings based on identity percentage
- Enhanced logging with counts of probes and sequences
"""

from baitUtils._version import __version__
import argparse
import logging
import os
import sys
import subprocess
from typing import Set, Dict
from Bio import SeqIO

def add_arguments(parser):
    """
    Add command-line arguments for the 'map' subcommand.

    Args:
        parser (argparse.ArgumentParser): The argument parser to which arguments will be added.
    """
    parser.add_argument('--version', action='version', version=f"baitUtils map {__version__}")
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input baits FASTA file'
    )
    parser.add_argument(
        '-q', '--query',
        required=True,
        help='FASTA genome file to map against'
    )
    parser.add_argument(
        '-o', '--outprefix',
        default='out',
        help='Output file prefix (Default = out)'
    )
    parser.add_argument(
        '-Z', '--outdir',
        default='.',
        help='Output directory path (Default is ./)'
    )
    parser.add_argument(
        '--mapper',
        choices=['pblat'],
        default='pblat',
        help='Mapping tool to use (only pblat is currently supported)'
    )
    parser.add_argument(
        '-X', '--threads',
        type=int,
        default=1,
        help='Number of threads (Default: 1)'
    )
    parser.add_argument(
        '--minMatch',
        type=int,
        default=2,
        help='Sets the number of tile matches. Usually set from 2 to 4. Default is 2 for nucleotide.'
    )
    parser.add_argument(
        '--minScore',
        type=int,
        default=30,
        help='Sets minimum score. Default is 30.'
    )
    parser.add_argument(
        '--minIdentity',
        type=int,
        default=90,
        help='Sets minimum sequence identity (in percent). Default is 90 for nucleotide searches.',
        choices=range(0, 101)  # Restrict to 0-100
    )
    parser.add_argument(
        '--filterIdentity',
        type=int,
        default=90,
        help='Filter mappings with identity percentage less than this value. Must be greater than or equal to minIdentity.',
        choices=range(0, 101)  # Restrict to 0-100
    )
    parser.add_argument(
        '--minMatchCount',  # Using minMatchCount to avoid confusion with pblat's minMatch
        type=int,
        default=0,
        help='Minimum number of matching bases required (corresponds to first column in PSL). Default is 0 (no filtering).'
    )
    parser.add_argument(
        '--fasta-output',
        choices=['mapped', 'unmapped', 'both', 'none'],
        default='mapped',
        help='Which probes to include in the FASTA output file: mapped, unmapped, both, or none. Default is mapped.'
    )
    parser.add_argument(
        '-l', '--log',
        action='store_true',
        help='Enable detailed logging'
    )

def main(args):
    """
    Main function for the 'map' subcommand.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    # Set up logging
    setup_logging(args.log)

    # Check that filterIdentity >= minIdentity
    if args.filterIdentity < args.minIdentity:
        logging.error(f"filterIdentity ({args.filterIdentity}) must be greater than or equal to minIdentity ({args.minIdentity}).")
        sys.exit(1)

    # Create output directory if it doesn't exist
    create_output_directory(args.outdir)

    # Read all bait IDs from input FASTA file
    logging.info("Reading input FASTA file...")
    try:
        seq_records = SeqIO.to_dict(SeqIO.parse(args.input, 'fasta'))
        all_bait_ids = set(seq_records.keys())
        logging.info(f"Loaded {len(all_bait_ids)} probes from input FASTA file.")
    except Exception as e:
        logging.error(f"Error reading input FASTA file {args.input}: {e}")
        sys.exit(1)

    # Count number of sequences in genome FASTA file
    logging.info("Reading genome FASTA file...")
    try:
        genome_seq_count = sum(1 for _ in SeqIO.parse(args.query, 'fasta'))
        logging.info(f"Genome FASTA file contains {genome_seq_count} sequences.")
    except Exception as e:
        logging.error(f"Error reading genome FASTA file {args.query}: {e}")
        sys.exit(1)

    # Map baits against query genome
    logging.info("Mapping baits against query genome...")
    mapping_output = map_baits(
        args.input,
        args.query,
        args.mapper,
        args.threads,
        args.outdir,
        args.outprefix,
        args.minMatch,
        args.minScore,
        args.minIdentity
    )

    # Determine filtered output path
    mapping_filtered_output = os.path.join(args.outdir, f"{args.outprefix}-mapping_filtered.psl")
    
    # Parse mapping output to get baits that mapped
    logging.info("Parsing mapping output...")
    mapped_baits = parse_psl_file(mapping_output, args.filterIdentity, args.minMatchCount, 
                                filtered_output=mapping_filtered_output)

    # Compute number of mapped baits
    logging.info(f"Number of mapped probes: {len(mapped_baits)}")

    # Compute unmapped baits
    unmapped_baits = all_bait_ids - mapped_baits
    logging.info(f"Number of unmapped probes: {len(unmapped_baits)}")

    # Write text files with mapped and unmapped bait IDs
    mapped_ids_file = os.path.join(args.outdir, f"{args.outprefix}-mapped-bait-ids.txt")
    unmapped_ids_file = os.path.join(args.outdir, f"{args.outprefix}-unmapped-bait-ids.txt")
    write_bait_ids(mapped_baits, mapped_ids_file)
    write_bait_ids(unmapped_baits, unmapped_ids_file)

    # Output selected baits to FASTA file
    if args.fasta_output != 'none':
        if args.fasta_output == 'mapped':
            selected_baits = mapped_baits
            fasta_file = os.path.join(args.outdir, f"{args.outprefix}-mapped-baits.fa")
            logging.info(f"Writing mapped baits to FASTA file: {fasta_file}")
            write_selected_baits(selected_baits, seq_records, fasta_file)
        elif args.fasta_output == 'unmapped':
            selected_baits = unmapped_baits
            fasta_file = os.path.join(args.outdir, f"{args.outprefix}-unmapped-baits.fa")
            logging.info(f"Writing unmapped baits to FASTA file: {fasta_file}")
            write_selected_baits(selected_baits, seq_records, fasta_file)
        elif args.fasta_output == 'both':
            # Write mapped baits
            selected_baits = mapped_baits
            fasta_file_mapped = os.path.join(args.outdir, f"{args.outprefix}-mapped-baits.fa")
            logging.info(f"Writing mapped baits to FASTA file: {fasta_file_mapped}")
            write_selected_baits(selected_baits, seq_records, fasta_file_mapped)
            # Write unmapped baits
            selected_baits = unmapped_baits
            fasta_file_unmapped = os.path.join(args.outdir, f"{args.outprefix}-unmapped-baits.fa")
            logging.info(f"Writing unmapped baits to FASTA file: {fasta_file_unmapped}")
            write_selected_baits(selected_baits, seq_records, fasta_file_unmapped)
        else:
            logging.error(f"Invalid option for --fasta-output: {args.fasta_output}")
            sys.exit(1)

    logging.info("Mapping complete.")

def setup_logging(enable_debug: bool) -> None:
    """
    Configure logging settings.

    Args:
        enable_debug (bool): If True, set logging level to DEBUG. Otherwise, INFO.
    """
    log_level = logging.DEBUG if enable_debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s %(levelname)s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def create_output_directory(directory: str) -> None:
    """
    Create the output directory if it doesn't exist.

    Args:
        directory (str): Path to the output directory.
    """
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
            logging.info(f"Created output directory: {directory}")
        except Exception as e:
            logging.error(f"Failed to create output directory '{directory}': {e}")
            sys.exit(1)
    else:
        logging.info(f"Output directory already exists: {directory}")

def map_baits(baits_file: str, target_file: str, mapper: str, threads: int, outdir: str, outprefix: str, minMatch: int, minScore: int, minIdentity: int) -> str:
    """
    Map baits against target sequences using pblat.

    Args:
        baits_file (str): Path to baits FASTA file.
        target_file (str): Path to genome FASTA file.
        mapper (str): Mapping tool to use ('pblat').
        threads (int): Number of threads to use.
        outdir (str): Output directory.
        outprefix (str): Output file prefix.
        minMatch (int): Minimum number of tile matches.
        minScore (int): Minimum score.
        minIdentity (int): Minimum sequence identity (percent).

    Returns:
        str: Path to mapping output file.
    """
    # Determine output file path
    mapping_output = os.path.join(outdir, f"{outprefix}-mapping.psl")

    if mapper == 'pblat':
        run_pblat(baits_file, target_file, threads, mapping_output, minMatch, minScore, minIdentity)
    else:
        logging.error(f"Unsupported mapper: {mapper}")
        sys.exit(1)

    return mapping_output

def run_pblat(baits_file: str, target_file: str, threads: int, output_file: str, minMatch: int, minScore: int, minIdentity: int) -> None:
    """
    Run pblat to map baits against the genome sequence.

    Args:
        baits_file (str): Path to baits FASTA file.
        target_file (str): Path to genome FASTA file.
        threads (int): Number of threads to use.
        output_file (str): Path to output PSL file.
        minMatch (int): Minimum number of tile matches.
        minScore (int): Minimum score.
        minIdentity (int): Minimum sequence identity (percent).
    """
    # Build command
    cmd = [
        'pblat',
        f'-threads={threads}',
        f'-minMatch={minMatch}',
        f'-minScore={minScore}',
        f'-minIdentity={minIdentity}',
        target_file,
        baits_file,
        output_file
    ]

    # Run command with stdout and stderr suppressed
    logging.info(f"Running pblat: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        logging.error(f"pblat failed with return code {e.returncode}.")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error running pblat: {e}")
        sys.exit(1)

def parse_psl_file(psl_file: str, filterIdentity: int, minMatchCount: int = 0, filtered_output: str = None) -> Set[str]:
    """
    Parse PSL file and return set of bait IDs that had hits with identity >= filterIdentity
    and matches >= minMatchCount. Optionally save filtered hits to a new PSL file.

    Args:
        psl_file (str): Path to PSL file.
        filterIdentity (int): Identity percentage threshold for filtering.
        minMatchCount (int): Minimum number of matching bases required (Default: 0).
        filtered_output (str): Optional path to save filtered PSL entries.

    Returns:
        Set[str]: Set of bait IDs that had hits meeting both thresholds.
    """
    mapped_baits = set()
    header_lines = []
    filtered_hits = []
    
    try:
        with open(psl_file, 'r') as f:
            for line in f:
                line = line.strip()
                # Store header lines if we're going to write filtered output
                if filtered_output and (line.startswith('psLayout') or line.startswith('-') or line.startswith('match')):
                    header_lines.append(line)
                    continue
                if not line or line.startswith('psLayout') or line.startswith('-') or line.startswith('match') or line.startswith('no matches'):
                    continue  # Skip header and separator lines
                
                # Split on whitespace
                cols = line.split()
                if len(cols) < 21:
                    continue  # Skip incomplete lines
                try:
                    matches = int(cols[0])
                    misMatches = int(cols[1])
                    repMatches = int(cols[2])
                    nCount = int(cols[3])
                    qNumInsert = int(cols[4])
                    qBaseInsert = int(cols[5])
                    tNumInsert = int(cols[6])
                    tBaseInsert = int(cols[7])
                    qName = cols[9]
                except ValueError:
                    logging.warning(f"Non-integer value encountered in PSL file: {line}")
                    continue

                # Check both filters
                if matches < minMatchCount:
                    continue

                size = matches + misMatches + repMatches + qNumInsert
                if size == 0:
                    continue
                    
                identity_percentage = 100.0 * (matches + repMatches) / size
                if identity_percentage >= filterIdentity:
                    mapped_baits.add(qName)
                    if filtered_output:
                        filtered_hits.append(line)

        # Write filtered hits if output path provided
        if filtered_output and filtered_hits:
            try:
                with open(filtered_output, 'w') as f:
                    # Write header
                    for header_line in header_lines:
                        f.write(header_line + '\n')
                    # Write filtered hits
                    for hit in filtered_hits:
                        f.write(hit + '\n')
                logging.info(f"Wrote {len(filtered_hits)} filtered hits to {filtered_output}")
            except Exception as e:
                logging.error(f"Error writing filtered PSL file {filtered_output}: {e}")
                
    except Exception as e:
        logging.error(f"Error parsing PSL file {psl_file}: {e}")
        sys.exit(1)
        
    return mapped_baits

def write_bait_ids(bait_ids: Set[str], output_file: str) -> None:
    """
    Write bait IDs to a text file, one per line.

    Args:
        bait_ids (Set[str]): Set of bait IDs to write.
        output_file (str): Path to output text file.
    """
    try:
        with open(output_file, 'w') as f:
            for bait_id in sorted(bait_ids):
                f.write(bait_id + '\n')
    except Exception as e:
        logging.error(f"Error writing bait IDs to file {output_file}: {e}")
        sys.exit(1)

def write_selected_baits(selected_baits: Set[str], seq_records: Dict[str, SeqIO.SeqRecord], output_fasta: str) -> None:
    """
    Write selected baits to a FASTA file.

    Args:
        selected_baits (Set[str]): Set of bait IDs to write.
        seq_records (dict): Dictionary of SeqRecord objects from input FASTA.
        output_fasta (str): Path to output FASTA file.
    """
    try:
        with open(output_fasta, 'w') as f:
            for bait_id in selected_baits:
                if bait_id in seq_records:
                    SeqIO.write(seq_records[bait_id], f, 'fasta')
                else:
                    logging.warning(f"Bait ID {bait_id} not found in input FASTA.")
    except Exception as e:
        logging.error(f"Error writing selected baits to FASTA file {output_fasta}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Map baits against a genome sequence using pblat.')
    add_arguments(parser)
    args = parser.parse_args()
    main(args)