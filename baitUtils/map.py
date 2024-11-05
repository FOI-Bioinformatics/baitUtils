#!/usr/bin/env python3

"""
map.py

Map baits against target (whitelist) and off-target (blacklist) sequences.

Features:
- Map baits using pblat or mmseqs2
- Generate a parameters file similar to the one from stats.py
- Filter baits that map better than a cutoff percentage (default 50%) against blacklist
- Output filtered baits and parameters file
"""

from baitUtils._version import __version__
import argparse
import logging
import os
import sys
import subprocess
import shutil
from typing import List, Dict, Any
from Bio import SeqIO

def add_arguments(parser):
    """
    Add command-line arguments specific to the 'map' subcommand.

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
        '-w', '--whitelist',
        required=True,
        help='Target sequences FASTA file (whitelist)'
    )
    parser.add_argument(
        '-b', '--blacklist',
        required=True,
        help='Off-target sequences FASTA file (blacklist)'
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
        choices=['pblat', 'mmseqs2'],
        default='pblat',
        help='Mapping tool to use (Default: pblat)'
    )
    parser.add_argument(
        '--cutoff',
        type=float,
        default=50.0,
        help='Cutoff percentage to filter baits that map better than this against blacklist (Default: 50.0)'
    )
    parser.add_argument(
        '-X', '--threads',
        type=int,
        default=1,
        help='Number of threads (Default: 1)'
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

    # Create output directory if it doesn't exist
    create_output_directory(args.outdir)

    # Map baits against whitelist
    logging.info("Mapping baits against whitelist...")
    whitelist_mappings = map_baits(args.input, args.whitelist, args.mapper, args.threads, 'whitelist', args.outdir, args.outprefix)

    # Map baits against blacklist
    logging.info("Mapping baits against blacklist...")
    blacklist_mappings = map_baits(args.input, args.blacklist, args.mapper, args.threads, 'blacklist', args.outdir, args.outprefix)

    # Parse mappings and generate parameters
    logging.info("Parsing mappings and generating parameters...")
    params = parse_mappings(whitelist_mappings, blacklist_mappings, args.cutoff, args.mapper)

    # Output parameters file
    params_file = os.path.join(args.outdir, f"{args.outprefix}-mapping-params.txt")
    write_params(params, params_file)

    # Filter baits based on cutoff
    logging.info("Filtering baits based on cutoff...")
    filtered_baits = filter_baits(params, args.cutoff)

    # Output filtered baits
    filtered_baits_file = os.path.join(args.outdir, f"{args.outprefix}-filtered-baits.fa")
    write_filtered_baits(filtered_baits, args.input, filtered_baits_file)

    logging.info("Mapping and filtering complete.")

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

def map_baits(baits_file: str, target_file: str, mapper: str, threads: int, label: str, outdir: str, outprefix: str) -> str:
    """
    Map baits against target sequences using specified mapping tool.

    Args:
        baits_file (str): Path to baits fasta file.
        target_file (str): Path to target fasta file.
        mapper (str): Mapping tool to use ('pblat' or 'mmseqs2').
        threads (int): Number of threads to use.
        label (str): Label for the mapping ('whitelist' or 'blacklist').
        outdir (str): Output directory.
        outprefix (str): Output file prefix.

    Returns:
        str: Path to mapping output file.
    """
    # Determine output file path
    mapping_output = os.path.join(outdir, f"{outprefix}-{label}-mapping.tsv")

    if mapper == 'pblat':
        run_pblat(baits_file, target_file, threads, mapping_output)
    elif mapper == 'mmseqs2':
        run_mmseqs2(baits_file, target_file, threads, mapping_output)
    else:
        logging.error(f"Unsupported mapper: {mapper}")
        sys.exit(1)

    return mapping_output

def run_pblat(baits_file: str, target_file: str, threads: int, output_file: str) -> None:
    """
    Run pblat to map baits against target sequences.

    Args:
        baits_file (str): Path to baits fasta file.
        target_file (str): Path to target fasta file.
        threads (int): Number of threads to use.
        output_file (str): Path to output file.
    """
    # pblat options
    # Build command
    cmd = [
        'pblat',
        '-threads={}'.format(threads),
        target_file,
        baits_file,
        output_file
    ]

    # Run command
    logging.info(f"Running pblat: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        logging.error(f"Error running pblat: {e}")
        sys.exit(1)

def run_mmseqs2(baits_file: str, target_file: str, threads: int, output_file: str) -> None:
    """
    Run mmseqs2 to map baits against target sequences.

    Args:
        baits_file (str): Path to baits fasta file.
        target_file (str): Path to target fasta file.
        threads (int): Number of threads to use.
        output_file (str): Path to output file.
    """
    # mmseqs2 commands
    # Create temporary directories
    tmp_dir = os.path.join(os.path.dirname(output_file), f'mmseqs_tmp_{os.getpid()}')
    os.makedirs(tmp_dir, exist_ok=True)

    # Define database and result paths
    query_db = os.path.join(tmp_dir, 'queryDB')
    target_db = os.path.join(tmp_dir, 'targetDB')
    result_db = os.path.join(tmp_dir, 'resultDB')

    # Build and run commands
    try:
        # Create query database
        cmd = ['mmseqs', 'createdb', baits_file, query_db]
        logging.info(f"Running mmseqs createdb for query: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

        # Create target database
        cmd = ['mmseqs', 'createdb', target_file, target_db]
        logging.info(f"Running mmseqs createdb for target: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

        # Run search
        cmd = ['mmseqs', 'search', query_db, target_db, result_db, tmp_dir, '--threads', str(threads)]
        logging.info(f"Running mmseqs search: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

        # Convert alignment results to tsv
        cmd = ['mmseqs', 'convertalis', query_db, target_db, result_db, output_file, '--format-output', 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits']
        logging.info(f"Running mmseqs convertalis: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    except Exception as e:
        logging.error(f"Error running mmseqs2: {e}")
        sys.exit(1)
    finally:
        # Clean up temporary directory
        shutil.rmtree(tmp_dir, ignore_errors=True)

def parse_mappings(whitelist_mappings: str, blacklist_mappings: str, cutoff: float, mapper: str) -> Dict[str, Dict[str, Any]]:
    """
    Parse mappings and generate parameters for each bait.

    Args:
        whitelist_mappings (str): Path to whitelist mapping output file.
        blacklist_mappings (str): Path to blacklist mapping output file.
        cutoff (float): Cutoff percentage to filter baits.
        mapper (str): Mapping tool used ('pblat' or 'mmseqs2')

    Returns:
        Dict[str, Dict[str, Any]]: Dictionary of parameters for each bait.
    """
    # Initialize parameters dictionary
    params: Dict[str, Dict[str, Any]] = {}

    # Parse whitelist mappings
    logging.info("Parsing whitelist mappings...")
    if mapper == 'pblat':
        whitelist_hits = parse_psl_file(whitelist_mappings)
    else:
        whitelist_hits = parse_mapping_file(whitelist_mappings)

    # Parse blacklist mappings
    logging.info("Parsing blacklist mappings...")
    if mapper == 'pblat':
        blacklist_hits = parse_psl_file(blacklist_mappings)
    else:
        blacklist_hits = parse_mapping_file(blacklist_mappings)

    # For each bait, collect mapping info
    all_baits = set(whitelist_hits.keys()).union(blacklist_hits.keys())

    for bait in all_baits:
        params[bait] = {}
        # Get the best hit against whitelist
        whitelist_best = whitelist_hits.get(bait, {'pident': 0.0})
        # Get the best hit against blacklist
        blacklist_best = blacklist_hits.get(bait, {'pident': 0.0})

        params[bait]['Whitelist_Pident'] = whitelist_best['pident']
        params[bait]['Blacklist_Pident'] = blacklist_best['pident']
        # Determine whether to keep or filter the bait
        if blacklist_best['pident'] >= cutoff:
            params[bait]['Kept'] = 'No'
        else:
            params[bait]['Kept'] = 'Yes'

    return params

def parse_mapping_file(mapping_file: str) -> Dict[str, Dict[str, Any]]:
    """
    Parse mapping output file and return best hit per bait.

    Args:
        mapping_file (str): Path to mapping output file.

    Returns:
        Dict[str, Dict[str, Any]]: Dictionary with bait IDs as keys and best hit info as values.
    """
    hits: Dict[str, Dict[str, Any]] = {}
    try:
        with open(mapping_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue  # Skip header or comments
                cols = line.strip().split('\t')
                if len(cols) < 3:
                    continue  # Skip incomplete lines
                query_id = cols[0]
                target_id = cols[1]
                pident = float(cols[2])

                # Update best hit if pident is higher
                if query_id not in hits or pident > hits[query_id]['pident']:
                    hits[query_id] = {
                        'target_id': target_id,
                        'pident': pident
                    }
    except Exception as e:
        logging.error(f"Error parsing mapping file {mapping_file}: {e}")
        sys.exit(1)

    return hits

def parse_psl_file(psl_file: str) -> Dict[str, Dict[str, Any]]:
    """
    Parse PSL file and return best hit per bait.

    Args:
        psl_file (str): Path to PSL file.

    Returns:
        Dict[str, Dict[str, Any]]: Dictionary with bait IDs as keys and best hit info as values.
    """
    hits: Dict[str, Dict[str, Any]] = {}
    try:
        with open(psl_file, 'r') as f:
            for line in f:
                if line.startswith('psLayout') or line.startswith('---'):
                    continue  # Skip header and separator lines
                cols = line.strip().split('\t')
                if len(cols) < 21:
                    continue  # Skip incomplete lines

                matches = int(cols[0])
                misMatches = int(cols[1])
                repMatches = int(cols[2])
                qBaseInsert = int(cols[5])
                qName = cols[9]
                tName = cols[13]

                # Compute alignment length
                alignment_length = matches + misMatches + qBaseInsert

                # Compute percent identity
                if alignment_length > 0:
                    pident = (matches + repMatches) / alignment_length * 100.0
                else:
                    pident = 0.0

                # Update best hit if pident is higher
                if qName not in hits or pident > hits[qName]['pident']:
                    hits[qName] = {
                        'target_id': tName,
                        'pident': pident
                    }
    except Exception as e:
        logging.error(f"Error parsing PSL file {psl_file}: {e}")
        sys.exit(1)

    return hits

def write_params(params: Dict[str, Dict[str, Any]], params_file: str) -> None:
    """
    Write parameters to file.

    Args:
        params (Dict[str, Dict[str, Any]]): Dictionary of parameters.
        params_file (str): Path to parameters file.
    """
    with open(params_file, 'w') as f:
        # Write header
        f.write('\t'.join(['Bait', 'Whitelist_Pident', 'Blacklist_Pident', 'Kept']) + '\n')
        for bait, param_dict in params.items():
            line = '\t'.join([
                bait,
                str(round(param_dict.get('Whitelist_Pident', 0.0), 2)),
                str(round(param_dict.get('Blacklist_Pident', 0.0), 2)),
                param_dict.get('Kept', 'No')
            ])
            f.write(line + '\n')

def filter_baits(params: Dict[str, Dict[str, Any]], cutoff: float) -> List[str]:
    """
    Filter baits based on blacklist percent identity cutoff.

    Args:
        params (Dict[str, Dict[str, Any]]): Dictionary of parameters.
        cutoff (float): Cutoff percentage.

    Returns:
        List[str]: List of bait IDs that are kept.
    """
    filtered_baits = [bait for bait, param_dict in params.items() if param_dict.get('Kept') == 'Yes']
    return filtered_baits

def write_filtered_baits(filtered_baits: List[str], input_fasta: str, output_fasta: str) -> None:
    """
    Write filtered baits to fasta file.

    Args:
        filtered_baits (List[str]): List of bait IDs to keep.
        input_fasta (str): Path to input baits fasta file.
        output_fasta (str): Path to output fasta file.
    """
    # Read input fasta file
    seq_records = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))

    # Write filtered sequences
    with open(output_fasta, 'w') as f:
        for bait_id in filtered_baits:
            if bait_id in seq_records:
                SeqIO.write(seq_records[bait_id], f, 'fasta')
            else:
                logging.warning(f"Bait ID {bait_id} not found in input fasta.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Map baits against target and off-target sequences.')
    add_arguments(parser)
    args = parser.parse_args()
    main(args)