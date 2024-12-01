#!/usr/bin/env python3

import argparse
import logging
import sys
from typing import List, Tuple
import matplotlib.pyplot as plt

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a coverage plot from a PSL file and optionally save sorted probe mappings to a BED file."
    )
    parser.add_argument(
        "psl_file", help="Input PSL file from pblat."
    )
    parser.add_argument(
        "-o", "--output", help="Output plot file (PNG format). If not specified, the plot will be displayed."
    )
    parser.add_argument(
        "-b", "--bed", help="Output BED file to save sorted probe mappings. If specified, probe mappings will be saved to this file."
    )
    parser.add_argument(
        "-l", "--log", help="Log file to write logs. Defaults to stderr."
    )
    return parser.parse_args()

def setup_logging(log_file: str = None) -> None:
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

def read_psl_file(psl_file: str) -> Tuple[str, int, List[Tuple[int, int]]]:
    """
    Reads the PSL file and extracts the target name, target length, and alignment positions.
    """
    target_name = ""
    target_length = 0
    alignments = []

    try:
        with open(psl_file, "r") as f:
            lines = f.readlines()
    except IOError as e:
        logging.error(f"Error reading PSL file: {e}")
        sys.exit(1)

    data_started = False
    for line in lines:
        line = line.strip()
        if line.startswith("psLayout"):
            continue  # Skip header lines
        if line.startswith("---"):
            data_started = True
            continue
        if not data_started or not line:
            continue  # Skip until data starts
        fields = line.split()
        if len(fields) < 17:
            continue  # Skip incomplete lines

        # Extract target name and length from the first alignment
        if not target_name:
            target_name = fields[13]
            target_length = int(fields[14])

        t_start = int(fields[15])
        t_end = int(fields[16])
        alignments.append((t_start, t_end))
        logging.debug(f"Alignment: t_start={t_start}, t_end={t_end}")

    if not target_name or target_length == 0:
        logging.error("Failed to extract target name or length from PSL file.")
        sys.exit(1)

    logging.info(f"Extracted target name: {target_name}")
    logging.info(f"Extracted target length: {target_length}")
    logging.info(f"Total alignments: {len(alignments)}")

    return target_name, target_length, alignments

def save_to_bed(bed_file: str, target_name: str, alignments: List[Tuple[int, int]]) -> None:
    """
    Saves the probe mappings to a BED file, sorted based on coordinates.
    """
    try:
        # Sort alignments based on t_start and t_end
        sorted_alignments = sorted(alignments, key=lambda x: (x[0], x[1]))
        with open(bed_file, "w") as f:
            for t_start, t_end in sorted_alignments:
                f.write(f"{target_name}\t{t_start}\t{t_end}\n")
        logging.info(f"Probe mappings saved to sorted BED file: {bed_file}")
    except IOError as e:
        logging.error(f"Error writing to BED file: {e}")
        sys.exit(1)

def generate_coverage_array(target_length: int, alignments: List[Tuple[int, int]]) -> List[int]:
    """
    Generates a coverage array for the target sequence based on alignments.
    """
    coverage = [0] * target_length
    for t_start, t_end in alignments:
        # Ensure indices are within bounds
        t_start = max(0, t_start)
        t_end = min(target_length, t_end)
        for pos in range(t_start, t_end):
            coverage[pos] += 1
    logging.info("Coverage array generated.")
    return coverage

def plot_coverage(
    coverage: List[int], target_name: str, output_file: str = None
) -> None:
    """
    Plots the coverage array.
    """
    positions = list(range(1, len(coverage) + 1))

    plt.figure(figsize=(12, 6))
    plt.plot(positions, coverage, color="blue")
    plt.xlabel("Position along target sequence (bp)")
    plt.ylabel("Coverage depth")
    plt.title(f"Probe Coverage along Target Sequence {target_name}")
    plt.grid(True)
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300)
        logging.info(f"Coverage plot saved to {output_file}")
    else:
        plt.show()
        logging.info("Coverage plot displayed.")

def main():
    args = parse_arguments()
    setup_logging(args.log)

    logging.info("Starting coverage plot generation.")

    target_name, target_length, alignments = read_psl_file(args.psl_file)

    # Save probe mappings to sorted BED file if --bed is specified
    if args.bed:
        save_to_bed(args.bed, target_name, alignments)

    coverage = generate_coverage_array(target_length, alignments)
    plot_coverage(coverage, target_name, args.output)

    logging.info("Coverage plot generation completed.")

if __name__ == "__main__":
    main()