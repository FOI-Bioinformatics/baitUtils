#!/usr/bin/env python3

"""
check.py

Subcommand: check
Evaluate coverage of a set of oligos (PSL-derived BED) and report uncovered regions.
Optionally write uncovered stretches above a threshold, and export them as FASTA.

Usage:
    baitUtils check --psl input.psl [options]
"""

import sys
import math
import logging
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict
from dataclasses import dataclass

import pybedtools
from pybedtools import BedTool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ---------------------------
# Utility Functions
# ---------------------------

def setup_logging(log_level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S"
    )

def write_uncovered_regions(
    uncovered_regions: List[Tuple[str, int, int, float]], 
    output_path: Path,
    length_cutoff: int = 0
) -> int:
    """
    Write uncovered regions above a certain length cutoff to file.
    Returns the number of regions written.
    """
    filtered_regions = []
    for ref, s, e, cov in uncovered_regions:
        length = e - s
        if length >= length_cutoff:
            filtered_regions.append((ref, s, e, cov, length))
    
    filtered_regions.sort(key=lambda x: x[4], reverse=True)
    
    with open(output_path, "w") as f:
        f.write("Reference\tStart\tEnd\tCoverage\tLength\n")
        for ref, s, e, cov, length in filtered_regions:
            f.write(f"{ref}\t{s}\t{e}\t{cov}\t{length}\n")
    return len(filtered_regions)

def split_sequence_at_n(sequence: str) -> List[Tuple[int, str]]:
    """Split a sequence at 'N' bases. Return list of (start_offset, subsequence)."""
    parts = []
    current_start = 0
    current_seq = []
    
    for i, base in enumerate(sequence):
        if base.upper() == 'N':
            if current_seq:
                parts.append((current_start, ''.join(current_seq)))
                current_seq = []
            current_start = i + 1
        else:
            current_seq.append(base)
    
    if current_seq:
        parts.append((current_start, ''.join(current_seq)))
    
    return parts

def export_uncovered_fasta(
    uncovered_regions: Dict[str, List[Tuple[int, int]]],
    fasta_path: Path,
    output_path: Path,
    n_split_path: Path,
    extend_bp: int = 0,
    length_cutoff: int = 0,
    oligo_length: int = 120
):
    """
    Export uncovered regions as FASTA (optionally extended on both sides),
    and split regions at 'N' into smaller segments.
    """
    from Bio import SeqIO
    reference_seqs = {record.id: record.seq 
                      for record in SeqIO.parse(fasta_path, "fasta")}
    
    uncovered_records = []
    n_split_records = []
    total_regions = 0
    n_split_regions = 0
    
    for ref_id, regions in uncovered_regions.items():
        if ref_id not in reference_seqs:
            logging.warning(f"Reference {ref_id} not found in FASTA.")
            continue
            
        ref_seq = reference_seqs[ref_id]
        for i, (start, end) in enumerate(regions):
            if end - start < length_cutoff:
                continue
                
            ext_start = max(0, start - extend_bp)
            ext_end = min(len(ref_seq), end + extend_bp)
            
            seq = str(ref_seq[ext_start:ext_end])
            total_regions += 1
            
            base_record = SeqRecord(
                Seq(seq),
                id=f"{ref_id}_uncovered_{i+1}",
                description=f"pos={ext_start}-{ext_end} original={start}-{end}"
            )
            uncovered_records.append(base_record)
            
            split_parts = split_sequence_at_n(seq)
            for j, (offset, subseq) in enumerate(split_parts):
                if len(subseq) >= oligo_length:
                    n_split_regions += 1
                    split_record = SeqRecord(
                        Seq(subseq),
                        id=f"{ref_id}_uncovered_{i+1}_split_{j+1}",
                        description=f"pos={ext_start+offset}-{ext_start+offset+len(subseq)} original={start}-{end}"
                    )
                    n_split_records.append(split_record)
    
    SeqIO.write(uncovered_records, output_path, "fasta")
    SeqIO.write(n_split_records, n_split_path, "fasta")
    
    logging.info(f"Wrote {total_regions} uncovered regions to {output_path}")
    logging.info(f"Wrote {n_split_regions} N-split regions to {n_split_path}")

def read_forced_oligos(path: Path | None):
    if not path:
        return set()
    try:
        with open(path) as f:
            return {line.strip().upper() for line in f if line.strip()}
    except Exception as e:
        logging.error(f"Error reading forced oligos: {e}")
        sys.exit(1)

def parse_psl_to_bed(
    psl_path: Path,
    min_length: int,
    min_similarity: float,
    temp_dir: Path | None = None
) -> BedTool:
    """
    Parse PSL file into BED format (keeping only hits above min_length and min_similarity).
    """
    from tqdm import tqdm
    from collections import defaultdict

    def parse_psl_lines():
        with open(psl_path) as f:
            for line in tqdm(f, desc="Parsing PSL", unit=" lines"):
                if line.startswith(("psLayout", "match", "-", "#")):
                    continue
                parts = line.strip().split()
                if len(parts) < 17:
                    continue
                try:
                    matches = int(parts[0])
                    qName = parts[9].strip().upper()
                    qSize = float(parts[10])
                    ref_id = parts[13]
                    tStart = int(parts[15])
                    tEnd = int(parts[16])
                    if (tEnd - tStart) < min_length:
                        continue
                    if (matches / qSize * 100.0) < min_similarity:
                        continue
                    yield f"{ref_id}\t{tStart}\t{tEnd}\t{qName}\t1.0\n"
                except (ValueError, IndexError):
                    continue

    temp_bed = "temp_check.bed"
    if temp_dir:
        temp_bed = str(temp_dir / "temp_check.bed")

    with open(temp_bed, "w") as f:
        for record in parse_psl_lines():
            f.write(record)
    return BedTool(temp_bed).sort()

def create_genome_file(bed: BedTool, temp_dir: Path | None = None) -> str:
    """
    Build a genome file from the sorted BED intervals 
    (required by bedtools genomecov).
    """
    from collections import defaultdict
    chrom_sizes = defaultdict(int)
    for interval in bed:
        chrom_sizes[interval.chrom] = max(chrom_sizes[interval.chrom], interval.end)

    genome_path = "genome_check.txt" if temp_dir is None else str(temp_dir / "genome_check.txt")
    with open(genome_path, "w") as gf:
        for chrom, size in chrom_sizes.items():
            gf.write(f"{chrom}\t{size}\n")
    return genome_path

def compute_uncovered_regions(
    coverage_bed: BedTool,
    min_coverage: float,
    max_coverage: float | None,
    genome_file: str
) -> Tuple[Dict[str, List[Tuple[int,int]]], int, BedTool]:
    """
    Compute uncovered regions below min_coverage using bedtools genomecov output.
    Returns a dictionary of uncovered intervals by chromosome, total uncovered bases,
    and the coverage bedtool object (for optional saving).
    """
    coverage = coverage_bed.genomecov(bga=True, g=genome_file)

    uncovered_regions = defaultdict(list)
    total_uncovered = 0

    curr_chrom = None
    curr_start = None
    prev_end = None

    for interval in coverage:
        chrom, start, end, cov = interval
        if chrom == "genome":
            continue
        cov = float(cov)
        start = int(start)
        end = int(end)
        length = end - start

        if cov < min_coverage:
            total_uncovered += length
            if curr_chrom != chrom:
                if curr_chrom is not None and curr_start is not None:
                    uncovered_regions[curr_chrom].append((curr_start, prev_end))
                curr_chrom = chrom
                curr_start = start
                prev_end = end
            else:
                if curr_start is None:
                    curr_start = start
                    prev_end = end
                else:
                    if prev_end == start:
                        prev_end = end
                    else:
                        uncovered_regions[curr_chrom].append((curr_start, prev_end))
                        curr_start = start
                        prev_end = end
        else:
            # Coverage >= min_coverage: close out any active uncovered region
            if curr_chrom == chrom and curr_start is not None:
                uncovered_regions[curr_chrom].append((curr_start, prev_end))
            curr_chrom = None
            curr_start = None
            prev_end = None

    # If we ended on an uncovered segment, finalize it
    if curr_chrom is not None and curr_start is not None:
        uncovered_regions[curr_chrom].append((curr_start, prev_end))

    return uncovered_regions, total_uncovered, coverage

def check_coverage_mode(
    bed: BedTool,
    forced_oligos: set,
    min_coverage: float,
    max_coverage: float | None,
    coverage_out: Path | None,
    longest_uncovered_out: Path | None,
    temp_dir: Path | None,
    uncovered_length_cutoff: int,
    args
):
    """
    Evaluate coverage of forced oligos (if any), or the entire BED if none forced,
    then write optional coverage data and uncovered intervals.
    """
    # Filter BED to forced oligos if forced_oligos is non-empty
    if forced_oligos:
        forced_bed_entries = []
        for interval in bed:
            if interval.name in forced_oligos:
                forced_bed_entries.append(str(interval))
        if not forced_bed_entries:
            logging.error("No forced oligos matched the BED. Exiting.")
            sys.exit(1)
        tmp_file = "forced_temp_check.bed"
        if temp_dir:
            tmp_file = str(temp_dir / "forced_temp_check.bed")
        with open(tmp_file, "w") as f:
            f.write("\n".join(forced_bed_entries) + "\n")
        bed_filtered = BedTool(tmp_file).sort()
    else:
        bed_filtered = bed

    # Create genome file & compute coverage
    genome_file = create_genome_file(bed_filtered, temp_dir)
    uncovered_regions, total_uncovered, coverage = compute_uncovered_regions(
        bed_filtered,
        min_coverage,
        max_coverage,
        genome_file
    )

    # Summarize coverage stats
    total_bases = sum(
        int(line.strip().split()[1])
        for line in open(genome_file)
    )
    if total_bases > 0:
        logging.info(f"Total bases: {total_bases:,}")
        logging.info(f"Bases under min coverage: {total_uncovered:,} ({total_uncovered/total_bases*100:.1f}%)")

    # Optionally write coverage out
    if coverage_out:
        with open(coverage_out, "w") as f:
            f.write("Reference\tStart\tEnd\tCoverage\n")
            for interval in coverage:
                if interval.chrom != "genome":
                    f.write(f"{interval.chrom}\t{interval.start}\t{interval.end}\t{interval.name}\n")
        logging.info(f"Coverage data written to {coverage_out}")

    # Optionally write uncovered intervals
    if longest_uncovered_out:
        zero_cov_list = []
        for chrom, intervals in uncovered_regions.items():
            for (s, e) in intervals:
                zero_cov_list.append((chrom, s, e, 0.0))
        count = write_uncovered_regions(zero_cov_list, longest_uncovered_out, uncovered_length_cutoff)
        logging.info(f"Wrote {count} uncovered stretches >= {uncovered_length_cutoff}bp to {longest_uncovered_out}")

    # Optionally export uncovered FASTA
    if args.uncovered_fasta and args.fasta_reference:
        if not args.n_split_fasta:
            logging.warning("--n_split_fasta not specified, skipping N-split output")
        else:
            export_uncovered_fasta(
                uncovered_regions,
                args.fasta_reference,
                args.uncovered_fasta,
                args.n_split_fasta,
                args.extend_region,
                uncovered_length_cutoff,
                args.min_oligo_length
            )

# ---------------------------
# Subcommand Interface
# ---------------------------

def add_arguments(parser):
    parser.add_argument(
        "--psl",
        type=Path,
        required=True,
        help="Path to PSL-like file"
    )
    parser.add_argument(
        "--min_coverage",
        type=float,
        default=10.0,
        help="Minimum coverage required per base (default=10.0)"
    )
    parser.add_argument(
        "--max_coverage",
        type=float,
        help="Maximum coverage allowed per base"
    )
    parser.add_argument(
        "--forced_oligos",
        type=Path,
        help="File with oligo IDs that must be included in coverage check"
    )
    parser.add_argument(
        "--coverage_out",
        type=Path,
        help="Output file for coverage data"
    )
    parser.add_argument(
        "--longest_uncovered_out",
        type=Path,
        default=Path("longest_uncovered.txt"),
        help="Output file for uncovered regions"
    )
    parser.add_argument(
        "--log_level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set logging level"
    )
    parser.add_argument(
        "--temp_dir",
        type=Path,
        help="Directory for temporary files"
    )
    parser.add_argument(
        "--uncovered_length_cutoff",
        type=int,
        default=0,
        help="Save uncovered stretches longer than this cutoff (default=0)"
    )
    parser.add_argument(
        "--min_length",
        type=int,
        default=100,
        help="Minimum mapping length to consider (default=100)"
    )
    parser.add_argument(
        "--min_similarity",
        type=float,
        default=95.0,
        help="Minimum percent identity to consider (default=95.0)"
    )
    # FASTA options for uncovered regions
    parser.add_argument(
        "--fasta_reference",
        type=Path,
        help="Reference FASTA file for exporting uncovered regions"
    )
    parser.add_argument(
        "--uncovered_fasta",
        type=Path,
        help="Output FASTA file for uncovered regions"
    )
    parser.add_argument(
        "--n_split_fasta",
        type=Path,
        help="Output FASTA file for N-split regions"
    )
    parser.add_argument(
        "--extend_region",
        type=int,
        default=0,
        help="Number of bases to extend uncovered regions on each side"
    )
    parser.add_argument(
        "--min_oligo_length",
        type=int,
        default=120,
        help="Minimum sub-region length after N-splitting (default=120)"
    )

def main(args):
    setup_logging(args.log_level)

    forced_oligos = read_forced_oligos(args.forced_oligos)
    
    bed = parse_psl_to_bed(
        args.psl,
        args.min_length,
        args.min_similarity,
        args.temp_dir
    )
    if bed.count() == 0:
        logging.error("No intervals found in PSL after filtering. Exiting.")
        sys.exit(1)

    check_coverage_mode(
        bed,
        forced_oligos,
        args.min_coverage,
        args.max_coverage,
        args.coverage_out,
        args.longest_uncovered_out,
        args.temp_dir,
        args.uncovered_length_cutoff,
        args
    )