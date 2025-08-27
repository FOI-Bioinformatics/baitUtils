#!/usr/bin/env python3

"""
fill.py

Subcommand: fill
Multi-pass selection of oligos to maximize coverage.
"""

import sys
import math
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Set
from collections import defaultdict
from dataclasses import dataclass

try:
    import pybedtools
    from pybedtools import BedTool
    HAS_PYBEDTOOLS = True
except ImportError:
    HAS_PYBEDTOOLS = False
    BedTool = None
    logging.warning("pybedtools not available - gap filling functionality may be limited")
from Bio import SeqIO
from Bio.Seq import Seq

@dataclass
class OligoMapping:
    """Represents a single oligo mapping."""
    oligo_id: str
    ref_id: str
    start: int
    end: int
    coverage: float = 1.0

# ----------------------------------------------------------------
# Argument parser setup
# ----------------------------------------------------------------
def add_arguments(parser):
    parser.add_argument(
        "--psl",
        type=Path,
        required=True,
        help="Path to PSL-like file"
    )
    parser.add_argument(
        "--min_coverage",
        type=int,
        default=1,
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
        help="File with oligo IDs that must be included"
    )
    parser.add_argument(
        "--spacing_distance",
        type=int,
        default=30,
        help="Min distance between start positions of selected oligos (default=30)"
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
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("selected_oligos.txt"),
        help="Output file path for selected oligos"
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
        "--min_contribution",
        type=int,
        default=5,
        help="Minimum coverage contribution for selection (default=5)"
    )
    parser.add_argument(
        "--max_passes",
        type=int,
        default=5,
        help="Maximum number of selection passes (default=5)"
    )
    parser.add_argument(
        "--uncovered_length_cutoff",
        type=int,
        default=0,
        help="Save uncovered stretches longer than this cutoff (default=0)"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Continue selection even if uncovered regions don't decrease"
    )
    parser.add_argument(
        "--stall_rounds",
        type=int,
        default=3,
        help="Number of rounds without improvement before stopping (default=3)"
    )
    parser.add_argument(
        "--max_oligos_per_pass",
        type=int,
        help="Maximum number of oligos to select in each pass"
    )
    parser.add_argument(
        "--reference_sequence",
        type=Path,
        help="Reference sequence file for sequence-based scoring"
    )
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
        help="Number of bases to extend uncovered regions by on each side"
    )
    parser.add_argument(
        "--min_oligo_length",
        type=int,
        default=120,
        help="Minimum oligo length for N-split regions (default=120)"
    )

# ----------------------------------------------------------------
# The main logic
# (the rest of your multi-pass selection functions remain the same)
# ----------------------------------------------------------------

def setup_logging(log_level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S"
    )

def read_forced_oligos(path: Path | None) -> Set[str]:
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

    temp_bed = "temp_fill.bed"
    if temp_dir:
        temp_bed = str(temp_dir / "temp_fill.bed")

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

    genome_path = "genome_fill.txt" if temp_dir is None else str(temp_dir / "genome_fill.txt")
    with open(genome_path, "w") as gf:
        for chrom, size in chrom_sizes.items():
            gf.write(f"{chrom}\t{size}\n")
    return genome_path

def build_mappings(bed: BedTool) -> Dict[str, List[OligoMapping]]:
    """Build dictionary of oligo mappings from a sorted BedTool."""
    ref_dict = defaultdict(list)
    for interval in bed:
        ref_id = interval.chrom
        start = interval.start
        end = interval.end
        oligo_id = interval.name
        ref_dict[ref_id].append(OligoMapping(oligo_id, ref_id, start, end, 1.0))
    
    # Sort each list by start position
    for ref_id in ref_dict:
        ref_dict[ref_id].sort(key=lambda x: x.start)
    return ref_dict

def coverage_from_selected(
    selected_oligos: Set[str],
    all_mappings: Dict[str, List[OligoMapping]],
    genome_file: str,
    temp_dir: Path | None = None
) -> BedTool:
    """Build BedTool from selected oligos for coverage calculation."""
    temp_bed = "selected_temp_fill.bed"
    if temp_dir:
        temp_bed = str(temp_dir / "selected_temp_fill.bed")

    with open(temp_bed, "w") as outbed:
        for ref_id, mapping_list in all_mappings.items():
            for m in mapping_list:
                if m.oligo_id in selected_oligos:
                    outbed.write(f"{ref_id}\t{m.start}\t{m.end}\t{m.oligo_id}\t1.0\n")

    return BedTool(temp_bed).sort()

def compute_uncovered_regions(
    coverage_bed: BedTool,
    min_coverage: float,
    max_coverage: float | None,
    genome_file: str
) -> Tuple[Dict[str, List[Tuple[int,int]]], int, BedTool]:
    """
    Compute uncovered regions below min_coverage using bedtools genomecov output.
    Returns dictionary of uncovered intervals, total uncovered bases, and the coverage tool.
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
            if curr_chrom == chrom and curr_start is not None:
                uncovered_regions[curr_chrom].append((curr_start, prev_end))
            curr_chrom = None
            curr_start = None
            prev_end = None

    if curr_chrom is not None and curr_start is not None:
        uncovered_regions[curr_chrom].append((curr_start, prev_end))

    return uncovered_regions, total_uncovered, coverage

# ---------------------------
# Multi-pass Coverage Selection
# ---------------------------

def analyze_coverage_patterns(
    uncovered: Dict[str, List[Tuple[int,int]]],
    all_mappings: Dict[str, List[OligoMapping]]
) -> Dict[str, float]:
    """Analyze coverage patterns to identify difficult or under-covered regions."""
    difficulty_scores = {}
    coverage_options = defaultdict(int)
    
    for ref_id, mappings in all_mappings.items():
        for mapping in mappings:
            for ustart, uend in uncovered.get(ref_id, []):
                if min(mapping.end, uend) - max(mapping.start, ustart) > 0:
                    key = (ref_id, ustart, uend)
                    coverage_options[key] += 1
    
    for (ref_id, ustart, uend), options in coverage_options.items():
        difficulty = 1.0 / (1.0 + math.log(1 + options))
        region_key = f"{ref_id}:{ustart}-{uend}"
        difficulty_scores[region_key] = difficulty
    
    return difficulty_scores

def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of a sequence."""
    if not sequence:
        return 0.0
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)

def score_oligo(
    m: OligoMapping,
    uncovered: Dict[str, List[Tuple[int,int]]],
    min_contribution: int,
    unique_region_bonus: float = 2.0,
    coverage_threshold: float = 0.8,
    difficulty_scores: Dict[str, float] = None,
    sequence: str = None,
    gc_weight: float = 0.2,
    optimal_gc: float = 0.5,
    length_weight: float = 1.5
) -> float:
    """
    Enhanced oligo scoring that considers uncovered region lengths, coverage efficiency,
    region difficulty, sequence-based GC content, and unique coverage bonuses.
    """
    overlap_regions = []
    total_overlap = 0
    difficulty_bonus = 0
    max_region_length = 0
    length_weighted_overlap = 0
    
    # Calculate overlaps
    for (ustart, uend) in uncovered.get(m.ref_id, []):
        region_length = uend - ustart
        overlap = min(m.end, uend) - max(m.start, ustart)
        
        if overlap > 0:
            length_weighted_overlap += overlap * math.log2(region_length)
            max_region_length = max(max_region_length, region_length)
            overlap_regions.append((overlap, region_length))
            total_overlap += overlap
            
            if difficulty_scores:
                region_key = f"{m.ref_id}:{ustart}-{uend}"
                difficulty_bonus += overlap * difficulty_scores.get(region_key, 0)
    
    if total_overlap < min_contribution:
        return 0.0
    
    oligo_length = m.end - m.start
    coverage_efficiency = total_overlap / oligo_length
    length_bonus = math.log2(max_region_length) * length_weight
    
    base_score = (total_overlap + length_weighted_overlap) * (1.0 + coverage_efficiency + length_bonus)

    if difficulty_bonus > 0:
        base_score *= (1.0 + difficulty_bonus)
    
    # Sequence-based GC scoring
    if sequence:
        gc_content = calculate_gc_content(sequence)
        gc_penalty = abs(gc_content - optimal_gc)
        sequence_score = 1.0 - (gc_weight * gc_penalty)
        base_score *= sequence_score
    
    # Unique coverage bonus
    unique_coverage_bonus = 0
    for overlap, region_length in overlap_regions:
        if region_length < 2 * oligo_length:
            coverage_ratio = overlap / region_length
            if coverage_ratio > coverage_threshold:
                bonus_scale = coverage_ratio * unique_region_bonus * math.log2(region_length)
                unique_coverage_bonus += overlap * bonus_scale

    return base_score + unique_coverage_bonus

def can_select_oligo(
    selected_positions: Dict[str, List[Tuple[int,int]]],
    candidate: OligoMapping,
    spacing: int
) -> bool:
    """Check if candidate oligo respects the required spacing from previously selected oligos."""
    import bisect
    positions = selected_positions.get(candidate.ref_id, [])
    if not positions:
        return True

    idx = bisect.bisect_left(positions, (candidate.start, candidate.end))

    if idx > 0:
        prev_start, _ = positions[idx - 1]
        if abs(candidate.start - prev_start) < spacing:
            return False

    if idx < len(positions):
        next_start, _ = positions[idx]
        if abs(next_start - candidate.start) < spacing:
            return False

    return True

def single_pass_greedy(
    all_mappings: Dict[str, List[OligoMapping]],
    uncovered: Dict[str, List[Tuple[int,int]]],
    selected_oligos: Set[str],
    spacing_distance: int,
    min_contribution: int,
    sequences: Dict[str, str] = None,
    max_oligos_per_pass: int = None
) -> Set[str]:
    """Perform a single pass of greedy selection based on scoring."""
    import bisect
    from collections import defaultdict

    selected_positions = defaultdict(list)
    for ref_id, mlist in all_mappings.items():
        for m in mlist:
            if m.oligo_id in selected_oligos:
                bisect.insort(selected_positions[ref_id], (m.start, m.end))

    difficulty_scores = analyze_coverage_patterns(uncovered, all_mappings)

    candidate_list = []
    for ref_id, mlist in all_mappings.items():
        for m in mlist:
            if m.oligo_id in selected_oligos:
                continue

            seq_slice = None
            if sequences and ref_id in sequences:
                seq_slice = sequences[ref_id][m.start:m.end]

            sc = score_oligo(
                m,
                uncovered,
                min_contribution,
                difficulty_scores=difficulty_scores,
                sequence=seq_slice
            )

            if sc > 0:
                candidate_list.append((sc, m))

    candidate_list.sort(key=lambda x: x[0], reverse=True)

    new_selected = set()
    changes = False

    for sc, cand_oligo in candidate_list:
        if max_oligos_per_pass and len(new_selected) >= max_oligos_per_pass:
            break

        if can_select_oligo(selected_positions, cand_oligo, spacing_distance):
            new_selected.add(cand_oligo.oligo_id)
            bisect.insort(
                selected_positions[cand_oligo.ref_id],
                (cand_oligo.start, cand_oligo.end)
            )
            changes = True

    if changes:
        logging.debug(f"Selected {len(new_selected)} new oligos in this pass")
    else:
        logging.debug("No new oligos selected in this pass")

    return new_selected

def multi_pass_selection(
    all_mappings: Dict[str, List[OligoMapping]],
    forced_oligos: Set[str],
    genome_file: str,
    min_coverage: float,
    max_coverage: float | None,
    spacing_distance: int,
    min_contribution: int,
    max_passes: int,
    max_oligos_per_pass: int = None,
    sequences: Dict[str, str] = None,
    force: bool = False,
    uncovered_length_cutoff: int = 0,
    stall_rounds: int = 3
) -> Set[str]:
    """
    Perform multi-pass coverage selection until coverage 
    stops improving or max passes is reached.
    """
    all_oligo_ids = set()
    for ref_id, mapping_list in all_mappings.items():
        for m in mapping_list:
            all_oligo_ids.add(m.oligo_id)

    selected_oligos = set(forced_oligos).intersection(all_oligo_ids)
    prev_uncovered = float('inf')
    prev_uncovered_count = float('inf')
    stall_count = 0

    pass_num = 1
    while pass_num <= max_passes:
        logging.info(f"Multi-pass iteration {pass_num}/{max_passes}")

        coverage_bed = coverage_from_selected(selected_oligos, all_mappings, genome_file)
        uncovered_regions, total_uncovered, _ = compute_uncovered_regions(
            coverage_bed,
            min_coverage,
            max_coverage,
            genome_file
        )
        
        uncovered_count = sum(
            1 for regions in uncovered_regions.values()
            for start, end in regions
            if end - start >= uncovered_length_cutoff
        )
        
        logging.info(f"Uncovered bases: {total_uncovered:,}")
        logging.info(f"Uncovered regions >= {uncovered_length_cutoff}bp: {uncovered_count}")

        # Check improvement in uncovered region count
        if not force:
            if uncovered_count >= prev_uncovered_count:
                stall_count += 1
                logging.info(f"No improvement in uncovered region count. Stall count: {stall_count}/{stall_rounds}")
                if stall_count >= stall_rounds:
                    logging.info(f"Stopping after {stall_rounds} rounds without improvement.")
                    break
            else:
                stall_count = 0
            
        if total_uncovered >= prev_uncovered:
            logging.info("No improvement in uncovered base count. Stopping.")
            break
            
        prev_uncovered = total_uncovered
        prev_uncovered_count = uncovered_count

        newly_selected = single_pass_greedy(
            all_mappings,
            uncovered_regions,
            selected_oligos,
            spacing_distance,
            min_contribution,
            sequences,
            max_oligos_per_pass
        )

        if not newly_selected:
            logging.info("No new oligos selected; coverage can't be improved further.")
            break

        old_count = len(selected_oligos)
        selected_oligos.update(newly_selected)
        new_count = len(selected_oligos)
        logging.info(f"Pass {pass_num}: selected {new_count - old_count} new oligos. Total now: {new_count}")

        pass_num += 1

    return selected_oligos

def calculate_coverage(
    bed: BedTool,
    min_coverage: float,
    max_coverage: float | None,
    temp_dir: Path | None,
    genome_file: str
) -> Tuple[
    Dict[str, List[Tuple[int,int]]], 
    int, 
    int, 
    int, 
    List[Tuple[str, int, int, float]], 
    Dict[str,int]
]:
    """
    Compute coverage stats (using bedtools genomecov). 
    Return uncovered regions, bases under/over coverage, total bases, 
    the raw coverage data, and chrom sizes.
    """
    chrom_sizes = {}
    with open(genome_file) as gf:
        total_bases = 0
        for line in gf:
            chrom, size = line.strip().split()
            size = int(size)
            chrom_sizes[chrom] = size
            total_bases += size

    coverage = bed.genomecov(bga=True, g=genome_file)
    uncovered_regions = defaultdict(list)
    coverage_data = []
    bases_under = 0
    bases_over = 0
    curr_chrom = None
    curr_start = None
    prev_end = None

    for interval in coverage:
        chrom, start, end, cov = interval
        if chrom == "genome":
            continue
        start = int(start)
        end = int(end)
        cov = float(cov)
        length = end - start
        coverage_data.append((chrom, start, end, cov))
        
        if cov < min_coverage:
            bases_under += length
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
            if max_coverage and cov > max_coverage:
                bases_over += length
            if curr_chrom == chrom and curr_start is not None:
                uncovered_regions[curr_chrom].append((curr_start, prev_end))
            curr_chrom = None
            curr_start = None
            prev_end = None

    if curr_chrom is not None and curr_start is not None:
        uncovered_regions[curr_chrom].append((curr_start, prev_end))

    return uncovered_regions, bases_under, bases_over, total_bases, coverage_data, chrom_sizes

def main(args):
    """Main entry point for fill subcommand."""
    setup_logging(args.log_level)

    forced_oligos = read_forced_oligos(args.forced_oligos)
    bed = parse_psl_to_bed(
        args.psl,
        args.min_length,
        args.min_similarity,
        args.temp_dir
    )
    count = bed.count()
    if count == 0:
        logging.error("No intervals found in PSL after filtering. Exiting.")
        sys.exit(1)

    logging.info(f"Total oligos for selection: {count}")

    mappings_dict = build_mappings(bed)
    genome_file = create_genome_file(bed, args.temp_dir)

    # Optional: load reference sequences if provided
    sequences = None
    if args.reference_sequence:
        sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.reference_sequence, "fasta")}

    # Multi-pass selection
    selected_oligos = multi_pass_selection(
        mappings_dict,
        forced_oligos,
        genome_file,
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

    # Write out final selections
    final_count = len(selected_oligos)
    logging.info(f"Final selection contains {final_count} oligos.")
    
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        for oligo_id in sorted(selected_oligos):
            f.write(f"{oligo_id}\n")
    logging.info(f"Selected oligos written to {args.output}")

    # Optional coverage check of final selection
    coverage_bed = coverage_from_selected(selected_oligos, mappings_dict, genome_file, args.temp_dir)
    uncovered, bases_under, bases_over, total_bases, coverage_data, _ = calculate_coverage(
        coverage_bed,
        args.min_coverage,
        args.max_coverage,
        args.temp_dir,
        genome_file
    )
    if total_bases > 0:
        logging.info(f"Final coverage check: Bases under min coverage: {bases_under:,} "
                     f"({bases_under/total_bases*100:.1f}%)")

    # Optionally write coverage data
    if args.coverage_out:
        with open(args.coverage_out, "w") as f:
            f.write("Reference\tStart\tEnd\tCoverage\n")
            for ref, start, end, cov in coverage_data:
                f.write(f"{ref}\t{start}\t{end}\t{cov}\n")
        logging.info(f"Coverage data written to {args.coverage_out}")

    # Optionally write uncovered region stats
    if args.longest_uncovered_out:
        from baitUtils.check import write_uncovered_regions  # Reuse from check.py or define locally
        uncovered_0 = [
            (ref, s, e, 0.0) 
            for ref, intervals in uncovered.items() 
            for (s, e) in intervals
        ]
        count_uncovered = write_uncovered_regions(uncovered_0, args.longest_uncovered_out, args.uncovered_length_cutoff)
        logging.info(f"Wrote {count_uncovered} uncovered stretches >= {args.uncovered_length_cutoff}bp to {args.longest_uncovered_out}")

    # Optionally export uncovered FASTA
    if args.uncovered_fasta and args.fasta_reference:
        if not args.n_split_fasta:
            logging.warning("--n_split_fasta not specified, skipping N-split output")
        else:
            from baitUtils.check import export_uncovered_fasta  # Reuse from check.py or define locally
            export_uncovered_fasta(
                uncovered,
                args.fasta_reference,
                args.uncovered_fasta,
                args.n_split_fasta,
                args.extend_region,
                args.uncovered_length_cutoff,
                args.min_oligo_length
            )