#!/usr/bin/env python3

"""
coverage_analysis.py

Coverage analysis utilities for calculating coverage statistics, 
handling BED files, and analyzing uncovered regions.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional, Any
from collections import defaultdict

try:
    from pybedtools import BedTool
    HAS_PYBEDTOOLS = True
except ImportError:
    HAS_PYBEDTOOLS = False
    BedTool = None

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from baitUtils.gap_filling_algorithm import OligoMapping
from baitUtils.mapping_utils import parse_alignments, filter_hits


def merge_uncovered_intervals(
    intervals: Any, min_coverage: float
) -> Tuple[Dict[str, List[Tuple[int, int]]], int]:
    """
    Merge adjacent genomecov (bga) intervals whose depth is below min_coverage.

    intervals yields (chrom, start, end, depth) tuples sorted by chromosome and
    position, as bedtools genomecov -bga produces. Returns a dict of merged
    (start, end) runs per chromosome and the total number of uncovered bases.
    """
    uncovered_regions: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    total_uncovered = 0
    run_chrom = None
    run_start = None
    run_end = None

    def flush():
        if run_chrom is not None:
            uncovered_regions[run_chrom].append((run_start, run_end))

    for interval in intervals:
        chrom, start, end, cov = interval[0], int(interval[1]), int(interval[2]), float(interval[3])
        if chrom == "genome":
            continue
        if cov < min_coverage:
            total_uncovered += end - start
            if run_chrom == chrom and run_end == start:
                run_end = end
            else:
                flush()
                run_chrom, run_start, run_end = chrom, start, end
        else:
            flush()
            run_chrom = run_start = run_end = None

    flush()
    return dict(uncovered_regions), total_uncovered


class PSLParser:
    """Converts PSL hits into a sorted BedTool of aligned blocks."""
    
    @staticmethod
    def parse_psl_to_bed(
        psl_path: Path,
        min_length: int,
        min_similarity: float,
        temp_dir: Optional[Path] = None,
        basename: str = "temp_fill.bed"
    ) -> Any:
        """
        Write one BED row per aligned block for hits with aligned length >= min_length
        and BLAT identity >= min_similarity, and return the sorted BedTool.
        """
        temp_bed = basename if temp_dir is None else str(temp_dir / basename)
        kept = 0
        with open(temp_bed, "w") as f:
            for hit in filter_hits(parse_alignments(psl_path), min_similarity, min_length):
                kept += 1
                for start, end in hit.target_blocks:
                    f.write(f"{hit.t_name}\t{start}\t{end}\t{hit.q_name}\t{kept}\n")
        logging.info(f"Kept {kept} PSL hits after filtering")
        return BedTool(temp_bed).sort()


class GenomeFileHandler:
    """Handles genome file creation and management."""
    
    @staticmethod
    def create_genome_file(bed: Any, temp_dir: Optional[Path] = None,
                           reference_fasta: Optional[Path] = None) -> str:
        """
        Build a genome file (chromosome sizes) for bedtools genomecov.

        Sizes are read from the reference FASTA when given. Without it the
        size of each reference is the end of its last mapped interval, which
        hides unmapped tails and references with no hits.
        """
        chrom_sizes = defaultdict(int)
        if reference_fasta is not None:
            for record in SeqIO.parse(str(reference_fasta), "fasta"):
                chrom_sizes[record.id] = len(record.seq)
        else:
            logging.warning("No reference FASTA given; reference sizes are taken from the last mapped base")
            for interval in bed:
                chrom_sizes[interval.chrom] = max(chrom_sizes[interval.chrom], interval.end)

        genome_path = "genome_fill.txt" if temp_dir is None else str(temp_dir / "genome_fill.txt")
        with open(genome_path, "w") as gf:
            for chrom, size in chrom_sizes.items():
                gf.write(f"{chrom}\t{size}\n")
        return genome_path


class MappingBuilder:
    """Builds oligo mappings from BED files."""
    
    @staticmethod
    def build_mappings(bed: Any) -> Dict[str, List[OligoMapping]]:
        """
        Build one OligoMapping per hit from a sorted BedTool whose rows are
        aligned blocks. Rows of the same hit share name and score (hit index)
        and are grouped back into one mapping spanning all its blocks.
        """
        grouped: Dict[Tuple[str, str, str], List[Tuple[int, int]]] = defaultdict(list)
        for interval in bed:
            grouped[(interval.chrom, interval.name, str(interval.score))].append((interval.start, interval.end))
        ref_dict = defaultdict(list)
        for (ref_id, oligo_id, hit_id), blocks in grouped.items():
            blocks.sort()
            ref_dict[ref_id].append(OligoMapping(
                oligo_id, ref_id, blocks[0][0], blocks[-1][1], 1.0, hit_id=hit_id, blocks=blocks))
        for ref_id in ref_dict:
            ref_dict[ref_id].sort(key=lambda x: x.start)
        return ref_dict


class CoverageCalculator:
    """Calculates coverage statistics and identifies uncovered regions."""
    
    @staticmethod
    def coverage_from_selected(
        selected: Set,
        all_mappings: Dict[str, List[OligoMapping]],
        genome_file: str,
        temp_dir: Optional[Path] = None
    ) -> Any:
        """
        Build a BedTool from the selected mappings. selected may hold mapping
        keys (see OligoMapping.key) or plain oligo IDs; an ID selects every
        mapping of that oligo.
        """
        temp_bed = "selected_temp_fill.bed" if temp_dir is None else str(temp_dir / "selected_temp_fill.bed")
        with open(temp_bed, "w") as outbed:
            for ref_id, mapping_list in all_mappings.items():
                for m in mapping_list:
                    if m.key in selected or m.oligo_id in selected:
                        for start, end in m.blocks:
                            outbed.write(f"{ref_id}\t{start}\t{end}\t{m.oligo_id}\t{m.hit_id}\n")
        return BedTool(temp_bed).sort()

    @staticmethod
    def compute_uncovered_regions(
        coverage_bed: Any,
        min_coverage: float,
        genome_file: str
    ) -> Tuple[Dict[str, List[Tuple[int, int]]], int, Any]:
        """
        Compute uncovered regions below min_coverage using bedtools genomecov output.
        Returns dictionary of uncovered intervals, total uncovered bases, and the coverage tool.
        """
        coverage = coverage_bed.genomecov(bga=True, g=genome_file)
        uncovered_regions, total_uncovered = merge_uncovered_intervals(coverage, min_coverage)
        return uncovered_regions, total_uncovered, coverage

    @staticmethod
    def calculate_coverage(
        bed: Any,
        min_coverage: float,
        genome_file: str,
        max_coverage: Optional[float] = None
    ) -> Tuple[Dict[str, List[Tuple[int, int]]], int, int, int, List[Tuple[str, int, int, float]], Dict[str, int]]:
        """
        Compute coverage with bedtools genomecov. Returns uncovered regions,
        bases under min_coverage, bases over max_coverage, total bases, the raw
        (chrom, start, end, depth) intervals and the chromosome sizes.
        """
        chrom_sizes: Dict[str, int] = {}
        with open(genome_file) as gf:
            for line in gf:
                chrom, size = line.strip().split()
                chrom_sizes[chrom] = int(size)
        total_bases = sum(chrom_sizes.values())

        coverage = bed.genomecov(bga=True, g=genome_file)
        coverage_data = [(iv[0], int(iv[1]), int(iv[2]), float(iv[3])) for iv in coverage if iv[0] != "genome"]
        uncovered_regions, bases_under = merge_uncovered_intervals(coverage_data, min_coverage)
        bases_over = 0
        if max_coverage is not None:
            bases_over = sum(e - st for _, st, e, cov in coverage_data if cov > max_coverage)
        return uncovered_regions, bases_under, bases_over, total_bases, coverage_data, chrom_sizes


def write_uncovered_regions(
    uncovered_regions: Dict[str, List[Tuple[int, int]]],
    output_path: Path,
    length_cutoff: int = 0
) -> int:
    """Write uncovered runs of at least length_cutoff bases, longest first."""
    rows = [(ref, s, e, e - s) for ref, intervals in uncovered_regions.items()
            for (s, e) in intervals if e - s >= length_cutoff]
    rows.sort(key=lambda x: x[3], reverse=True)
    with open(output_path, "w") as f:
        f.write("Reference\tStart\tEnd\tLength\n")
        for ref, s, e, length in rows:
            f.write(f"{ref}\t{s}\t{e}\t{length}\n")
    return len(rows)


class SequenceProcessor:
    """Processes sequences for FASTA export of uncovered regions."""
    
    @staticmethod
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

    @staticmethod
    def export_uncovered_fasta(
        uncovered_regions: Dict[str, List[Tuple[int, int]]],
        fasta_path: Path,
        output_path: Path,
        n_split_path: Path,
        extend_bp: int = 0,
        length_cutoff: int = 0,
        oligo_length: int = 120
    ) -> None:
        """
        Export uncovered regions as FASTA (optionally extended on both sides),
        and split regions at 'N' into smaller segments.
        """
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
                
                split_parts = SequenceProcessor.split_sequence_at_n(seq)
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


class SequenceLoader:
    """Loads reference sequences from FASTA files."""
    
    @staticmethod
    def load_reference_sequences(fasta_path: Path) -> Dict[str, str]:
        """Load reference sequences from FASTA file."""
        sequences = {}
        try:
            for record in SeqIO.parse(fasta_path, "fasta"):
                sequences[record.id] = str(record.seq)
            logging.info(f"Loaded {len(sequences)} reference sequences")
        except Exception as e:
            logging.error(f"Error loading reference sequences: {e}")
            raise
        return sequences


class ForcedOligoHandler:
    """Handles forced oligo lists."""
    
    @staticmethod
    def read_forced_oligos(path: Optional[Path]) -> Set[str]:
        """Read forced oligo IDs (one per line, case-sensitive)."""
        if not path:
            return set()
        with open(path) as f:
            return {line.strip() for line in f if line.strip()}

    @staticmethod
    def filter_bed_for_forced_oligos(
        bed: Any,
        forced_oligos: Set[str],
        temp_dir: Optional[Path] = None
    ) -> Any:
        """Restrict a BedTool to rows whose name is in forced_oligos."""
        if not forced_oligos:
            return bed
        entries = [str(interval) for interval in bed if interval.name in forced_oligos]
        if not entries:
            raise ValueError("No forced oligos matched the mapped oligos")
        tmp_file = "forced_temp.bed" if temp_dir is None else str(temp_dir / "forced_temp.bed")
        with open(tmp_file, "w") as f:
            f.write("".join(entries))
        return BedTool(tmp_file).sort()


class CoverageAnalysisOrchestrator:
    """Orchestrates coverage analysis workflow."""
    
    def __init__(self):
        """Initialize coverage analysis orchestrator."""
        self.psl_parser = PSLParser()
        self.genome_handler = GenomeFileHandler()
        self.mapping_builder = MappingBuilder()
        self.coverage_calculator = CoverageCalculator()
        self.sequence_loader = SequenceLoader()
        self.forced_handler = ForcedOligoHandler()
    
    def setup_analysis(
        self,
        psl_path: Path,
        min_length: int,
        min_similarity: float,
        temp_dir: Optional[Path] = None,
        reference_fasta: Optional[Path] = None
    ) -> Tuple[Any, Dict[str, List[OligoMapping]], str]:
        """Set up the analysis by parsing PSL and creating necessary files."""
        # Parse PSL to BED
        bed = self.psl_parser.parse_psl_to_bed(
            psl_path, min_length, min_similarity, temp_dir
        )
        
        # Build mappings
        mappings_dict = self.mapping_builder.build_mappings(bed)
        
        # Create genome file
        genome_file = self.genome_handler.create_genome_file(bed, temp_dir, reference_fasta)
        
        return bed, mappings_dict, genome_file
    
    def create_coverage_calculator(self, genome_file: str, temp_dir: Optional[Path] = None):
        """Create a coverage calculator function for use with multi-pass selection."""
        def calculate_coverage_for_selection(
            selected: Set,
            all_mappings: Dict[str, List[OligoMapping]],
            min_coverage: float
        ) -> Tuple[Dict[str, List[Tuple[int, int]]], int]:
            """Calculate coverage for the multi-pass selection algorithm."""
            coverage_bed = self.coverage_calculator.coverage_from_selected(
                selected, all_mappings, genome_file, temp_dir
            )
            uncovered_regions, total_uncovered, _ = self.coverage_calculator.compute_uncovered_regions(
                coverage_bed, min_coverage, genome_file
            )
            return uncovered_regions, total_uncovered
        
        return calculate_coverage_for_selection


class CoverageChecker:
    """Coverage check for the check command."""
    
    def __init__(self):
        self.psl_converter = PSLParser()
        self.genome_builder = GenomeFileHandler()
        self.coverage_calculator = CoverageCalculator()
        self.sequence_processor = SequenceProcessor()
        self.forced_filter = ForcedOligoHandler()
    
    def check_coverage(
        self,
        bed: Any,
        forced_oligos: Set[str],
        min_coverage: float,
        coverage_out: Optional[Path],
        longest_uncovered_out: Optional[Path],
        temp_dir: Optional[Path],
        uncovered_length_cutoff: int,
        args
    ) -> None:
        """
        Evaluate coverage of forced oligos, or of the whole BED when none are
        forced, then write optional coverage and uncovered-region outputs.
        """
        bed_filtered = self.forced_filter.filter_bed_for_forced_oligos(bed, forced_oligos, temp_dir)
        genome_file = self.genome_builder.create_genome_file(bed_filtered, temp_dir, args.reference)
        uncovered_regions, total_uncovered, coverage = \
            self.coverage_calculator.compute_uncovered_regions(bed_filtered, min_coverage, genome_file)
        
        total_bases = sum(int(line.strip().split()[1]) for line in open(genome_file))
        if total_bases > 0:
            logging.info(f"Total bases: {total_bases:,}")
            logging.info(f"Bases under min coverage: {total_uncovered:,} "
                        f"({total_uncovered/total_bases*100:.1f}%)")
        
        if coverage_out:
            with open(coverage_out, "w") as f:
                f.write("Reference\tStart\tEnd\tCoverage\n")
                for interval in coverage:
                    if interval.chrom != "genome":
                        f.write(f"{interval.chrom}\t{interval.start}\t{interval.end}\t{interval.name}\n")
            logging.info(f"Coverage data written to {coverage_out}")
        
        if longest_uncovered_out:
            count = write_uncovered_regions(uncovered_regions, longest_uncovered_out, uncovered_length_cutoff)
            logging.info(f"Wrote {count} uncovered stretches >= "
                        f"{uncovered_length_cutoff}bp to {longest_uncovered_out}")
        
        if args.uncovered_fasta and args.reference:
            if not args.n_split_fasta:
                logging.warning("--n_split_fasta not specified, skipping N-split output")
            else:
                self.sequence_processor.export_uncovered_fasta(
                    uncovered_regions, args.reference, args.uncovered_fasta, args.n_split_fasta,
                    args.extend_region, uncovered_length_cutoff, args.min_oligo_length)
