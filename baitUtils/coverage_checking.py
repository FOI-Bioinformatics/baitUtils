#!/usr/bin/env python3

"""
coverage_checking.py

Core coverage checking utilities for evaluating oligo coverage and 
identifying uncovered regions. Refactored from check.py for better organization.
"""

import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional, Any
from collections import defaultdict

try:
    import pybedtools
    from pybedtools import BedTool
    HAS_PYBEDTOOLS = True
except ImportError:
    HAS_PYBEDTOOLS = False
    BedTool = None

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


def _merge_uncovered_intervals(intervals, min_coverage):
    from baitUtils.coverage_analysis import merge_uncovered_intervals
    return merge_uncovered_intervals(intervals, min_coverage)


class PSLToBedConverter:
    """Converts PSL files to BED format for coverage analysis (see coverage_analysis.PSLParser)."""
    
    @staticmethod
    def parse_psl_to_bed(
        psl_path: Path,
        min_length: int,
        min_similarity: float,
        temp_dir: Optional[Path] = None
    ) -> Any:
        from baitUtils.coverage_analysis import PSLParser
        return PSLParser.parse_psl_to_bed(psl_path, min_length, min_similarity, temp_dir, "temp_check.bed")


class GenomeFileBuilder:
    """Builds genome files for bedtools operations."""
    
    @staticmethod
    def create_genome_file(bed: Any, temp_dir: Optional[Path] = None,
                           reference_fasta: Optional[Path] = None) -> str:
        """Build a genome file for bedtools genomecov (see GenomeFileHandler)."""
        from baitUtils.coverage_analysis import GenomeFileHandler
        return GenomeFileHandler.create_genome_file(bed, temp_dir, reference_fasta)


class UncoveredRegionAnalyzer:
    """Analyzes and identifies uncovered regions from coverage data."""
    
    @staticmethod
    def compute_uncovered_regions(
        coverage_bed: Any,
        min_coverage: float,
        max_coverage: Optional[float],
        genome_file: str
    ) -> Tuple[Dict[str, List[Tuple[int, int]]], int, Any]:
        """
        Compute uncovered regions below min_coverage using bedtools genomecov output.
        Returns dictionary of uncovered intervals, total uncovered bases, and the coverage tool.
        """
        coverage = coverage_bed.genomecov(bga=True, g=genome_file)
        uncovered_regions, total_uncovered = _merge_uncovered_intervals(coverage, min_coverage)
        return uncovered_regions, total_uncovered, coverage

    @staticmethod
    def write_uncovered_regions(
        uncovered_regions: List[Tuple[str, int, int, float]], 
        output_path: Path,
        length_cutoff: int = 0
    ) -> int:
        """Write uncovered regions above a certain length cutoff to file."""
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


class ForcedOligoFilter:
    """Handles filtering for forced oligos."""
    
    @staticmethod
    def read_forced_oligos(path: Optional[Path]) -> Set[str]:
        """Read forced oligos from file."""
        if not path:
            return set()
        try:
            with open(path) as f:
                return {line.strip() for line in f if line.strip()}
        except Exception as e:
            logging.error(f"Error reading forced oligos: {e}")
            sys.exit(1)

    @staticmethod
    def filter_bed_for_forced_oligos(
        bed: Any,
        forced_oligos: Set[str],
        temp_dir: Optional[Path] = None
    ) -> Any:
        """Filter BED to include only forced oligos."""
        if not forced_oligos:
            return bed
        
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
        
        return BedTool(tmp_file).sort()


class CoverageChecker:
    """Main coverage checking orchestrator."""
    
    def __init__(self):
        """Initialize coverage checker."""
        self.psl_converter = PSLToBedConverter()
        self.genome_builder = GenomeFileBuilder()
        self.uncovered_analyzer = UncoveredRegionAnalyzer()
        self.sequence_processor = SequenceProcessor()
        self.forced_filter = ForcedOligoFilter()
    
    def check_coverage(
        self,
        bed: Any,
        forced_oligos: Set[str],
        min_coverage: float,
        max_coverage: Optional[float],
        coverage_out: Optional[Path],
        longest_uncovered_out: Optional[Path],
        temp_dir: Optional[Path],
        uncovered_length_cutoff: int,
        args
    ) -> None:
        """
        Evaluate coverage of forced oligos or entire BED if none forced,
        then write optional coverage data and uncovered intervals.
        """
        # Filter BED to forced oligos if specified
        bed_filtered = self.forced_filter.filter_bed_for_forced_oligos(
            bed, forced_oligos, temp_dir
        )
        
        # Create genome file & compute coverage
        genome_file = self.genome_builder.create_genome_file(bed_filtered, temp_dir, args.reference)
        uncovered_regions, total_uncovered, coverage = \
            self.uncovered_analyzer.compute_uncovered_regions(
                bed_filtered, min_coverage, max_coverage, genome_file
            )
        
        # Summarize coverage stats
        total_bases = sum(
            int(line.strip().split()[1])
            for line in open(genome_file)
        )
        
        if total_bases > 0:
            logging.info(f"Total bases: {total_bases:,}")
            logging.info(f"Bases under min coverage: {total_uncovered:,} "
                        f"({total_uncovered/total_bases*100:.1f}%)")
        
        # Write optional outputs
        self._write_coverage_outputs(
            coverage, coverage_out, uncovered_regions, 
            longest_uncovered_out, uncovered_length_cutoff, args
        )
    
    def _write_coverage_outputs(
        self,
        coverage,
        coverage_out: Optional[Path],
        uncovered_regions: Dict[str, List[Tuple[int, int]]],
        longest_uncovered_out: Optional[Path],
        uncovered_length_cutoff: int,
        args
    ) -> None:
        """Write coverage and uncovered region outputs."""
        # Coverage data
        if coverage_out:
            with open(coverage_out, "w") as f:
                f.write("Reference\tStart\tEnd\tCoverage\n")
                for interval in coverage:
                    if interval.chrom != "genome":
                        f.write(f"{interval.chrom}\t{interval.start}\t{interval.end}\t{interval.name}\n")
            logging.info(f"Coverage data written to {coverage_out}")
        
        # Uncovered intervals
        if longest_uncovered_out:
            zero_cov_list = []
            for chrom, intervals in uncovered_regions.items():
                for (s, e) in intervals:
                    zero_cov_list.append((chrom, s, e, 0.0))
            
            count = self.uncovered_analyzer.write_uncovered_regions(
                zero_cov_list, longest_uncovered_out, uncovered_length_cutoff
            )
            logging.info(f"Wrote {count} uncovered stretches >= "
                        f"{uncovered_length_cutoff}bp to {longest_uncovered_out}")
        
        # Uncovered FASTA
        if args.uncovered_fasta and args.reference:
            if not args.n_split_fasta:
                logging.warning("--n_split_fasta not specified, skipping N-split output")
            else:
                self.sequence_processor.export_uncovered_fasta(
                    uncovered_regions,
                    args.reference,
                    args.uncovered_fasta,
                    args.n_split_fasta,
                    args.extend_region,
                    uncovered_length_cutoff,
                    args.min_oligo_length
                )