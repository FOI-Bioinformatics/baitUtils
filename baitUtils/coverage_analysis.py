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
    import pybedtools
    from pybedtools import BedTool
    HAS_PYBEDTOOLS = True
except ImportError:
    HAS_PYBEDTOOLS = False
    BedTool = None
    logging.warning("pybedtools not available - coverage analysis may be limited")

from tqdm import tqdm
from Bio import SeqIO

from baitUtils.gap_filling_algorithm import OligoMapping


class PSLParser:
    """Parses PSL files into BED format."""
    
    @staticmethod
    def parse_psl_to_bed(
        psl_path: Path,
        min_length: int,
        min_similarity: float,
        temp_dir: Optional[Path] = None
    ) -> Any:
        """
        Parse PSL file into BED format (keeping only hits above min_length and min_similarity).
        """
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


class GenomeFileHandler:
    """Handles genome file creation and management."""
    
    @staticmethod
    def create_genome_file(bed: Any, temp_dir: Optional[Path] = None) -> str:
        """
        Build a genome file from the sorted BED intervals 
        (required by bedtools genomecov).
        """
        chrom_sizes = defaultdict(int)
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


class CoverageCalculator:
    """Calculates coverage statistics and identifies uncovered regions."""
    
    @staticmethod
    def coverage_from_selected(
        selected_oligos: Set[str],
        all_mappings: Dict[str, List[OligoMapping]],
        genome_file: str,
        temp_dir: Optional[Path] = None
    ) -> Any:
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

    @staticmethod
    def calculate_coverage(
        bed: Any,
        min_coverage: float,
        max_coverage: Optional[float],
        temp_dir: Optional[Path],
        genome_file: str
    ) -> Tuple[
        Dict[str, List[Tuple[int, int]]], 
        int, 
        int, 
        int, 
        List[Tuple[str, int, int, float]], 
        Dict[str, int]
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
        """Read forced oligos from file."""
        if not path:
            return set()
        try:
            with open(path) as f:
                return {line.strip().upper() for line in f if line.strip()}
        except Exception as e:
            logging.error(f"Error reading forced oligos: {e}")
            raise


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
        temp_dir: Optional[Path] = None
    ) -> Tuple[Any, Dict[str, List[OligoMapping]], str]:
        """Set up the analysis by parsing PSL and creating necessary files."""
        # Parse PSL to BED
        bed = self.psl_parser.parse_psl_to_bed(
            psl_path, min_length, min_similarity, temp_dir
        )
        
        # Build mappings
        mappings_dict = self.mapping_builder.build_mappings(bed)
        
        # Create genome file
        genome_file = self.genome_handler.create_genome_file(bed, temp_dir)
        
        return bed, mappings_dict, genome_file
    
    def create_coverage_calculator(self, genome_file: str, temp_dir: Optional[Path] = None):
        """Create a coverage calculator function for use with multi-pass selection."""
        def calculate_coverage_for_selection(
            selected_oligos: Set[str],
            all_mappings: Dict[str, List[OligoMapping]],
            min_coverage: float,
            max_coverage: Optional[float]
        ) -> Tuple[Dict[str, List[Tuple[int, int]]], int]:
            """Calculate coverage for the multi-pass selection algorithm."""
            coverage_bed = self.coverage_calculator.coverage_from_selected(
                selected_oligos, all_mappings, genome_file, temp_dir
            )
            uncovered_regions, total_uncovered, _ = self.coverage_calculator.compute_uncovered_regions(
                coverage_bed, min_coverage, max_coverage, genome_file
            )
            return uncovered_regions, total_uncovered
        
        return calculate_coverage_for_selection