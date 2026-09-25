#!/usr/bin/env python3

"""
mapping_utils.py

Core mapping utilities for sequence mapping operations using external tools.
Provides functionality for running mappers, parsing results, and managing outputs.
"""

import logging
import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Set, Dict, Optional, List, Iterator, Iterable, Tuple, Union
from pathlib import Path

import pandas as pd
from Bio import SeqIO


PSL_HEADER_PREFIXES = ("psLayout", "match", "-", "#", "no matches")


@dataclass
class PSLHit:
    """One alignment from a BLAT/pblat PSL file (all 21 columns)."""
    matches: int
    mismatches: int
    rep_matches: int
    n_count: int
    q_num_insert: int
    q_base_insert: int
    t_num_insert: int
    t_base_insert: int
    strand: str
    q_name: str
    q_size: int
    q_start: int
    q_end: int
    t_name: str
    t_size: int
    t_start: int
    t_end: int
    block_sizes: List[int] = field(default_factory=list)
    q_starts: List[int] = field(default_factory=list)
    t_starts: List[int] = field(default_factory=list)
    line: str = ""

    @property
    def aligned_length(self) -> int:
        """Number of query bases inside aligned blocks."""
        return sum(self.block_sizes) if self.block_sizes else (self.q_end - self.q_start)

    @property
    def target_span(self) -> int:
        """Target bases from first to last aligned base, including target inserts."""
        return self.t_end - self.t_start

    @property
    def target_blocks(self) -> List[Tuple[int, int]]:
        """Half-open target intervals of the aligned blocks."""
        if not self.block_sizes:
            return [(self.t_start, self.t_end)]
        return [(ts, ts + size) for ts, size in zip(self.t_starts, self.block_sizes)]

    def milli_bad(self, is_mrna: bool = True) -> float:
        """
        BLAT's pslCalcMilliBad: mismatches per thousand aligned bases,
        penalising query inserts and (unless is_mrna) target inserts and
        alignment size differences. pblat applies -minIdentity with
        is_mrna=True, so that is the default here.
        """
        q_ali = self.q_end - self.q_start
        t_ali = self.t_end - self.t_start
        ali_size = min(q_ali, t_ali)
        if ali_size <= 0:
            return 0.0
        size_dif = q_ali - t_ali
        if size_dif < 0:
            size_dif = 0 if is_mrna else -size_dif
        insert_factor = self.q_num_insert
        if not is_mrna:
            insert_factor += self.t_num_insert
        total = self.matches + self.rep_matches + self.mismatches
        if total == 0:
            return 0.0
        penalty = self.mismatches + insert_factor + round(3 * math.log(1 + size_dif))
        return 1000.0 * penalty / total

    @property
    def identity(self) -> float:
        """Percent identity as BLAT reports it (100 - milliBad / 10)."""
        return 100.0 - self.milli_bad() / 10.0


def _int_list(text: str) -> List[int]:
    return [int(x) for x in text.strip().split(",") if x]


def parse_psl_line(line: str) -> Optional[PSLHit]:
    """Parse one PSL data line; return None when it is not a valid 21-column row."""
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 21:
        cols = line.split()
    if len(cols) < 21:
        return None
    try:
        return PSLHit(
            matches=int(cols[0]), mismatches=int(cols[1]), rep_matches=int(cols[2]),
            n_count=int(cols[3]), q_num_insert=int(cols[4]), q_base_insert=int(cols[5]),
            t_num_insert=int(cols[6]), t_base_insert=int(cols[7]), strand=cols[8],
            q_name=cols[9], q_size=int(cols[10]), q_start=int(cols[11]), q_end=int(cols[12]),
            t_name=cols[13], t_size=int(cols[14]), t_start=int(cols[15]), t_end=int(cols[16]),
            block_sizes=_int_list(cols[18]), q_starts=_int_list(cols[19]),
            t_starts=_int_list(cols[20]), line=line.rstrip("\n"),
        )
    except (ValueError, IndexError):
        return None


def is_psl_header(line: str) -> bool:
    return line.startswith(PSL_HEADER_PREFIXES)


def parse_psl(psl_path: Union[str, Path]) -> Iterator[PSLHit]:
    """
    Yield PSLHit records from a PSL file, skipping header lines.

    Malformed data lines are skipped and counted; a single warning with the
    count is logged at the end so that problems are visible without flooding
    the log.
    """
    skipped = 0
    with open(psl_path) as fh:
        for line in fh:
            if not line.strip() or is_psl_header(line):
                continue
            hit = parse_psl_line(line)
            if hit is None:
                skipped += 1
                continue
            yield hit
    if skipped:
        logging.warning(f"Skipped {skipped} malformed line(s) in {psl_path}")


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_paf_line(line: str) -> Optional[PSLHit]:
    """
    Convert one PAF line (minimap2) into a PSLHit.

    Aligned blocks come from the cg:Z cigar when present (minimap2 -c);
    without it the hit is one block spanning the target interval and the
    mismatch count is taken from the alignment block length minus matches.
    Query block starts are given in forward query coordinates.
    """
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 12:
        return None
    try:
        q_name, q_size, q_start, q_end, strand = cols[0], int(cols[1]), int(cols[2]), int(cols[3]), cols[4]
        t_name, t_size, t_start, t_end = cols[5], int(cols[6]), int(cols[7]), int(cols[8])
        matches, block_len = int(cols[9]), int(cols[10])
    except ValueError:
        return None
    if strand not in ("+", "-"):
        return None

    tags = {}
    for tag in cols[12:]:
        parts = tag.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]

    block_sizes: List[int] = []
    q_starts: List[int] = []
    t_starts: List[int] = []
    q_num_insert = q_base_insert = t_num_insert = t_base_insert = 0
    cigar = tags.get("cg")
    if cigar:
        q_pos, t_pos = q_start, t_start
        for length, op in _CIGAR_RE.findall(cigar):
            length = int(length)
            if op in ("M", "=", "X"):
                if block_sizes and q_starts[-1] + block_sizes[-1] == q_pos and t_starts[-1] + block_sizes[-1] == t_pos:
                    block_sizes[-1] += length
                else:
                    block_sizes.append(length)
                    q_starts.append(q_pos)
                    t_starts.append(t_pos)
                q_pos += length
                t_pos += length
            elif op == "I":
                q_num_insert += 1
                q_base_insert += length
                q_pos += length
            elif op in ("D", "N"):
                t_num_insert += 1
                t_base_insert += length
                t_pos += length
            # S, H and P do not occur in PAF cigars from minimap2
    if not block_sizes:
        block_sizes, q_starts, t_starts = [min(q_end - q_start, t_end - t_start)], [q_start], [t_start]

    aligned = sum(block_sizes)
    if "NM" in tags:
        mismatches = max(0, int(tags["NM"]) - q_base_insert - t_base_insert)
    else:
        mismatches = max(0, block_len - matches - q_base_insert - t_base_insert)
    matches = min(matches, aligned)

    return PSLHit(
        matches=matches, mismatches=mismatches, rep_matches=0, n_count=0,
        q_num_insert=q_num_insert, q_base_insert=q_base_insert,
        t_num_insert=t_num_insert, t_base_insert=t_base_insert,
        strand=strand, q_name=q_name, q_size=q_size, q_start=q_start, q_end=q_end,
        t_name=t_name, t_size=t_size, t_start=t_start, t_end=t_end,
        block_sizes=block_sizes, q_starts=q_starts, t_starts=t_starts, line=line.rstrip("\n"),
    )


def parse_paf(paf_path: Union[str, Path]) -> Iterator[PSLHit]:
    """Yield PSLHit records from a PAF file, skipping comment and malformed lines."""
    skipped = 0
    with open(paf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            hit = parse_paf_line(line)
            if hit is None:
                skipped += 1
                continue
            yield hit
    if skipped:
        logging.warning(f"Skipped {skipped} malformed line(s) in {paf_path}")


def alignment_format(path: Union[str, Path]) -> str:
    """'paf' for .paf files (optionally gzipped name), otherwise 'psl'."""
    name = str(path).lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return "paf" if name.endswith(".paf") else "psl"


def parse_alignments(path: Union[str, Path]) -> Iterator[PSLHit]:
    """Yield PSLHit records from a PSL or PAF file, chosen by file extension."""
    if alignment_format(path) == "paf":
        return parse_paf(path)
    return parse_psl(path)


def run_minimap2(reference: str, query: str, output_paf: str, threads: int = 1,
                 preset: str = "sr", max_secondary: int = 10) -> None:
    """
    Run minimap2 with base-level alignment (-c) so that cigars and NM tags
    are available, keeping secondary hits for off-target assessment.
    """
    cmd = ["minimap2", "-c", "-x", preset, "--secondary=yes", f"-N{max_secondary}",
           "-t", str(threads), str(reference), str(query)]
    with open(output_paf, "w") as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"minimap2 failed with return code {result.returncode}: {result.stderr.strip()}")
    if result.stderr:
        logging.debug(f"minimap2 stderr: {result.stderr.strip()}")


MAPPERS = ("pblat", "minimap2")


def check_mapper_available(mapper: str) -> bool:
    """True when the mapper binary is on PATH."""
    return shutil.which(mapper) is not None


def run_mapper(mapper: str, reference: str, query: str, output_dir: str, basename: str,
               threads: int = 1, min_identity: float = 90.0, minimap2_preset: str = "sr",
               min_score: int = 30, min_match: int = 2) -> str:
    """Run the chosen mapper and return the alignment file path (PSL for pblat, PAF for minimap2)."""
    if mapper == "pblat":
        output = os.path.join(output_dir, f"{basename}.psl")
        run_pblat(reference, query, output, threads=threads, min_identity=min_identity,
                  min_score=min_score, min_match=min_match)
    elif mapper == "minimap2":
        output = os.path.join(output_dir, f"{basename}.paf")
        run_minimap2(reference, query, output, threads=threads, preset=minimap2_preset)
    else:
        raise ValueError(f"Unsupported mapper: {mapper}")
    return output


def build_hit_table(hits: Iterable[PSLHit]) -> pd.DataFrame:
    """
    Summarise hits per query (bait) for off-target assessment.

    Columns: oligo_id, n_hits, n_targets, best_identity, second_best_identity,
    best_target, best_start, best_end, best_strand, best_aligned_length.
    Hits are ranked by identity, then aligned length.
    """
    per_query: Dict[str, List[PSLHit]] = {}
    for hit in hits:
        per_query.setdefault(hit.q_name, []).append(hit)

    rows = []
    for q_name, q_hits in per_query.items():
        ranked = sorted(q_hits, key=lambda h: (h.identity, h.aligned_length), reverse=True)
        best = ranked[0]
        rows.append({
            'oligo_id': q_name,
            'n_hits': len(ranked),
            'n_targets': len({h.t_name for h in ranked}),
            'best_identity': round(best.identity, 2),
            'second_best_identity': round(ranked[1].identity, 2) if len(ranked) > 1 else float('nan'),
            'best_target': best.t_name,
            'best_start': best.t_start,
            'best_end': best.t_end,
            'best_strand': best.strand,
            'best_aligned_length': best.aligned_length,
        })
    columns = ['oligo_id', 'n_hits', 'n_targets', 'best_identity', 'second_best_identity',
               'best_target', 'best_start', 'best_end', 'best_strand', 'best_aligned_length']
    return pd.DataFrame(rows, columns=columns).sort_values('oligo_id').reset_index(drop=True)


STRAND_CHOICES = ("both", "plus", "minus")
_STRAND_SYMBOL = {"plus": "+", "minus": "-"}


def filter_hits(hits: Iterable[PSLHit], min_identity: float = 0.0, min_length: int = 0,
                min_matches: int = 0, strand: str = "both") -> Iterator[PSLHit]:
    """
    Keep hits with identity, aligned length and match count at or above the
    thresholds; strand restricts to hits on the given target strand.
    """
    wanted = _STRAND_SYMBOL.get(strand)
    for hit in hits:
        if wanted is not None and hit.strand != wanted:
            continue
        if hit.aligned_length < min_length or hit.matches < min_matches:
            continue
        if hit.identity < min_identity:
            continue
        yield hit


def run_pblat(reference: str, query: str, output_psl: str, threads: int = 1,
              min_identity: float = 90.0, min_score: int = 30, min_match: int = 2) -> None:
    """
    Run pblat and raise RuntimeError on failure.

    When pblat is killed by a signal (negative return code) and more than one
    thread was requested, the run is retried once with a single thread, since
    multi-threaded pblat has been observed to die sporadically on small inputs.
    """
    def command(n_threads: int) -> List[str]:
        return ['pblat', f'-threads={n_threads}', f'-minIdentity={min_identity}',
                f'-minScore={min_score}', f'-minMatch={min_match}',
                str(reference), str(query), str(output_psl)]

    result = subprocess.run(command(threads), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode < 0 and threads > 1:
        logging.warning(f"pblat was terminated by signal {-result.returncode} with {threads} threads; "
                        "retrying with a single thread")
        result = subprocess.run(command(1), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"pblat failed with return code {result.returncode}: {result.stderr.strip()}")
    if result.stderr:
        logging.debug(f"pblat stderr: {result.stderr.strip()}")


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
    """Runs pblat for the map command (see run_pblat for the retry behaviour)."""
    
    @staticmethod
    def run_pblat(
        input_file: str,
        target_file: str,
        output_file: str,
        threads: int = 1,
        min_match: int = 2,
        min_score: int = 30,
        min_identity: int = 90
    ) -> None:
        """Map input_file (baits) against target_file (reference) into output_file."""
        logging.info(f"Running pblat with {threads} thread(s)...")
        run_pblat(target_file, input_file, output_file, threads=threads,
                  min_identity=min_identity, min_score=min_score, min_match=min_match)
        logging.info(f"pblat mapping completed: {output_file}")


class PSLParser:
    last_hit_table = None  # hit table from the most recent parse_psl_file call
    """Parses PSL files and extracts mapping information."""
    
    @staticmethod
    def parse_psl_file(
        psl_file: str,
        min_identity: float = 90.0,
        min_match_count: int = 0,
        filtered_output: Optional[str] = None,
        min_length: int = 0
    ) -> Set[str]:
        """
        Return the set of query IDs with at least one hit passing the filters.

        Args:
            psl_file: Path to PSL file
            min_identity: Minimum BLAT identity percentage
            min_match_count: Minimum number of matching bases
            filtered_output: Optional path for the passing PSL rows
            min_length: Minimum aligned length in bases
        """
        mapped_sequences = set()
        filtered_hits = []
        kept = []
        for hit in filter_hits(parse_alignments(psl_file), min_identity, min_length, min_match_count):
            mapped_sequences.add(hit.q_name)
            kept.append(hit)
            if filtered_output:
                filtered_hits.append(hit.line)
        PSLParser.last_hit_table = build_hit_table(kept)

        if filtered_output:
            header_lines = []
            if alignment_format(psl_file) == "psl":
                with open(psl_file) as fh:
                    for line in fh:
                        if is_psl_header(line):
                            header_lines.append(line.rstrip("\n"))
                        elif line.strip():
                            break
            PSLParser._write_filtered_psl(filtered_output, header_lines, filtered_hits)
            logging.info(f"Wrote {len(filtered_hits)} filtered hits to {filtered_output}")

        return mapped_sequences

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
        output_file: str,
        reverse_complement: Optional[Set[str]] = None
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
            reverse_complement = reverse_complement or set()
            with open(output_file, 'w') as f:
                for seq_id in sorted(sequence_ids):
                    if seq_id in seq_records:
                        record = seq_records[seq_id]
                        if seq_id in reverse_complement:
                            record = record.reverse_complement(id=record.id, name=record.name,
                                                               description=f"{record.description} reverse_complemented")
                        SeqIO.write(record, f, 'fasta')
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
        min_identity: int = 90,
        minimap2_preset: str = 'sr'
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
        return run_mapper(mapper, target_file, input_file, output_dir, f"{output_prefix}-mapping",
                          threads=threads, min_identity=min_identity, minimap2_preset=minimap2_preset,
                          min_score=min_score, min_match=min_match)


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