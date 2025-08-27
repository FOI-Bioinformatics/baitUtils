#!/usr/bin/env python3

"""
coverage_stats.py

Comprehensive coverage statistics analysis for oligo set evaluation.
Analyzes PSL mapping files and reference sequences to compute detailed
coverage metrics including breadth, depth, uniformity, and distributions.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from datetime import datetime

from Bio import SeqIO

try:
    import pybedtools
    from pybedtools import BedTool
    HAS_PYBEDTOOLS = True
except ImportError:
    HAS_PYBEDTOOLS = False
    logging.warning("pybedtools not available - some coverage analysis features may be limited")


class CoverageAnalyzer:
    """Comprehensive coverage statistics analyzer."""
    
    def __init__(
        self,
        psl_file: Path,
        reference_file: Path,
        min_coverage: float = 1.0,
        target_coverage: float = 10.0,
        min_identity: float = 90.0,
        min_length: int = 100
    ):
        """
        Initialize the coverage analyzer.
        
        Args:
            psl_file: Path to PSL mapping file
            reference_file: Path to reference FASTA file
            min_coverage: Minimum coverage threshold
            target_coverage: Target coverage depth for analysis
            min_identity: Minimum mapping identity to consider
            min_length: Minimum mapping length to consider
        """
        self.psl_file = Path(psl_file)
        self.reference_file = Path(reference_file)
        self.min_coverage = min_coverage
        self.target_coverage = target_coverage
        self.min_identity = min_identity
        self.min_length = min_length
        
        # Data storage
        self.reference_sequences = {}
        self.mappings = []
        self.coverage_arrays = {}
        self.stats = {}
    
    def analyze(self) -> Dict[str, Any]:
        """
        Perform comprehensive coverage analysis.
        
        Returns:
            Dictionary containing all coverage statistics
        """
        logging.info("Loading reference sequences...")
        self._load_reference_sequences()
        
        logging.info("Parsing mapping data...")
        self._parse_psl_file()
        
        logging.info("Computing coverage arrays...")
        self._compute_coverage_arrays()
        
        logging.info("Calculating statistics...")
        self._calculate_statistics()
        
        return self.stats
    
    def _load_reference_sequences(self) -> None:
        """Load reference sequences from FASTA file."""
        try:
            with open(self.reference_file, 'r') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    self.reference_sequences[record.id] = {
                        'sequence': str(record.seq),
                        'length': len(record.seq)
                    }
            
            total_length = sum(ref['length'] for ref in self.reference_sequences.values())
            logging.info(f"Loaded {len(self.reference_sequences)} reference sequences "
                        f"({total_length:,} bp total)")
            
        except Exception as e:
            logging.error(f"Error loading reference sequences: {e}")
            raise
    
    def _parse_psl_file(self) -> None:
        """Parse PSL file and extract valid mappings."""
        valid_mappings = 0
        total_lines = 0
        
        try:
            with open(self.psl_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    
                    # Skip header and empty lines
                    if (line.startswith(('psLayout', 'match', '-', '#')) or 
                        not line or line.startswith('no matches')):
                        continue
                    
                    total_lines += 1
                    parts = line.split()
                    
                    if len(parts) < 21:
                        continue
                    
                    try:
                        # Parse PSL fields
                        matches = int(parts[0])
                        mismatches = int(parts[1])
                        rep_matches = int(parts[2])
                        n_count = int(parts[3])
                        q_name = parts[9]
                        q_size = int(parts[10])
                        t_name = parts[13]
                        t_size = int(parts[14])
                        t_start = int(parts[15])
                        t_end = int(parts[16])
                        
                        # Apply filters
                        mapping_length = t_end - t_start
                        if mapping_length < self.min_length:
                            continue
                        
                        # Calculate identity
                        total_aligned = matches + mismatches + rep_matches
                        if total_aligned == 0:
                            continue
                        
                        identity = (matches + rep_matches) / total_aligned * 100.0
                        if identity < self.min_identity:
                            continue
                        
                        # Store valid mapping
                        mapping = {
                            'query_name': q_name,
                            'target_name': t_name,
                            'target_start': t_start,
                            'target_end': t_end,
                            'length': mapping_length,
                            'identity': identity,
                            'matches': matches
                        }
                        
                        self.mappings.append(mapping)
                        valid_mappings += 1
                        
                    except (ValueError, IndexError) as e:
                        logging.debug(f"Error parsing PSL line: {e}")
                        continue
            
            logging.info(f"Parsed {valid_mappings} valid mappings from {total_lines} total lines")
            
        except Exception as e:
            logging.error(f"Error parsing PSL file: {e}")
            raise
    
    def _compute_coverage_arrays(self) -> None:
        """Compute coverage depth arrays for each reference sequence."""
        # Initialize coverage arrays
        for ref_id, ref_data in self.reference_sequences.items():
            self.coverage_arrays[ref_id] = np.zeros(ref_data['length'], dtype=np.int32)
        
        # Add coverage from mappings
        for mapping in self.mappings:
            ref_id = mapping['target_name']
            start = mapping['target_start']
            end = mapping['target_end']
            
            if ref_id in self.coverage_arrays:
                # Ensure indices are within bounds
                start = max(0, start)
                end = min(len(self.coverage_arrays[ref_id]), end)
                
                if start < end:
                    self.coverage_arrays[ref_id][start:end] += 1
        
        logging.info("Coverage arrays computed for all reference sequences")
    
    def _calculate_statistics(self) -> None:
        """Calculate comprehensive coverage statistics."""
        self.stats = {
            'analysis_date': datetime.now().isoformat(),
            'parameters': {
                'min_coverage': self.min_coverage,
                'target_coverage': self.target_coverage,
                'min_identity': self.min_identity,
                'min_length': self.min_length
            }
        }
        
        # Basic mapping statistics
        unique_queries = len(set(m['query_name'] for m in self.mappings))
        self.stats.update({
            'total_mappings': len(self.mappings),
            'mapped_oligos': unique_queries,
            'total_oligos': self._count_total_oligos(),
            'mapping_efficiency': (unique_queries / max(1, self._count_total_oligos())) * 100
        })
        
        # Reference statistics
        total_ref_length = sum(ref['length'] for ref in self.reference_sequences.values())
        self.stats.update({
            'reference_sequences': len(self.reference_sequences),
            'reference_length': total_ref_length
        })
        
        # Coverage depth statistics
        self._calculate_depth_statistics()
        
        # Coverage breadth statistics
        self._calculate_breadth_statistics()
        
        # Coverage uniformity statistics
        self._calculate_uniformity_statistics()
        
        # Per-reference statistics
        self._calculate_per_reference_statistics()
        
        logging.info("Comprehensive statistics calculated")
    
    def _count_total_oligos(self) -> int:
        """Count total number of oligos from original input (estimate from mappings)."""
        # This is an approximation since we only see mapped oligos
        # In practice, this would come from counting the input FASTA
        return len(set(m['query_name'] for m in self.mappings))
    
    def _calculate_depth_statistics(self) -> None:
        """Calculate coverage depth statistics."""
        # Concatenate all coverage arrays
        all_coverage = np.concatenate(list(self.coverage_arrays.values()))
        
        # Basic depth statistics
        self.stats.update({
            'mean_depth': float(np.mean(all_coverage)),
            'median_depth': float(np.median(all_coverage)),
            'depth_std': float(np.std(all_coverage)),
            'max_depth': int(np.max(all_coverage)),
            'min_depth': int(np.min(all_coverage))
        })
        
        # Depth distribution
        depth_thresholds = [1, 5, 10, 20, 50, 100]
        depth_dist = {}
        
        for threshold in depth_thresholds:
            bases_above = np.sum(all_coverage >= threshold)
            percentage = (bases_above / len(all_coverage)) * 100.0
            depth_dist[threshold] = percentage
        
        self.stats['depth_distribution'] = depth_dist
        
        # Coverage histogram data
        hist, bin_edges = np.histogram(all_coverage, bins=50)
        self.stats['depth_histogram'] = {
            'counts': hist.tolist(),
            'bin_edges': bin_edges.tolist()
        }
    
    def _calculate_breadth_statistics(self) -> None:
        """Calculate coverage breadth statistics."""
        total_bases = sum(len(cov_array) for cov_array in self.coverage_arrays.values())
        covered_bases = sum(
            np.sum(cov_array >= self.min_coverage) 
            for cov_array in self.coverage_arrays.values()
        )
        
        self.stats.update({
            'total_bases': total_bases,
            'covered_bases': int(covered_bases),
            'uncovered_bases': total_bases - int(covered_bases),
            'coverage_breadth': (covered_bases / total_bases * 100.0) if total_bases > 0 else 0.0
        })
        
        # Target coverage breadth
        target_covered = sum(
            np.sum(cov_array >= self.target_coverage)
            for cov_array in self.coverage_arrays.values()
        )
        
        self.stats['target_coverage_breadth'] = (
            (target_covered / total_bases * 100.0) if total_bases > 0 else 0.0
        )
    
    def _calculate_uniformity_statistics(self) -> None:
        """Calculate coverage uniformity statistics."""
        all_coverage = np.concatenate(list(self.coverage_arrays.values()))
        covered_positions = all_coverage[all_coverage >= self.min_coverage]
        
        if len(covered_positions) > 0:
            # Coefficient of variation for covered positions
            cv = np.std(covered_positions) / np.mean(covered_positions)
            
            # Gini coefficient (measure of inequality)
            gini = self._calculate_gini_coefficient(covered_positions)
            
            self.stats.update({
                'coverage_cv': float(cv),
                'coverage_gini': float(gini),
                'uniformity_score': max(0.0, 1.0 - cv)  # Simple uniformity metric
            })
        else:
            self.stats.update({
                'coverage_cv': float('inf'),
                'coverage_gini': 1.0,
                'uniformity_score': 0.0
            })
    
    def _calculate_gini_coefficient(self, values: np.ndarray) -> float:
        """Calculate Gini coefficient for coverage uniformity."""
        if len(values) == 0:
            return 1.0
        
        sorted_values = np.sort(values)
        n = len(sorted_values)
        
        # Calculate Gini coefficient
        index = np.arange(1, n + 1)
        gini = (2 * np.sum(index * sorted_values)) / (n * np.sum(sorted_values)) - (n + 1) / n
        
        return gini
    
    def _calculate_per_reference_statistics(self) -> None:
        """Calculate statistics for each reference sequence."""
        per_ref_stats = {}
        
        for ref_id, cov_array in self.coverage_arrays.items():
            ref_length = len(cov_array)
            covered_bases = np.sum(cov_array >= self.min_coverage)
            
            stats = {
                'length': ref_length,
                'mean_depth': float(np.mean(cov_array)),
                'max_depth': int(np.max(cov_array)),
                'covered_bases': int(covered_bases),
                'coverage_breadth': (covered_bases / ref_length * 100.0) if ref_length > 0 else 0.0,
                'gaps': self._count_gaps_in_sequence(cov_array)
            }
            
            per_ref_stats[ref_id] = stats
        
        self.stats['per_reference'] = per_ref_stats
    
    def _count_gaps_in_sequence(self, coverage_array: np.ndarray) -> int:
        """Count number of coverage gaps in a sequence."""
        # Find positions below minimum coverage
        below_threshold = coverage_array < self.min_coverage
        
        # Count transitions from covered to uncovered
        if len(below_threshold) <= 1:
            return 0
        
        # Find gap boundaries
        gap_starts = np.where(np.diff(below_threshold.astype(int)) == 1)[0] + 1
        gap_ends = np.where(np.diff(below_threshold.astype(int)) == -1)[0] + 1
        
        # Handle edge cases
        if below_threshold[0]:
            gap_starts = np.concatenate([[0], gap_starts])
        if below_threshold[-1]:
            gap_ends = np.concatenate([gap_ends, [len(below_threshold)]])
        
        return len(gap_starts)
    
    def export_coverage_data(self, output_dir: Path) -> None:
        """Export detailed coverage data to files."""
        data_dir = output_dir / "data"
        data_dir.mkdir(exist_ok=True)
        
        # Export per-position coverage
        coverage_data = []
        for ref_id, cov_array in self.coverage_arrays.items():
            for pos, depth in enumerate(cov_array):
                coverage_data.append({
                    'chromosome': ref_id,
                    'position': pos + 1,  # 1-based coordinates
                    'coverage': depth
                })
        
        coverage_df = pd.DataFrame(coverage_data)
        coverage_file = data_dir / "coverage_per_position.csv"
        coverage_df.to_csv(coverage_file, index=False)
        
        logging.info(f"Coverage data exported to {coverage_file}")
        
        # Export gap regions in BED format
        gap_regions = []
        for ref_id, cov_array in self.coverage_arrays.items():
            below_threshold = cov_array < self.min_coverage
            
            if not np.any(below_threshold):
                continue
            
            # Find gap boundaries
            gap_starts = []
            gap_ends = []
            
            in_gap = False
            gap_start = 0
            
            for i, is_gap in enumerate(below_threshold):
                if is_gap and not in_gap:
                    gap_start = i
                    in_gap = True
                elif not is_gap and in_gap:
                    gap_regions.append({
                        'chromosome': ref_id,
                        'start': gap_start,
                        'end': i,
                        'length': i - gap_start
                    })
                    in_gap = False
            
            # Handle gap at end
            if in_gap:
                gap_regions.append({
                    'chromosome': ref_id,
                    'start': gap_start,
                    'end': len(cov_array),
                    'length': len(cov_array) - gap_start
                })
        
        if gap_regions:
            gap_df = pd.DataFrame(gap_regions)
            gap_file = data_dir / "gap_regions.bed"
            
            # Save in BED format (0-based)
            with open(gap_file, 'w') as f:
                for _, row in gap_df.iterrows():
                    f.write(f"{row['chromosome']}\t{row['start']}\t{row['end']}\n")
            
            logging.info(f"Gap regions exported to {gap_file}")
        
        return coverage_df, gap_regions if gap_regions else pd.DataFrame()