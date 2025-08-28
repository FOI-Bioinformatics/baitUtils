#!/usr/bin/env python3

"""
gap_filling_algorithm.py

Core gap filling algorithms for multi-pass oligo selection to maximize coverage.
Provides greedy selection strategies and coverage optimization logic.
"""

import logging
import math
import bisect
from typing import Dict, List, Tuple, Set, Optional
from collections import defaultdict
from dataclasses import dataclass


@dataclass
class OligoMapping:
    """Represents a single oligo mapping."""
    oligo_id: str
    ref_id: str
    start: int
    end: int
    coverage: float = 1.0


class CoveragePatternAnalyzer:
    """Analyzes coverage patterns to identify difficult regions."""
    
    @staticmethod
    def analyze_coverage_patterns(
        uncovered: Dict[str, List[Tuple[int, int]]],
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


class OligoScorer:
    """Scores oligos for selection based on multiple criteria."""
    
    @staticmethod
    def calculate_gc_content(sequence: str) -> float:
        """Calculate GC content of a sequence."""
        if not sequence:
            return 0.0
        gc_count = sequence.count('G') + sequence.count('C')
        return gc_count / len(sequence)
    
    @staticmethod
    def score_oligo(
        m: OligoMapping,
        uncovered: Dict[str, List[Tuple[int, int]]],
        min_contribution: int,
        unique_region_bonus: float = 2.0,
        coverage_threshold: float = 0.8,
        difficulty_scores: Optional[Dict[str, float]] = None,
        sequence: Optional[str] = None,
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
        
        base_score = (total_overlap + length_weighted_overlap) * \
                    (1.0 + coverage_efficiency + length_bonus)

        if difficulty_bonus > 0:
            base_score *= (1.0 + difficulty_bonus)
        
        # Sequence-based GC scoring
        if sequence:
            gc_content = OligoScorer.calculate_gc_content(sequence)
            gc_penalty = abs(gc_content - optimal_gc)
            sequence_score = 1.0 - (gc_weight * gc_penalty)
            base_score *= sequence_score
        
        # Unique coverage bonus
        unique_coverage_bonus = 0
        for overlap, region_length in overlap_regions:
            if region_length < 2 * oligo_length:
                coverage_ratio = overlap / region_length
                if coverage_ratio > coverage_threshold:
                    bonus_scale = coverage_ratio * unique_region_bonus * \
                                 math.log2(region_length)
                    unique_coverage_bonus += overlap * bonus_scale

        return base_score + unique_coverage_bonus


class SpacingConstraintChecker:
    """Checks spacing constraints for oligo selection."""
    
    @staticmethod
    def can_select_oligo(
        selected_positions: Dict[str, List[Tuple[int, int]]],
        candidate: OligoMapping,
        spacing: int
    ) -> bool:
        """Check if candidate oligo respects the required spacing from previously selected oligos."""
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


class GreedySelector:
    """Implements greedy selection algorithm for gap filling."""
    
    def __init__(self):
        """Initialize greedy selector."""
        self.pattern_analyzer = CoveragePatternAnalyzer()
        self.oligo_scorer = OligoScorer()
        self.spacing_checker = SpacingConstraintChecker()
    
    def single_pass_greedy(
        self,
        all_mappings: Dict[str, List[OligoMapping]],
        uncovered: Dict[str, List[Tuple[int, int]]],
        selected_oligos: Set[str],
        spacing_distance: int,
        min_contribution: int,
        sequences: Optional[Dict[str, str]] = None,
        max_oligos_per_pass: Optional[int] = None
    ) -> Set[str]:
        """Perform a single pass of greedy selection based on scoring."""
        selected_positions = defaultdict(list)
        for ref_id, mlist in all_mappings.items():
            for m in mlist:
                if m.oligo_id in selected_oligos:
                    bisect.insort(selected_positions[ref_id], (m.start, m.end))

        difficulty_scores = self.pattern_analyzer.analyze_coverage_patterns(
            uncovered, all_mappings
        )

        candidate_list = []
        for ref_id, mlist in all_mappings.items():
            for m in mlist:
                if m.oligo_id in selected_oligos:
                    continue

                seq_slice = None
                if sequences and ref_id in sequences:
                    seq_slice = sequences[ref_id][m.start:m.end]

                sc = self.oligo_scorer.score_oligo(
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

            if self.spacing_checker.can_select_oligo(
                selected_positions, cand_oligo, spacing_distance
            ):
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


class MultiPassSelector:
    """Implements multi-pass selection algorithm for gap filling."""
    
    def __init__(self):
        """Initialize multi-pass selector."""
        self.greedy_selector = GreedySelector()
    
    def multi_pass_selection(
        self,
        all_mappings: Dict[str, List[OligoMapping]],
        forced_oligos: Set[str],
        coverage_calculator,  # Function to calculate coverage
        min_coverage: float,
        max_coverage: Optional[float],
        spacing_distance: int,
        min_contribution: int,
        max_passes: int,
        max_oligos_per_pass: Optional[int] = None,
        sequences: Optional[Dict[str, str]] = None,
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

            uncovered_regions, total_uncovered = coverage_calculator(
                selected_oligos, all_mappings, min_coverage, max_coverage
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
                    logging.info(f"No improvement in uncovered region count. "
                               f"Stall count: {stall_count}/{stall_rounds}")
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

            newly_selected = self.greedy_selector.single_pass_greedy(
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
            logging.info(f"Pass {pass_num}: selected {new_count - old_count} new oligos. "
                        f"Total now: {new_count}")

            pass_num += 1

        return selected_oligos