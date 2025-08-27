#!/usr/bin/env python3

"""
reference_analyzer.py

Reference sequence feature analysis for enhanced coverage evaluation.
Analyzes reference sequence characteristics to identify challenging regions
and correlate with coverage performance.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
import re

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis


class ReferenceAnalyzer:
    """Analyze reference sequence features and their correlation with coverage."""
    
    def __init__(
        self,
        reference_file: Path,
        coverage_data: Dict[str, Any],
        window_size: int = 1000,
        overlap: int = 500
    ):
        """
        Initialize the reference analyzer.
        
        Args:
            reference_file: Path to reference FASTA file
            coverage_data: Coverage statistics from CoverageAnalyzer
            window_size: Window size for sliding window analysis
            overlap: Overlap between adjacent windows
        """
        self.reference_file = Path(reference_file)
        self.coverage_data = coverage_data
        self.window_size = window_size
        self.overlap = overlap
        self.step_size = window_size - overlap
        
        # Data storage
        self.reference_sequences = {}
        self.sequence_features = {}
        self.window_features = {}
        self.coverage_correlations = {}
        self.challenging_regions = {}
    
    def analyze(self) -> Dict[str, Any]:
        """
        Perform comprehensive reference sequence analysis.
        
        Returns:
            Dictionary containing reference analysis results
        """
        logging.info("Analyzing reference sequence features...")
        
        self._load_reference_sequences()
        self._analyze_sequence_features()
        self._perform_window_analysis()
        self._correlate_with_coverage()
        self._identify_challenging_regions()
        
        return {
            'sequence_features': self.sequence_features,
            'window_features': self.window_features,
            'coverage_correlations': self.coverage_correlations,
            'challenging_regions': self.challenging_regions,
            'analysis_summary': self._generate_summary()
        }
    
    def _load_reference_sequences(self) -> None:
        """Load and store reference sequences."""
        try:
            with open(self.reference_file, 'r') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    self.reference_sequences[record.id] = str(record.seq).upper()
            
            logging.info(f"Loaded {len(self.reference_sequences)} reference sequences")
            
        except Exception as e:
            logging.error(f"Error loading reference sequences: {e}")
            raise
    
    def _analyze_sequence_features(self) -> None:
        """Analyze comprehensive features for each reference sequence."""
        for ref_id, sequence in self.reference_sequences.items():
            features = self._calculate_comprehensive_features(sequence)
            self.sequence_features[ref_id] = features
        
        logging.info("Sequence-level features analyzed")
    
    def _calculate_comprehensive_features(self, sequence: str) -> Dict[str, Any]:
        """Calculate comprehensive sequence features."""
        if not sequence:
            return self._empty_features()
        
        features = {}
        
        # Basic composition
        features.update(self._analyze_composition(sequence))
        
        # Complexity metrics
        features.update(self._analyze_complexity(sequence))
        
        # Structural features
        features.update(self._analyze_structure(sequence))
        
        # Repetitive elements
        features.update(self._analyze_repeats(sequence))
        
        # Problematic regions
        features.update(self._analyze_problematic_regions(sequence))
        
        return features
    
    def _analyze_composition(self, sequence: str) -> Dict[str, float]:
        """Analyze nucleotide composition features."""
        if not sequence:
            return {'gc_content': 0, 'at_content': 0, 'n_content': 0, 'length': 0}
        
        length = len(sequence)
        
        # Base composition
        gc_content = gc_fraction(sequence) * 100
        at_content = 100 - gc_content
        n_content = (sequence.count('N') / length) * 100
        
        # Dinucleotide composition
        dinuc_counts = Counter()
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            if 'N' not in dinuc:
                dinuc_counts[dinuc] += 1
        
        # CpG content (important for many analyses)
        cpg_observed = dinuc_counts.get('CG', 0)
        cpg_expected = (sequence.count('C') * sequence.count('G')) / length if length > 0 else 0
        cpg_ratio = (cpg_observed / cpg_expected) if cpg_expected > 0 else 0
        
        # Dinucleotide skew
        dinuc_skew = self._calculate_dinucleotide_skew(sequence)
        
        return {
            'gc_content': gc_content,
            'at_content': at_content,
            'n_content': n_content,
            'length': length,
            'cpg_ratio': cpg_ratio,
            'cpg_observed': cpg_observed,
            'cpg_expected': cpg_expected,
            'dinuc_skew': dinuc_skew
        }
    
    def _analyze_complexity(self, sequence: str) -> Dict[str, float]:
        """Analyze sequence complexity metrics."""
        if not sequence:
            return {'shannon_entropy': 0, 'linguistic_complexity': 0, 
                   'effective_length': 0, 'compression_ratio': 0}
        
        # Shannon entropy
        shannon_entropy = self._calculate_shannon_entropy(sequence)
        
        # Linguistic complexity (approximate)
        linguistic_complexity = self._calculate_linguistic_complexity(sequence)
        
        # Effective sequence length (non-repetitive content)
        effective_length = self._estimate_effective_length(sequence)
        
        # Compression ratio (estimate of redundancy)
        compression_ratio = effective_length / len(sequence) if len(sequence) > 0 else 0
        
        return {
            'shannon_entropy': shannon_entropy,
            'linguistic_complexity': linguistic_complexity,
            'effective_length': effective_length,
            'compression_ratio': compression_ratio
        }
    
    def _analyze_structure(self, sequence: str) -> Dict[str, Any]:
        """Analyze structural features of the sequence."""
        if not sequence:
            return {'max_homopolymer': 0, 'homopolymer_regions': [], 
                   'tandem_repeats': 0, 'inverted_repeats': 0}
        
        # Homopolymer analysis
        max_homopolymer, homopolymer_regions = self._analyze_homopolymers(sequence)
        
        # Simple tandem repeat detection
        tandem_repeats = self._count_tandem_repeats(sequence)
        
        # Inverted repeat detection (palindromes)
        inverted_repeats = self._count_inverted_repeats(sequence)
        
        # Secondary structure potential (simplified)
        structure_potential = self._estimate_structure_potential(sequence)
        
        return {
            'max_homopolymer': max_homopolymer,
            'homopolymer_regions': len(homopolymer_regions),
            'tandem_repeats': tandem_repeats,
            'inverted_repeats': inverted_repeats,
            'structure_potential': structure_potential
        }
    
    def _analyze_repeats(self, sequence: str) -> Dict[str, Any]:
        """Analyze repetitive elements in the sequence."""
        if not sequence:
            return {'repeat_content': 0, 'simple_repeats': 0, 'complex_repeats': 0}
        
        # Simple repeat content
        simple_repeats = self._count_simple_repeats(sequence)
        
        # Complex repeat patterns
        complex_repeats = self._count_complex_repeats(sequence)
        
        # Overall repeat content estimate
        repeat_content = (simple_repeats + complex_repeats) / len(sequence) * 100
        
        return {
            'repeat_content': repeat_content,
            'simple_repeats': simple_repeats,
            'complex_repeats': complex_repeats
        }
    
    def _analyze_problematic_regions(self, sequence: str) -> Dict[str, Any]:
        """Identify potentially problematic regions for oligo design."""
        if not sequence:
            return {'extreme_gc_regions': 0, 'low_complexity_regions': 0,
                   'masked_regions': 0, 'problematic_score': 0}
        
        # Extreme GC content regions (sliding window)
        extreme_gc_regions = self._count_extreme_gc_regions(sequence)
        
        # Low complexity regions
        low_complexity_regions = self._count_low_complexity_regions(sequence)
        
        # Masked regions (lowercase or N's)
        masked_regions = self._count_masked_regions(sequence)
        
        # Overall problematic score
        total_problematic = extreme_gc_regions + low_complexity_regions + masked_regions
        problematic_score = total_problematic / len(sequence) * 100
        
        return {
            'extreme_gc_regions': extreme_gc_regions,
            'low_complexity_regions': low_complexity_regions,
            'masked_regions': masked_regions,
            'problematic_score': problematic_score
        }
    
    def _perform_window_analysis(self) -> None:
        """Perform sliding window analysis of sequence features."""
        for ref_id, sequence in self.reference_sequences.items():
            windows = []
            
            for start in range(0, len(sequence) - self.window_size + 1, self.step_size):
                end = start + self.window_size
                window_seq = sequence[start:end]
                
                window_features = {
                    'start': start,
                    'end': end,
                    'length': len(window_seq)
                }
                
                # Calculate features for this window
                window_features.update(self._calculate_window_features(window_seq))
                windows.append(window_features)
            
            self.window_features[ref_id] = windows
        
        logging.info("Window-based analysis completed")
    
    def _calculate_window_features(self, window_seq: str) -> Dict[str, float]:
        """Calculate features for a sequence window."""
        if not window_seq:
            return {'gc_content': 0, 'complexity': 0, 'repeat_density': 0}
        
        # Basic features for the window
        gc_content = gc_fraction(window_seq) * 100
        complexity = self._calculate_shannon_entropy(window_seq)
        
        # Repeat density in window
        repeat_density = self._calculate_repeat_density(window_seq)
        
        # Homopolymer density
        homopolymer_density = self._calculate_homopolymer_density(window_seq)
        
        return {
            'gc_content': gc_content,
            'complexity': complexity,
            'repeat_density': repeat_density,
            'homopolymer_density': homopolymer_density
        }
    
    def _correlate_with_coverage(self) -> None:
        """Correlate sequence features with coverage performance."""
        per_ref_stats = self.coverage_data.get('per_reference', {})
        
        for ref_id in self.reference_sequences.keys():
            if ref_id not in per_ref_stats:
                continue
            
            ref_coverage_stats = per_ref_stats[ref_id]
            ref_features = self.sequence_features.get(ref_id, {})
            
            correlations = self._calculate_feature_correlations(ref_features, ref_coverage_stats)
            self.coverage_correlations[ref_id] = correlations
        
        logging.info("Coverage-feature correlations calculated")
    
    def _calculate_feature_correlations(
        self, 
        features: Dict[str, Any], 
        coverage_stats: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Calculate correlations between sequence features and coverage."""
        correlations = {}
        
        # Coverage quality metrics
        coverage_breadth = coverage_stats.get('coverage_breadth', 0)
        mean_depth = coverage_stats.get('mean_depth', 0)
        gap_count = coverage_stats.get('gaps', 0)
        
        # Correlate with sequence features
        correlations.update({
            'gc_vs_coverage': self._correlate_gc_coverage(features.get('gc_content', 50), coverage_breadth),
            'complexity_vs_coverage': self._correlate_complexity_coverage(features.get('shannon_entropy', 1.5), coverage_breadth),
            'repeats_vs_gaps': self._correlate_repeats_gaps(features.get('repeat_content', 0), gap_count),
            'problematic_vs_coverage': self._correlate_problematic_coverage(features.get('problematic_score', 0), coverage_breadth)
        })
        
        return correlations
    
    def _identify_challenging_regions(self) -> None:
        """Identify regions that are challenging for oligo design."""
        for ref_id, features in self.sequence_features.items():
            challenging_score = 0
            reasons = []
            
            # Extreme GC content
            gc_content = features.get('gc_content', 50)
            if gc_content < 20 or gc_content > 80:
                challenging_score += 2
                reasons.append(f"Extreme GC content ({gc_content:.1f}%)")
            
            # Low complexity
            if features.get('shannon_entropy', 2) < 1.0:
                challenging_score += 2
                reasons.append("Low sequence complexity")
            
            # High repeat content
            if features.get('repeat_content', 0) > 30:
                challenging_score += 1
                reasons.append("High repetitive content")
            
            # Long homopolymers
            if features.get('max_homopolymer', 0) > 10:
                challenging_score += 1
                reasons.append("Long homopolymer runs")
            
            # High N content
            if features.get('n_content', 0) > 5:
                challenging_score += 1
                reasons.append("High N content")
            
            # Problematic regions
            if features.get('problematic_score', 0) > 20:
                challenging_score += 1
                reasons.append("Multiple problematic features")
            
            self.challenging_regions[ref_id] = {
                'challenging_score': challenging_score,
                'reasons': reasons,
                'difficulty_level': self._get_difficulty_level(challenging_score)
            }
    
    def _get_difficulty_level(self, score: int) -> str:
        """Convert challenging score to difficulty level."""
        if score >= 5:
            return "Very Difficult"
        elif score >= 3:
            return "Difficult"
        elif score >= 1:
            return "Moderate"
        else:
            return "Easy"
    
    def _generate_summary(self) -> Dict[str, Any]:
        """Generate summary of reference analysis."""
        if not self.sequence_features:
            return {}
        
        # Aggregate statistics
        all_features = list(self.sequence_features.values())
        
        summary = {
            'total_sequences': len(all_features),
            'mean_gc_content': np.mean([f.get('gc_content', 50) for f in all_features]),
            'mean_complexity': np.mean([f.get('shannon_entropy', 1.5) for f in all_features]),
            'sequences_with_extreme_gc': sum(1 for f in all_features 
                                           if f.get('gc_content', 50) < 20 or f.get('gc_content', 50) > 80),
            'sequences_with_high_repeats': sum(1 for f in all_features 
                                             if f.get('repeat_content', 0) > 30),
            'challenging_sequences': sum(1 for r in self.challenging_regions.values() 
                                       if r.get('challenging_score', 0) >= 3)
        }
        
        return summary
    
    # Helper methods for feature calculations
    def _calculate_shannon_entropy(self, sequence: str) -> float:
        """Calculate Shannon entropy of sequence."""
        if not sequence:
            return 0.0
        
        counts = Counter(sequence)
        total = len(sequence)
        
        entropy = 0.0
        for count in counts.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)
        
        return entropy
    
    def _calculate_dinucleotide_skew(self, sequence: str) -> float:
        """Calculate dinucleotide composition skew."""
        if len(sequence) < 2:
            return 0.0
        
        dinuc_counts = Counter()
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            if 'N' not in dinuc:
                dinuc_counts[dinuc] += 1
        
        if not dinuc_counts:
            return 0.0
        
        counts = list(dinuc_counts.values())
        return np.std(counts) / np.mean(counts) if np.mean(counts) > 0 else 0.0
    
    def _calculate_linguistic_complexity(self, sequence: str) -> float:
        """Calculate linguistic complexity (simplified)."""
        if len(sequence) < 10:
            return 0.0
        
        # Count unique k-mers of different lengths
        complexity_scores = []
        for k in [2, 3, 4]:
            if len(sequence) >= k:
                kmers = set()
                for i in range(len(sequence) - k + 1):
                    kmers.add(sequence[i:i+k])
                
                max_possible = min(4**k, len(sequence) - k + 1)
                complexity_scores.append(len(kmers) / max_possible)
        
        return np.mean(complexity_scores) if complexity_scores else 0.0
    
    def _estimate_effective_length(self, sequence: str) -> float:
        """Estimate effective (non-repetitive) sequence length."""
        if len(sequence) < 20:
            return len(sequence)
        
        # Use 6-mer frequency to estimate repetitiveness
        kmers = Counter()
        k = 6
        
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if 'N' not in kmer:
                kmers[kmer] += 1
        
        # Count unique positions
        unique_positions = sum(1 for count in kmers.values() if count == 1) * k
        return max(unique_positions, len(sequence) * 0.1)  # At least 10% unique
    
    def _analyze_homopolymers(self, sequence: str) -> Tuple[int, List[Tuple[int, int, str]]]:
        """Analyze homopolymer runs in sequence."""
        if not sequence:
            return 0, []
        
        max_length = 0
        regions = []
        
        current_base = sequence[0]
        current_length = 1
        start_pos = 0
        
        for i in range(1, len(sequence)):
            if sequence[i] == current_base:
                current_length += 1
            else:
                if current_length >= 5:  # Report homopolymers >= 5bp
                    regions.append((start_pos, start_pos + current_length, current_base))
                    max_length = max(max_length, current_length)
                
                current_base = sequence[i]
                current_length = 1
                start_pos = i
        
        # Handle last homopolymer
        if current_length >= 5:
            regions.append((start_pos, start_pos + current_length, current_base))
            max_length = max(max_length, current_length)
        
        return max_length, regions
    
    def _count_tandem_repeats(self, sequence: str) -> int:
        """Count simple tandem repeats."""
        if len(sequence) < 6:
            return 0
        
        repeat_count = 0
        
        # Look for 2-6 bp repeats
        for unit_size in range(2, 7):
            for i in range(len(sequence) - unit_size * 2 + 1):
                unit = sequence[i:i+unit_size]
                if 'N' in unit:
                    continue
                
                # Check if this unit repeats
                repeats = 1
                pos = i + unit_size
                
                while pos + unit_size <= len(sequence) and sequence[pos:pos+unit_size] == unit:
                    repeats += 1
                    pos += unit_size
                
                if repeats >= 3:  # At least 3 copies
                    repeat_count += repeats * unit_size
        
        return repeat_count
    
    def _count_inverted_repeats(self, sequence: str) -> int:
        """Count inverted repeats (palindromes)."""
        if len(sequence) < 6:
            return 0
        
        complement = str.maketrans('ATGC', 'TACG')
        inverted_count = 0
        
        # Look for palindromes of 6+ bp
        for length in range(6, min(21, len(sequence) // 2 + 1)):
            for i in range(len(sequence) - length + 1):
                subseq = sequence[i:i+length]
                if 'N' in subseq:
                    continue
                
                # Check if it's a palindrome
                reverse_complement = subseq.translate(complement)[::-1]
                if subseq == reverse_complement:
                    inverted_count += 1
        
        return inverted_count
    
    def _estimate_structure_potential(self, sequence: str) -> float:
        """Estimate secondary structure formation potential."""
        if len(sequence) < 20:
            return 0.0
        
        # Simple estimate based on GC content and palindromes
        gc_content = gc_fraction(sequence)
        palindromes = self._count_inverted_repeats(sequence)
        
        # Higher GC content and more palindromes = higher structure potential
        structure_score = (gc_content * 0.7) + (palindromes / len(sequence) * 100 * 0.3)
        
        return min(structure_score, 1.0)
    
    def _count_simple_repeats(self, sequence: str) -> int:
        """Count simple repetitive elements."""
        # This is a simplified implementation
        return self._count_tandem_repeats(sequence)
    
    def _count_complex_repeats(self, sequence: str) -> int:
        """Count complex repetitive patterns."""
        if len(sequence) < 20:
            return 0
        
        # Look for longer repeated sequences
        complex_repeats = 0
        
        for length in range(10, min(51, len(sequence) // 2)):
            seen = set()
            for i in range(len(sequence) - length + 1):
                subseq = sequence[i:i+length]
                if 'N' in subseq:
                    continue
                
                if subseq in seen:
                    complex_repeats += length
                else:
                    seen.add(subseq)
        
        return complex_repeats
    
    def _count_extreme_gc_regions(self, sequence: str, window=50) -> int:
        """Count regions with extreme GC content."""
        if len(sequence) < window:
            return 0
        
        extreme_regions = 0
        
        for i in range(len(sequence) - window + 1):
            subseq = sequence[i:i+window]
            gc_content = gc_fraction(subseq) * 100
            
            if gc_content < 20 or gc_content > 80:
                extreme_regions += 1
        
        return extreme_regions
    
    def _count_low_complexity_regions(self, sequence: str, window=50) -> int:
        """Count low complexity regions."""
        if len(sequence) < window:
            return 0
        
        low_complexity = 0
        
        for i in range(len(sequence) - window + 1):
            subseq = sequence[i:i+window]
            entropy = self._calculate_shannon_entropy(subseq)
            
            if entropy < 1.0:  # Very low complexity
                low_complexity += 1
        
        return low_complexity
    
    def _count_masked_regions(self, sequence: str) -> int:
        """Count masked or N regions."""
        return sequence.count('N') + sum(1 for c in sequence if c.islower())
    
    def _calculate_repeat_density(self, sequence: str) -> float:
        """Calculate repeat density for a sequence window."""
        if len(sequence) < 10:
            return 0.0
        
        repeats = self._count_simple_repeats(sequence) + self._count_complex_repeats(sequence)
        return repeats / len(sequence) * 100
    
    def _calculate_homopolymer_density(self, sequence: str) -> float:
        """Calculate homopolymer density for a sequence window."""
        if len(sequence) < 5:
            return 0.0
        
        homopolymer_bases = 0
        current_run = 1
        
        for i in range(1, len(sequence)):
            if sequence[i] == sequence[i-1]:
                current_run += 1
            else:
                if current_run >= 4:  # Count runs of 4+
                    homopolymer_bases += current_run
                current_run = 1
        
        # Handle final run
        if current_run >= 4:
            homopolymer_bases += current_run
        
        return homopolymer_bases / len(sequence) * 100
    
    def _correlate_gc_coverage(self, gc_content: float, coverage_breadth: float) -> Dict[str, float]:
        """Correlate GC content with coverage performance."""
        # Optimal GC range is typically 40-60%
        optimal_distance = min(abs(gc_content - 40), abs(gc_content - 60))
        if 40 <= gc_content <= 60:
            optimal_distance = 0
        
        return {
            'gc_content': gc_content,
            'coverage_breadth': coverage_breadth,
            'optimal_distance': optimal_distance,
            'correlation_strength': max(0, 1 - (optimal_distance / 30))  # Normalize
        }
    
    def _correlate_complexity_coverage(self, complexity: float, coverage_breadth: float) -> Dict[str, float]:
        """Correlate sequence complexity with coverage."""
        # Higher complexity generally better for oligo design
        return {
            'complexity': complexity,
            'coverage_breadth': coverage_breadth,
            'complexity_score': min(complexity / 2.0, 1.0)  # Normalize to 0-1
        }
    
    def _correlate_repeats_gaps(self, repeat_content: float, gap_count: int) -> Dict[str, float]:
        """Correlate repeat content with coverage gaps."""
        return {
            'repeat_content': repeat_content,
            'gap_count': gap_count,
            'negative_correlation': repeat_content / 100  # Higher repeats = more gaps
        }
    
    def _correlate_problematic_coverage(self, problematic_score: float, coverage_breadth: float) -> Dict[str, float]:
        """Correlate problematic features with coverage."""
        return {
            'problematic_score': problematic_score,
            'coverage_breadth': coverage_breadth,
            'impact_score': problematic_score / 100  # Higher problematic = lower coverage
        }
    
    def _empty_features(self) -> Dict[str, Any]:
        """Return empty feature dictionary for invalid sequences."""
        return {
            'gc_content': 0, 'at_content': 0, 'n_content': 0, 'length': 0,
            'shannon_entropy': 0, 'linguistic_complexity': 0,
            'max_homopolymer': 0, 'repeat_content': 0,
            'problematic_score': 0
        }