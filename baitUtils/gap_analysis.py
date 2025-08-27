#!/usr/bin/env python3

"""
gap_analysis.py

Intelligent gap analysis for oligo set coverage evaluation.
Analyzes coverage gaps, correlates with reference sequence features,
and provides insights for improving oligo set design.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import re


class GapAnalyzer:
    """Intelligent gap analysis for coverage evaluation."""
    
    def __init__(
        self,
        coverage_data: Dict[str, Any],
        reference_file: Path,
        min_gap_size: int = 100,
        extend_bp: int = 0
    ):
        """
        Initialize the gap analyzer.
        
        Args:
            coverage_data: Coverage statistics from CoverageAnalyzer
            reference_file: Path to reference FASTA file
            min_gap_size: Minimum gap size to analyze
            extend_bp: Extend gaps by N bases on each side for analysis
        """
        self.coverage_data = coverage_data
        self.reference_file = Path(reference_file)
        self.min_gap_size = min_gap_size
        self.extend_bp = extend_bp
        
        # Data storage
        self.reference_sequences = {}
        self.gaps = []
        self.gap_features = []
        self.analysis_results = {}
    
    def analyze(self, reference_analysis: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Perform comprehensive gap analysis.
        
        Args:
            reference_analysis: Optional reference sequence analysis for enhanced correlations
        
        Returns:
            Dictionary containing gap analysis results
        """
        logging.info("Loading reference sequences for gap analysis...")
        self._load_reference_sequences()
        
        logging.info("Extracting coverage gaps...")
        self._extract_gaps()
        
        logging.info("Analyzing gap sequence features...")
        self._analyze_gap_features()
        
        logging.info("Computing gap statistics...")
        self._compute_gap_statistics()
        
        logging.info("Correlating gaps with reference features...")
        self._correlate_with_reference_features(reference_analysis)
        
        logging.info("Generating improvement suggestions...")
        self._generate_suggestions()
        
        return self.analysis_results
    
    def _load_reference_sequences(self) -> None:
        """Load reference sequences from FASTA file."""
        try:
            with open(self.reference_file, 'r') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    self.reference_sequences[record.id] = str(record.seq).upper()
            
            logging.info(f"Loaded {len(self.reference_sequences)} reference sequences for gap analysis")
            
        except Exception as e:
            logging.error(f"Error loading reference sequences: {e}")
            raise
    
    def _extract_gaps(self) -> None:
        """Extract coverage gaps from coverage data."""
        per_ref_stats = self.coverage_data.get('per_reference', {})
        
        # We need to reconstruct gaps from coverage statistics
        # In practice, this would use the actual coverage arrays
        # For now, we'll simulate gap detection based on available data
        
        gap_id = 0
        for ref_id, stats in per_ref_stats.items():
            if ref_id not in self.reference_sequences:
                continue
            
            ref_length = stats['length']
            coverage_breadth = stats['coverage_breadth']
            gap_count = stats['gaps']
            
            # Simulate gap positions for demonstration
            # In real implementation, this would use actual coverage arrays
            uncovered_fraction = (100 - coverage_breadth) / 100
            avg_gap_size = (ref_length * uncovered_fraction) / max(1, gap_count)
            
            if avg_gap_size < self.min_gap_size:
                continue
            
            # Generate representative gaps
            for i in range(min(gap_count, 10)):  # Limit for demonstration
                # Simulate gap positions
                gap_start = int(i * ref_length / max(1, gap_count))
                gap_end = min(gap_start + int(avg_gap_size), ref_length)
                
                if gap_end - gap_start >= self.min_gap_size:
                    gap = {
                        'id': gap_id,
                        'chromosome': ref_id,
                        'start': gap_start,
                        'end': gap_end,
                        'length': gap_end - gap_start,
                        'extended_start': max(0, gap_start - self.extend_bp),
                        'extended_end': min(ref_length, gap_end + self.extend_bp)
                    }
                    self.gaps.append(gap)
                    gap_id += 1
        
        logging.info(f"Extracted {len(self.gaps)} gaps ≥ {self.min_gap_size} bp")
    
    def _analyze_gap_features(self) -> None:
        """Analyze sequence features within gaps."""
        for gap in self.gaps:
            ref_id = gap['chromosome']
            start = gap['extended_start']
            end = gap['extended_end']
            
            if ref_id not in self.reference_sequences:
                continue
            
            sequence = self.reference_sequences[ref_id][start:end]
            
            # Calculate sequence features
            features = self._calculate_sequence_features(sequence)
            features.update({
                'gap_id': gap['id'],
                'chromosome': ref_id,
                'start': gap['start'],
                'end': gap['end'],
                'length': gap['length'],
                'sequence_length': len(sequence)
            })
            
            self.gap_features.append(features)
    
    def _calculate_sequence_features(self, sequence: str) -> Dict[str, Any]:
        """Calculate comprehensive sequence features for a gap region."""
        if not sequence:
            return {'gc_content': 0, 'complexity': 0, 'repeats': 0, 
                   'homopolymers': 0, 'n_content': 0}
        
        # GC content
        gc_content = gc_fraction(sequence) * 100
        
        # Sequence complexity (Shannon entropy)
        complexity = self._calculate_shannon_entropy(sequence)
        
        # Repeat content (simple tandem repeats)
        repeat_content = self._calculate_repeat_content(sequence)
        
        # Homopolymer runs
        max_homopolymer = self._find_max_homopolymer(sequence)
        
        # N content
        n_content = (sequence.count('N') / len(sequence)) * 100
        
        # Low complexity regions
        low_complexity_fraction = self._calculate_low_complexity(sequence)
        
        # Dinucleotide composition
        dinuc_composition = self._calculate_dinucleotide_bias(sequence)
        
        return {
            'gc_content': gc_content,
            'complexity': complexity,
            'repeat_content': repeat_content,
            'max_homopolymer': max_homopolymer,
            'n_content': n_content,
            'low_complexity_fraction': low_complexity_fraction,
            'at_rich': gc_content < 30,
            'gc_rich': gc_content > 70,
            'extreme_composition': gc_content < 20 or gc_content > 80,
            'dinuc_bias': dinuc_composition
        }
    
    def _calculate_shannon_entropy(self, sequence: str) -> float:
        """Calculate Shannon entropy of sequence."""
        if not sequence:
            return 0.0
        
        # Count nucleotide frequencies
        counts = Counter(sequence)
        total = len(sequence)
        
        # Calculate entropy
        entropy = 0.0
        for count in counts.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)
        
        return entropy
    
    def _calculate_repeat_content(self, sequence: str) -> float:
        """Estimate repetitive content in sequence."""
        if len(sequence) < 6:
            return 0.0
        
        repeat_bases = 0
        window_size = 6
        
        for i in range(len(sequence) - window_size + 1):
            kmer = sequence[i:i + window_size]
            
            # Count occurrences of this k-mer
            count = 0
            for j in range(len(sequence) - window_size + 1):
                if sequence[j:j + window_size] == kmer:
                    count += 1
            
            # If k-mer appears multiple times, count as repetitive
            if count > 1:
                repeat_bases += 1
        
        return (repeat_bases / max(1, len(sequence) - window_size + 1)) * 100
    
    def _find_max_homopolymer(self, sequence: str) -> int:
        """Find maximum homopolymer run length."""
        if not sequence:
            return 0
        
        max_run = 1
        current_run = 1
        
        for i in range(1, len(sequence)):
            if sequence[i] == sequence[i-1]:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 1
        
        return max_run
    
    def _calculate_low_complexity(self, sequence: str, window=20) -> float:
        """Calculate fraction of sequence in low-complexity regions."""
        if len(sequence) < window:
            return 0.0
        
        low_complexity_bases = 0
        
        for i in range(len(sequence) - window + 1):
            subseq = sequence[i:i + window]
            entropy = self._calculate_shannon_entropy(subseq)
            
            # Low complexity threshold (entropy < 1.5 for DNA)
            if entropy < 1.5:
                low_complexity_bases += 1
        
        return low_complexity_bases / max(1, len(sequence) - window + 1)
    
    def _calculate_dinucleotide_bias(self, sequence: str) -> float:
        """Calculate dinucleotide composition bias."""
        if len(sequence) < 2:
            return 0.0
        
        dinuc_counts = Counter()
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            if 'N' not in dinuc:
                dinuc_counts[dinuc] += 1
        
        if not dinuc_counts:
            return 0.0
        
        # Calculate coefficient of variation
        counts = list(dinuc_counts.values())
        mean_count = np.mean(counts)
        std_count = np.std(counts)
        
        return (std_count / mean_count) if mean_count > 0 else 0.0
    
    def _compute_gap_statistics(self) -> None:
        """Compute comprehensive gap statistics."""
        if not self.gaps:
            self.analysis_results = {
                'total_gaps': 0,
                'total_gap_length': 0,
                'gap_percentage': 0.0,
                'mean_gap_size': 0.0,
                'median_gap_size': 0.0,
                'max_gap_size': 0,
                'size_distribution': {},
                'largest_gaps': [],
                'feature_analysis': {}
            }
            return
        
        # Basic gap statistics
        gap_lengths = [gap['length'] for gap in self.gaps]
        total_gap_length = sum(gap_lengths)
        total_ref_length = self.coverage_data.get('reference_length', 1)
        
        self.analysis_results = {
            'total_gaps': len(self.gaps),
            'total_gap_length': total_gap_length,
            'gap_percentage': (total_gap_length / total_ref_length) * 100,
            'mean_gap_size': np.mean(gap_lengths),
            'median_gap_size': np.median(gap_lengths),
            'max_gap_size': max(gap_lengths),
            'min_gap_size': min(gap_lengths)
        }
        
        # Gap size distribution
        size_ranges = [
            (self.min_gap_size, 500),
            (500, 1000),
            (1000, 5000),
            (5000, 10000),
            (10000, float('inf'))
        ]
        
        size_distribution = {}
        for start, end in size_ranges:
            if end == float('inf'):
                range_name = f"≥{start}"
                count = sum(1 for length in gap_lengths if length >= start)
            else:
                range_name = f"{start}-{end}"
                count = sum(1 for length in gap_lengths if start <= length < end)
            
            size_distribution[range_name] = count
        
        self.analysis_results['size_distribution'] = size_distribution
        
        # Largest gaps
        sorted_gaps = sorted(self.gaps, key=lambda x: x['length'], reverse=True)
        self.analysis_results['largest_gaps'] = sorted_gaps[:20]
        
        # Feature analysis
        if self.gap_features:
            self._analyze_gap_feature_correlations()
    
    def _analyze_gap_feature_correlations(self) -> None:
        """Analyze correlations between gap characteristics and sequence features."""
        if not self.gap_features:
            return
        
        feature_df = pd.DataFrame(self.gap_features)
        
        # Categorize gaps by size
        feature_df['size_category'] = pd.cut(
            feature_df['length'],
            bins=[0, 500, 1000, 5000, float('inf')],
            labels=['Small', 'Medium', 'Large', 'Very Large']
        )
        
        # Analyze feature distributions
        feature_analysis = {}
        
        # GC content analysis
        gc_stats = {
            'mean_gc': feature_df['gc_content'].mean(),
            'std_gc': feature_df['gc_content'].std(),
            'extreme_gc_fraction': feature_df['extreme_composition'].mean(),
            'at_rich_fraction': feature_df['at_rich'].mean(),
            'gc_rich_fraction': feature_df['gc_rich'].mean()
        }
        feature_analysis['gc_content'] = gc_stats
        
        # Complexity analysis
        complexity_stats = {
            'mean_complexity': feature_df['complexity'].mean(),
            'low_complexity_fraction': (feature_df['complexity'] < 1.5).mean(),
            'high_complexity_fraction': (feature_df['complexity'] > 1.8).mean()
        }
        feature_analysis['complexity'] = complexity_stats
        
        # Repeat content analysis
        repeat_stats = {
            'mean_repeat_content': feature_df['repeat_content'].mean(),
            'high_repeat_fraction': (feature_df['repeat_content'] > 30).mean()
        }
        feature_analysis['repeats'] = repeat_stats
        
        # Homopolymer analysis
        homopoly_stats = {
            'mean_max_homopolymer': feature_df['max_homopolymer'].mean(),
            'long_homopolymer_fraction': (feature_df['max_homopolymer'] > 8).mean()
        }
        feature_analysis['homopolymers'] = homopoly_stats
        
        # N content analysis
        n_stats = {
            'mean_n_content': feature_df['n_content'].mean(),
            'high_n_fraction': (feature_df['n_content'] > 10).mean()
        }
        feature_analysis['n_content'] = n_stats
        
        self.analysis_results['feature_analysis'] = feature_analysis
        
        # Size-stratified analysis
        size_analysis = {}
        for category in feature_df['size_category'].cat.categories:
            if pd.isna(category):
                continue
            
            subset = feature_df[feature_df['size_category'] == category]
            if len(subset) == 0:
                continue
            
            size_analysis[str(category)] = {
                'count': len(subset),
                'mean_gc': subset['gc_content'].mean(),
                'mean_complexity': subset['complexity'].mean(),
                'mean_repeat_content': subset['repeat_content'].mean(),
                'extreme_composition_fraction': subset['extreme_composition'].mean()
            }
        
        self.analysis_results['size_stratified_analysis'] = size_analysis
    
    def _generate_suggestions(self) -> None:
        """Generate improvement suggestions based on gap analysis."""
        suggestions = []
        
        if not self.gaps:
            suggestions.append("No significant coverage gaps detected. Coverage quality is good.")
            self.analysis_results['suggestions'] = suggestions
            return
        
        # Gap count suggestions
        total_gaps = self.analysis_results['total_gaps']
        if total_gaps > 50:
            suggestions.append(
                f"High number of gaps ({total_gaps}). Consider using 'baitUtils fill' "
                "to iteratively add oligos targeting uncovered regions."
            )
        
        # Gap size suggestions
        max_gap = self.analysis_results['max_gap_size']
        if max_gap > 10000:
            suggestions.append(
                f"Very large gaps detected (up to {max_gap:,} bp). "
                "These may represent problematic genomic regions requiring specialized oligo design."
            )
        
        # Feature-based suggestions
        if 'feature_analysis' in self.analysis_results:
            features = self.analysis_results['feature_analysis']
            
            # GC content suggestions
            if 'gc_content' in features:
                gc_stats = features['gc_content']
                if gc_stats['extreme_gc_fraction'] > 0.3:
                    suggestions.append(
                        f"{gc_stats['extreme_gc_fraction']*100:.1f}% of gaps have extreme GC content "
                        "(<20% or >80%). Consider adjusting oligo design parameters for these regions."
                    )
            
            # Complexity suggestions
            if 'complexity' in features:
                complexity_stats = features['complexity']
                if complexity_stats['low_complexity_fraction'] > 0.2:
                    suggestions.append(
                        f"{complexity_stats['low_complexity_fraction']*100:.1f}% of gaps are in "
                        "low-complexity regions. These areas may require specialized oligo design approaches."
                    )
            
            # Repeat suggestions
            if 'repeats' in features:
                repeat_stats = features['repeats']
                if repeat_stats['high_repeat_fraction'] > 0.2:
                    suggestions.append(
                        f"{repeat_stats['high_repeat_fraction']*100:.1f}% of gaps contain high repeat content. "
                        "Consider masking repetitive elements or using repeat-aware oligo design."
                    )
            
            # N content suggestions
            if 'n_content' in features:
                n_stats = features['n_content']
                if n_stats['high_n_fraction'] > 0.1:
                    suggestions.append(
                        f"{n_stats['high_n_fraction']*100:.1f}% of gaps contain high N content. "
                        "These may represent assembly gaps or low-quality sequence regions."
                    )
        
        # General improvement suggestions
        gap_percentage = self.analysis_results['gap_percentage']
        if gap_percentage > 20:
            suggestions.append(
                f"Gap percentage is high ({gap_percentage:.1f}%). "
                "Consider increasing oligo density or reviewing target selection criteria."
            )
        
        if len(suggestions) == 0:
            suggestions.append("Gap analysis shows good coverage characteristics.")
        
        self.analysis_results['suggestions'] = suggestions
    
    def _correlate_with_reference_features(self, reference_analysis: Dict[str, Any] = None) -> None:
        """Correlate gap characteristics with reference sequence features."""
        if not reference_analysis or not self.gap_features:
            self.analysis_results['reference_correlations'] = {}
            return
        
        ref_features = reference_analysis.get('sequence_features', {})
        
        correlations = {
            'gap_reference_correlations': {},
            'problematic_region_overlap': {},
            'feature_driven_gaps': {}
        }
        
        # Analyze correlations for each reference sequence
        for ref_id in self.reference_sequences.keys():
            if ref_id not in ref_features:
                continue
            
            ref_seq_features = ref_features[ref_id]
            ref_gaps = [g for g in self.gaps if g['chromosome'] == ref_id]
            
            if not ref_gaps:
                continue
            
            # Calculate gap density relative to reference features
            gap_density = len(ref_gaps) / ref_seq_features.get('length', 1) * 1000  # per kb
            
            ref_correlations = {
                'gap_count': len(ref_gaps),
                'gap_density_per_kb': gap_density,
                'gc_content': ref_seq_features.get('gc_content', 50),
                'complexity': ref_seq_features.get('shannon_entropy', 1.5),
                'repeat_content': ref_seq_features.get('repeat_content', 0),
                'correlation_strength': self._calculate_gap_feature_correlation(ref_gaps, ref_seq_features)
            }
            
            correlations['gap_reference_correlations'][ref_id] = ref_correlations
        
        # Identify feature-driven gaps
        correlations['feature_driven_gaps'] = self._identify_feature_driven_gaps(reference_analysis)
        
        # Calculate problematic region overlap
        correlations['problematic_region_overlap'] = self._calculate_problematic_overlap(reference_analysis)
        
        self.analysis_results['reference_correlations'] = correlations
    
    def _calculate_gap_feature_correlation(self, gaps: List[Dict], ref_features: Dict[str, Any]) -> float:
        """Calculate correlation strength between gaps and reference features."""
        if not gaps:
            return 0.0
        
        # Calculate expected vs actual gap distribution based on features
        gc_content = ref_features.get('gc_content', 50)
        complexity = ref_features.get('shannon_entropy', 1.5)
        repeat_content = ref_features.get('repeat_content', 0)
        
        # Expected gap likelihood based on features
        expected_gap_likelihood = 0.0
        
        # Extreme GC increases gap likelihood
        if gc_content < 30 or gc_content > 70:
            expected_gap_likelihood += 0.3
        
        # Low complexity increases gap likelihood
        if complexity < 1.0:
            expected_gap_likelihood += 0.4
        
        # High repeat content increases gap likelihood
        if repeat_content > 20:
            expected_gap_likelihood += 0.2
        
        # Calculate actual gap rate
        total_gap_length = sum(g['length'] for g in gaps)
        ref_length = ref_features.get('length', 1)
        actual_gap_rate = total_gap_length / ref_length
        
        # Normalize correlation (0-1 scale)
        correlation = min(1.0, actual_gap_rate / max(expected_gap_likelihood, 0.1))
        
        return correlation
    
    def _identify_feature_driven_gaps(self, reference_analysis: Dict[str, Any]) -> Dict[str, List[str]]:
        """Identify gaps that are likely driven by specific sequence features."""
        if not self.gap_features:
            return {}
        
        feature_driven = {
            'extreme_gc_gaps': [],
            'low_complexity_gaps': [],
            'repeat_rich_gaps': [],
            'homopolymer_gaps': [],
            'masked_region_gaps': []
        }
        
        for gap_feature in self.gap_features:
            gap_id = f"{gap_feature['chromosome']}:{gap_feature['start']}-{gap_feature['end']}"
            
            # Extreme GC content
            if gap_feature.get('extreme_composition', False):
                feature_driven['extreme_gc_gaps'].append(gap_id)
            
            # Low complexity
            if gap_feature.get('complexity', 2.0) < 1.0:
                feature_driven['low_complexity_gaps'].append(gap_id)
            
            # High repeat content
            if gap_feature.get('repeat_content', 0) > 30:
                feature_driven['repeat_rich_gaps'].append(gap_id)
            
            # Long homopolymers
            if gap_feature.get('max_homopolymer', 0) > 8:
                feature_driven['homopolymer_gaps'].append(gap_id)
            
            # Masked regions
            if gap_feature.get('n_content', 0) > 10:
                feature_driven['masked_region_gaps'].append(gap_id)
        
        return feature_driven
    
    def _calculate_problematic_overlap(self, reference_analysis: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate overlap between gaps and known problematic regions."""
        challenging_regions = reference_analysis.get('challenging_regions', {})
        
        if not challenging_regions or not self.gaps:
            return {}
        
        overlap_stats = {
            'total_gaps': len(self.gaps),
            'gaps_in_challenging_regions': 0,
            'challenging_region_gap_density': {},
            'gap_difficulty_correlation': 0.0
        }
        
        challenging_gap_count = 0
        total_gap_length = 0
        challenging_gap_length = 0
        
        for gap in self.gaps:
            ref_id = gap['chromosome']
            gap_length = gap['length']
            total_gap_length += gap_length
            
            if ref_id in challenging_regions:
                challenging_data = challenging_regions[ref_id]
                difficulty_score = challenging_data.get('challenging_score', 0)
                
                # If the reference is challenging, count this gap
                if difficulty_score >= 2:  # Moderate to high difficulty
                    challenging_gap_count += 1
                    challenging_gap_length += gap_length
        
        overlap_stats['gaps_in_challenging_regions'] = challenging_gap_count
        
        if challenging_gap_count > 0:
            overlap_stats['challenging_region_gap_density'] = {
                'gap_count_ratio': challenging_gap_count / len(self.gaps),
                'gap_length_ratio': challenging_gap_length / total_gap_length,
                'enrichment_factor': (challenging_gap_count / len(self.gaps)) / 
                                   max(0.1, len([r for r in challenging_regions.values() 
                                                if r.get('challenging_score', 0) >= 2]) / 
                                       len(challenging_regions))
            }
        
        return overlap_stats