#!/usr/bin/env python3

"""
sequence_analysis.py

Core sequence analysis functions for oligo/bait analysis.
Provides utilities for calculating sequence properties, statistics, and quality metrics.
"""

import logging
import math
from typing import Union, List, Tuple
from collections import Counter
from math import log2

from Bio import Seq
from Bio.SeqUtils import gc_fraction, MeltingTemp as mt
from Bio.Align import PairwiseAligner

try:
    import RNA  # ViennaRNA package for MFE calculation
    HAS_VIENNA_RNA = True
except ImportError:
    HAS_VIENNA_RNA = False


class SequenceAnalyzer:
    """
    Comprehensive sequence analysis toolkit for oligonucleotide analysis.
    
    Provides methods for calculating various sequence properties including
    GC content, melting temperature, minimum free energy, complexity metrics,
    and secondary structure predictions.
    """
    
    def __init__(self, na_equivalent: str = '50mM', dnac1_equivalent: float = 250.0, 
                 dnac2_equivalent: float = 250.0):
        """
        Initialize sequence analyzer with thermodynamic parameters.
        
        Args:
            na_equivalent: Salt concentration for Tm calculations
            dnac1_equivalent: DNA concentration 1 for Tm calculations (nM)
            dnac2_equivalent: DNA concentration 2 for Tm calculations (nM)
        """
        self.na_equivalent = na_equivalent
        self.dnac1_equivalent = dnac1_equivalent
        self.dnac2_equivalent = dnac2_equivalent
        
        # Initialize pairwise aligner for self-alignment
        self.aligner = PairwiseAligner()
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -2
        self.aligner.extend_gap_score = -0.5
    
    def calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content percentage."""
        return gc_fraction(sequence) * 100.0
    
    def calculate_melting_temperature(self, sequence: str) -> Union[float, str]:
        """
        Calculate melting temperature using BioPython's MeltingTemp module.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Melting temperature in Celsius or 'NA' if calculation fails
        """
        try:
            tm = mt.Tm_NN(
                sequence, 
                Na=self.na_equivalent, 
                dnac1=self.dnac1_equivalent, 
                dnac2=self.dnac2_equivalent
            )
            return round(tm, 2)
        except Exception as e:
            logging.debug(f"Tm calculation error for sequence: {e}")
            return 'NA'
    
    def calculate_mfe(self, sequence: str) -> Union[float, str]:
        """
        Calculate Minimum Free Energy using ViennaRNA (if available).
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            MFE value in kcal/mol or 'NA' if calculation fails or ViennaRNA unavailable
        """
        if not HAS_VIENNA_RNA:
            return 'NA'
        
        try:
            md = RNA.md()
            md.material = 'DNA'  # Set to DNA parameters
            fc = RNA.fold_compound(sequence, md)
            structure, mfe = fc.mfe()
            return mfe
        except Exception as e:
            logging.debug(f"MFE calculation error for sequence: {e}")
            return 'NA'
    
    def calculate_self_alignment_score(self, sequence: str) -> Union[float, str]:
        """
        Calculate self-alignment score between sequence and its reverse complement.
        
        Args:
            sequence: DNA or RNA sequence
            
        Returns:
            Alignment score or 'NA' if calculation fails
        """
        try:
            seq_obj = Seq.Seq(sequence)
            reverse_complement = str(seq_obj.reverse_complement())
            
            alignments = self.aligner.align(sequence, reverse_complement)
            if alignments:
                return alignments[0].score
            return 0.0
        except Exception as e:
            logging.debug(f"Self-alignment calculation error: {e}")
            return 'NA'
    
    def calculate_entropy(self, sequence: str) -> float:
        """
        Calculate Shannon entropy of the sequence.
        
        Args:
            sequence: Input sequence
            
        Returns:
            Shannon entropy value
        """
        if not sequence:
            return 0.0
        
        # Count frequency of each character
        counts = Counter(sequence.upper())
        length = len(sequence)
        
        # Calculate entropy
        entropy = 0.0
        for count in counts.values():
            if count > 0:
                probability = count / length
                entropy -= probability * log2(probability)
        
        return entropy
    
    def calculate_complexity(self, sequence: str, k: int = 2) -> float:
        """
        Calculate sequence complexity based on k-mer diversity.
        
        Args:
            sequence: Input sequence
            k: K-mer length for complexity calculation
            
        Returns:
            Complexity score (0-1, higher = more complex)
        """
        if len(sequence) < k:
            return 0.0
        
        # Generate all k-mers
        kmers = []
        for i in range(len(sequence) - k + 1):
            kmers.append(sequence[i:i+k])
        
        # Calculate complexity as unique k-mers / total k-mers
        if not kmers:
            return 0.0
            
        unique_kmers = len(set(kmers))
        total_kmers = len(kmers)
        
        return unique_kmers / total_kmers
    
    def count_homopolymer_runs(self, sequence: str, min_length: int = 3) -> List[Tuple[str, int, int]]:
        """
        Count homopolymer runs (consecutive identical bases).
        
        Args:
            sequence: Input sequence
            min_length: Minimum run length to report
            
        Returns:
            List of (base, start_pos, length) tuples
        """
        runs = []
        if not sequence:
            return runs
        
        current_base = sequence[0]
        current_start = 0
        current_length = 1
        
        for i in range(1, len(sequence)):
            if sequence[i] == current_base:
                current_length += 1
            else:
                if current_length >= min_length:
                    runs.append((current_base, current_start, current_length))
                current_base = sequence[i]
                current_start = i
                current_length = 1
        
        # Check final run
        if current_length >= min_length:
            runs.append((current_base, current_start, current_length))
        
        return runs
    
    def calculate_dinucleotide_bias(self, sequence: str) -> float:
        """
        Calculate dinucleotide composition bias.
        
        Args:
            sequence: Input sequence
            
        Returns:
            Bias score (0 = uniform, higher = more biased)
        """
        if len(sequence) < 2:
            return 0.0
        
        # Count dinucleotides
        dinucleotides = []
        for i in range(len(sequence) - 1):
            dinucleotides.append(sequence[i:i+2])
        
        if not dinucleotides:
            return 0.0
        
        # Calculate frequencies
        counts = Counter(dinucleotides)
        frequencies = [count / len(dinucleotides) for count in counts.values()]
        
        # Calculate bias as coefficient of variation
        if not frequencies:
            return 0.0
        
        mean_freq = sum(frequencies) / len(frequencies)
        if mean_freq == 0:
            return 0.0
        
        variance = sum((f - mean_freq) ** 2 for f in frequencies) / len(frequencies)
        return math.sqrt(variance) / mean_freq
    
    def count_n_bases(self, sequence: str) -> int:
        """Count number of N bases in sequence."""
        return sequence.upper().count('N')
    
    def count_masked_bases(self, sequence: str) -> int:
        """Count number of lowercase (masked) bases."""
        return sum(1 for c in sequence if c.islower())
    
    def analyze_sequence(self, sequence: str, sequence_id: str = "") -> dict:
        """
        Perform comprehensive sequence analysis.
        
        Args:
            sequence: DNA sequence to analyze
            sequence_id: Optional sequence identifier
            
        Returns:
            Dictionary containing all calculated metrics
        """
        results = {
            'sequence_id': sequence_id,
            'length': len(sequence),
            'gc_content': self.calculate_gc_content(sequence),
            'melting_temperature': self.calculate_melting_temperature(sequence),
            'mfe': self.calculate_mfe(sequence),
            'entropy': self.calculate_entropy(sequence),
            'complexity_2mer': self.calculate_complexity(sequence, k=2),
            'complexity_3mer': self.calculate_complexity(sequence, k=3),
            'self_alignment_score': self.calculate_self_alignment_score(sequence),
            'dinucleotide_bias': self.calculate_dinucleotide_bias(sequence),
            'n_count': self.count_n_bases(sequence),
            'masked_count': self.count_masked_bases(sequence),
            'masked_percentage': (self.count_masked_bases(sequence) / len(sequence) * 100) if sequence else 0
        }
        
        # Add homopolymer analysis
        homopolymer_runs = self.count_homopolymer_runs(sequence)
        results['homopolymer_runs'] = len(homopolymer_runs)
        results['max_homopolymer_length'] = max([run[2] for run in homopolymer_runs], default=0)
        results['homopolymer_details'] = homopolymer_runs
        
        return results


def reverse_complement(sequence: str) -> str:
    """
    Calculate reverse complement of DNA sequence.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Reverse complement sequence
    """
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    new_seq = []
    
    for base in sequence.upper()[::-1]:
        if base in complement_map:
            new_seq.append(complement_map[base])
        else:
            new_seq.append(base)
    
    return ''.join(new_seq)


def clean_sequence(sequence: str, remove_ns: bool = False) -> str:
    """
    Clean sequence by removing non-ATCGN characters.
    
    Args:
        sequence: Input sequence
        remove_ns: If True, also remove N characters
        
    Returns:
        Cleaned sequence
    """
    valid_chars = set('ATCGN')
    if remove_ns:
        valid_chars.remove('N')
    
    return ''.join(c for c in sequence.upper() if c in valid_chars)


def validate_dna_sequence(sequence: str) -> Tuple[bool, str]:
    """
    Validate DNA sequence for allowed characters.
    
    Args:
        sequence: Sequence to validate
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not sequence:
        return False, "Empty sequence"
    
    valid_chars = set('ATCGNRYSWKMBDHV')  # Include IUPAC ambiguous codes
    invalid_chars = set(sequence.upper()) - valid_chars
    
    if invalid_chars:
        return False, f"Invalid characters found: {', '.join(sorted(invalid_chars))}"
    
    return True, "Valid DNA sequence"