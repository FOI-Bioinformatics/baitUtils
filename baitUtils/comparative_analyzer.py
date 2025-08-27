#!/usr/bin/env python3

"""
comparative_analyzer.py

Comparative analysis framework for comparing multiple oligo sets against the same reference.
Enables side-by-side comparison of coverage performance, gap patterns, and quality metrics.

This module provides the core framework for Phase 3 comparative analysis functionality,
allowing users to evaluate multiple oligo set designs and identify the best performing options.
"""

import logging
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, NamedTuple
from pathlib import Path
from dataclasses import dataclass
import tempfile
import shutil

from baitUtils.coverage_stats import CoverageAnalyzer
from baitUtils.gap_analysis import GapAnalyzer
from baitUtils.reference_analyzer import ReferenceAnalyzer
from baitUtils.quality_scorer import QualityScorer, QualityScore
from baitUtils.benchmark import BenchmarkAnalyzer


@dataclass
class OligoSetResult:
    """Container for single oligo set analysis results."""
    name: str
    file_path: str
    coverage_stats: Dict
    gap_analysis: Dict
    quality_score: QualityScore
    benchmark_results: Optional[Dict] = None


@dataclass
class ComparisonMetrics:
    """Container for comparative metrics between oligo sets."""
    coverage_breadth_diff: float
    mean_depth_diff: float
    gap_count_diff: int
    quality_score_diff: float
    mapping_efficiency_diff: float
    uniformity_diff: float


class ComparativeAnalyzer:
    """
    Framework for comparing multiple oligo sets against the same reference sequence.
    
    Provides comprehensive comparison of coverage statistics, gap patterns,
    quality scores, and performance benchmarks across multiple oligo set designs.
    """
    
    def __init__(self, reference_file: str, output_dir: Path, 
                 min_identity: float = 90.0, min_length: int = 100,
                 min_coverage: float = 1.0, target_coverage: float = 10.0):
        """
        Initialize comparative analyzer.
        
        Args:
            reference_file: Path to reference FASTA file
            output_dir: Directory for output files
            min_identity: Minimum mapping identity threshold
            min_length: Minimum mapping length
            min_coverage: Minimum coverage depth threshold
            target_coverage: Target coverage depth
        """
        self.reference_file = reference_file
        self.output_dir = Path(output_dir)
        self.min_identity = min_identity
        self.min_length = min_length
        self.min_coverage = min_coverage
        self.target_coverage = target_coverage
        
        self.oligo_sets: List[OligoSetResult] = []
        self.reference_analysis: Optional[Dict] = None
        self.comparison_matrix: Optional[pd.DataFrame] = None
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def add_oligo_set(self, name: str, oligo_file: str, 
                      psl_file: Optional[str] = None) -> None:
        """
        Add an oligo set for comparative analysis.
        
        Args:
            name: Descriptive name for the oligo set
            oligo_file: Path to oligo FASTA file
            psl_file: Optional pre-computed PSL mapping file
        """
        logging.info(f"Adding oligo set '{name}' for comparison...")
        
        # Perform mapping if PSL file not provided
        if psl_file is None:
            with tempfile.TemporaryDirectory() as temp_dir:
                psl_file = self._perform_mapping(oligo_file, Path(temp_dir))
                
                # Analyze this oligo set
                result = self._analyze_oligo_set(name, oligo_file, psl_file)
        else:
            result = self._analyze_oligo_set(name, oligo_file, psl_file)
        
        self.oligo_sets.append(result)
        logging.info(f"Added oligo set '{name}' with quality score {result.quality_score.overall_score:.1f}")
    
    def _perform_mapping(self, oligo_file: str, temp_dir: Path) -> str:
        """Perform oligo mapping using pblat."""
        import subprocess
        
        psl_file = temp_dir / "mapping.psl"
        
        cmd = [
            'pblat',
            f'-minIdentity={self.min_identity}',
            f'-minScore=30',
            f'-minMatch=2',
            self.reference_file,
            oligo_file,
            str(psl_file)
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"pblat mapping failed: {e}")
        
        return str(psl_file)
    
    def _analyze_oligo_set(self, name: str, oligo_file: str, psl_file: str) -> OligoSetResult:
        """Perform complete analysis of a single oligo set."""
        
        # Coverage analysis
        coverage_analyzer = CoverageAnalyzer(
            psl_file=psl_file,
            reference_file=self.reference_file,
            min_coverage=self.min_coverage,
            target_coverage=self.target_coverage,
            min_identity=self.min_identity,
            min_length=self.min_length
        )
        coverage_stats = coverage_analyzer.analyze()
        
        # Gap analysis
        gap_analyzer = GapAnalyzer(
            coverage_data=coverage_stats,
            reference_file=self.reference_file
        )
        gap_analysis = gap_analyzer.analyze()
        
        # Reference analysis (shared across all oligo sets)
        if self.reference_analysis is None:
            ref_analyzer = ReferenceAnalyzer(self.reference_file)
            self.reference_analysis = ref_analyzer.analyze()
        
        # Quality scoring
        quality_scorer = QualityScorer(
            coverage_stats=coverage_stats,
            gap_analysis=gap_analysis,
            reference_analysis=self.reference_analysis
        )
        quality_score = quality_scorer.calculate_score()
        
        # Benchmarking
        benchmark_analyzer = BenchmarkAnalyzer(
            coverage_stats=coverage_stats,
            gap_analysis=gap_analysis,
            reference_analysis=self.reference_analysis,
            quality_score=quality_score
        )
        benchmark_results = benchmark_analyzer.run_full_benchmark()
        
        return OligoSetResult(
            name=name,
            file_path=oligo_file,
            coverage_stats=coverage_stats,
            gap_analysis=gap_analysis,
            quality_score=quality_score,
            benchmark_results=benchmark_results
        )
    
    def generate_comparison_matrix(self) -> pd.DataFrame:
        """Generate comparison matrix of key metrics across all oligo sets."""
        
        if len(self.oligo_sets) < 2:
            raise ValueError("Need at least 2 oligo sets for comparison")
        
        data = []
        for result in self.oligo_sets:
            row = {
                'Name': result.name,
                'Coverage_Breadth_%': result.coverage_stats.get('coverage_breadth', 0),
                'Mean_Depth_x': result.coverage_stats.get('mean_depth', 0),
                'Total_Gaps': result.gap_analysis.get('total_gaps', 0),
                'Largest_Gap_bp': result.gap_analysis.get('max_gap_size', 0),
                'Mapping_Efficiency_%': result.coverage_stats.get('mapping_efficiency', 0),
                'Gini_Coefficient': result.coverage_stats.get('gini_coefficient', 0),
                'Quality_Score': result.quality_score.overall_score,
                'Quality_Grade': result.quality_score.grade,
                'Total_Oligos': result.coverage_stats.get('total_oligos', 0),
                'Mapped_Oligos': result.coverage_stats.get('mapped_oligos', 0),
            }
            
            # Add benchmark efficiency ratios
            if result.benchmark_results:
                row['Coverage_Efficiency_%'] = result.benchmark_results['coverage_breadth'].efficiency_ratio * 100
                row['Depth_Efficiency_%'] = result.benchmark_results['depth_uniformity'].efficiency_ratio * 100
                row['Gap_Efficiency_%'] = result.benchmark_results['gap_reduction'].efficiency_ratio * 100
                row['Overall_Efficiency_%'] = result.benchmark_results['overall_quality'].efficiency_ratio * 100
            
            data.append(row)
        
        self.comparison_matrix = pd.DataFrame(data)
        return self.comparison_matrix
    
    def calculate_pairwise_comparisons(self) -> Dict[Tuple[str, str], ComparisonMetrics]:
        """Calculate pairwise comparison metrics between all oligo sets."""
        
        comparisons = {}
        
        for i, set1 in enumerate(self.oligo_sets):
            for j, set2 in enumerate(self.oligo_sets[i+1:], i+1):
                
                # Calculate differences
                coverage_diff = (set1.coverage_stats.get('coverage_breadth', 0) - 
                               set2.coverage_stats.get('coverage_breadth', 0))
                depth_diff = (set1.coverage_stats.get('mean_depth', 0) - 
                            set2.coverage_stats.get('mean_depth', 0))
                gap_diff = (set1.gap_analysis.get('total_gaps', 0) - 
                          set2.gap_analysis.get('total_gaps', 0))
                quality_diff = set1.quality_score.overall_score - set2.quality_score.overall_score
                mapping_diff = (set1.coverage_stats.get('mapping_efficiency', 0) - 
                              set2.coverage_stats.get('mapping_efficiency', 0))
                uniformity_diff = (set1.coverage_stats.get('gini_coefficient', 0) - 
                                 set2.coverage_stats.get('gini_coefficient', 0))
                
                comparison = ComparisonMetrics(
                    coverage_breadth_diff=coverage_diff,
                    mean_depth_diff=depth_diff,
                    gap_count_diff=gap_diff,
                    quality_score_diff=quality_diff,
                    mapping_efficiency_diff=mapping_diff,
                    uniformity_diff=uniformity_diff
                )
                
                comparisons[(set1.name, set2.name)] = comparison
        
        return comparisons
    
    def identify_best_performer(self, metric: str = 'quality_score') -> OligoSetResult:
        """
        Identify the best performing oligo set based on specified metric.
        
        Args:
            metric: Metric to optimize ('quality_score', 'coverage_breadth', 
                   'gap_count', 'mapping_efficiency', 'depth_uniformity')
        """
        if not self.oligo_sets:
            raise ValueError("No oligo sets added for comparison")
        
        if metric == 'quality_score':
            return max(self.oligo_sets, key=lambda x: x.quality_score.overall_score)
        elif metric == 'coverage_breadth':
            return max(self.oligo_sets, key=lambda x: x.coverage_stats.get('coverage_breadth', 0))
        elif metric == 'gap_count':
            return min(self.oligo_sets, key=lambda x: x.gap_analysis.get('total_gaps', float('inf')))
        elif metric == 'mapping_efficiency':
            return max(self.oligo_sets, key=lambda x: x.coverage_stats.get('mapping_efficiency', 0))
        elif metric == 'depth_uniformity':
            return min(self.oligo_sets, key=lambda x: x.coverage_stats.get('gini_coefficient', 1.0))
        else:
            raise ValueError(f"Unknown metric: {metric}")
    
    def generate_ranking(self, weights: Optional[Dict[str, float]] = None) -> List[Tuple[str, float]]:
        """
        Generate ranking of oligo sets based on weighted composite score.
        
        Args:
            weights: Custom weights for different metrics
        """
        if weights is None:
            weights = {
                'quality_score': 0.4,
                'coverage_breadth': 0.25,
                'mapping_efficiency': 0.2,
                'gap_penalty': 0.15  # Penalty for gaps (inverted)
            }
        
        scores = []
        for result in self.oligo_sets:
            # Calculate composite score
            score = (
                result.quality_score.overall_score * weights['quality_score'] +
                result.coverage_stats.get('coverage_breadth', 0) / 10 * weights['coverage_breadth'] +
                result.coverage_stats.get('mapping_efficiency', 0) / 10 * weights['mapping_efficiency'] +
                max(0, (100 - result.gap_analysis.get('total_gaps', 0)) / 10) * weights['gap_penalty']
            )
            scores.append((result.name, score))
        
        # Sort by score (descending)
        scores.sort(key=lambda x: x[1], reverse=True)
        return scores
    
    def analyze_gap_overlap(self) -> Dict[str, Dict]:
        """Analyze overlapping gap regions between oligo sets."""
        
        if len(self.oligo_sets) < 2:
            return {}
        
        gap_analysis = {}
        
        # Extract gap regions for each set
        all_gaps = {}
        for result in self.oligo_sets:
            gaps = result.gap_analysis.get('gaps', [])
            all_gaps[result.name] = gaps
        
        # Find overlapping gaps
        overlaps = {}
        set_names = list(all_gaps.keys())
        
        for i, name1 in enumerate(set_names):
            for name2 in set_names[i+1:]:
                overlapping_gaps = self._find_overlapping_gaps(
                    all_gaps[name1], all_gaps[name2]
                )
                overlaps[f"{name1}_vs_{name2}"] = overlapping_gaps
        
        # Find unique gaps (present in only one set)
        unique_gaps = {}
        for name in set_names:
            unique_gaps[name] = self._find_unique_gaps(name, all_gaps)
        
        gap_analysis['overlapping_gaps'] = overlaps
        gap_analysis['unique_gaps'] = unique_gaps
        gap_analysis['total_unique_regions'] = len(set(
            (gap['chromosome'], gap['start'], gap['end']) 
            for gaps in all_gaps.values() 
            for gap in gaps
        ))
        
        return gap_analysis
    
    def _find_overlapping_gaps(self, gaps1: List[Dict], gaps2: List[Dict]) -> List[Dict]:
        """Find gaps that overlap between two sets."""
        overlapping = []
        
        for gap1 in gaps1:
            for gap2 in gaps2:
                if (gap1['chromosome'] == gap2['chromosome'] and
                    not (gap1['end'] < gap2['start'] or gap2['end'] < gap1['start'])):
                    
                    # Calculate overlap
                    overlap_start = max(gap1['start'], gap2['start'])
                    overlap_end = min(gap1['end'], gap2['end'])
                    overlap_length = overlap_end - overlap_start
                    
                    overlapping.append({
                        'chromosome': gap1['chromosome'],
                        'start': overlap_start,
                        'end': overlap_end,
                        'length': overlap_length,
                        'gap1_length': gap1['length'],
                        'gap2_length': gap2['length']
                    })
        
        return overlapping
    
    def _find_unique_gaps(self, target_set: str, all_gaps: Dict[str, List[Dict]]) -> List[Dict]:
        """Find gaps unique to a specific set."""
        target_gaps = all_gaps[target_set]
        other_gaps = []
        
        for name, gaps in all_gaps.items():
            if name != target_set:
                other_gaps.extend(gaps)
        
        unique = []
        for gap in target_gaps:
            is_unique = True
            for other_gap in other_gaps:
                if (gap['chromosome'] == other_gap['chromosome'] and
                    not (gap['end'] < other_gap['start'] or other_gap['end'] < gap['start'])):
                    is_unique = False
                    break
            
            if is_unique:
                unique.append(gap)
        
        return unique
    
    def export_comparison_data(self) -> Dict[str, str]:
        """Export comparison data to various formats."""
        
        exported_files = {}
        
        # Comparison matrix
        if self.comparison_matrix is not None:
            matrix_file = self.output_dir / "comparison_matrix.csv"
            self.comparison_matrix.to_csv(matrix_file, index=False)
            exported_files['comparison_matrix'] = str(matrix_file)
        
        # Pairwise comparisons
        comparisons = self.calculate_pairwise_comparisons()
        comparison_file = self.output_dir / "pairwise_comparisons.csv"
        
        comparison_data = []
        for (set1, set2), metrics in comparisons.items():
            comparison_data.append({
                'Set1': set1,
                'Set2': set2,
                'Coverage_Breadth_Diff': metrics.coverage_breadth_diff,
                'Mean_Depth_Diff': metrics.mean_depth_diff,
                'Gap_Count_Diff': metrics.gap_count_diff,
                'Quality_Score_Diff': metrics.quality_score_diff,
                'Mapping_Efficiency_Diff': metrics.mapping_efficiency_diff,
                'Uniformity_Diff': metrics.uniformity_diff
            })
        
        pd.DataFrame(comparison_data).to_csv(comparison_file, index=False)
        exported_files['pairwise_comparisons'] = str(comparison_file)
        
        # Gap overlap analysis
        gap_analysis = self.analyze_gap_overlap()
        if gap_analysis:
            gap_file = self.output_dir / "gap_overlap_analysis.txt"
            with open(gap_file, 'w') as f:
                f.write("GAP OVERLAP ANALYSIS\n")
                f.write("=" * 50 + "\n\n")
                
                f.write(f"Total unique gap regions: {gap_analysis['total_unique_regions']}\n\n")
                
                f.write("UNIQUE GAPS BY SET:\n")
                for set_name, gaps in gap_analysis['unique_gaps'].items():
                    f.write(f"\n{set_name}: {len(gaps)} unique gaps\n")
                    for gap in gaps[:5]:  # Show first 5
                        f.write(f"  {gap['chromosome']}:{gap['start']}-{gap['end']} ({gap['length']} bp)\n")
                    if len(gaps) > 5:
                        f.write(f"  ... and {len(gaps) - 5} more\n")
            
            exported_files['gap_analysis'] = str(gap_file)
        
        logging.info(f"Exported comparison data to {len(exported_files)} files")
        return exported_files