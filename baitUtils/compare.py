#!/usr/bin/env python3

"""
compare.py

Main command for Phase 3 comparative analysis of multiple oligo sets.
Enables comprehensive comparison of oligo set performance with statistical analysis,
visualizations, and detailed reporting.

Usage:
    baitUtils compare -r reference.fasta -o comparison_report/ \
        --sets "Set1:oligos1.fasta" "Set2:oligos2.fasta" "Set3:oligos3.fasta"
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List

from baitUtils._version import __version__
from baitUtils.logging_utils import ensure_logging
from baitUtils.comparative_analyzer import ComparativeAnalyzer
from baitUtils.differential_analysis import DifferentialAnalyzer
from baitUtils.json_export import write_json, tool_versions, arguments_record
from baitUtils.mapping_utils import check_mapper_available
from baitUtils.comparative_visualizations import ComparativeVisualizer
from baitUtils.comparative_report_generator import ComparativeReportGenerator


def add_arguments(parser):
    """Add command-line arguments for the 'compare' subcommand."""
    parser.add_argument('--version', action='version', version=f"baitUtils compare {__version__}")
    
    # Required arguments
    parser.add_argument(
        '-r', '--reference',
        required=True,
        help='Reference genome/target FASTA file'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory for comparison results'
    )
    parser.add_argument(
        '--sets',
        nargs='+',
        required=True,
        help='Oligo sets to compare in format "Name:path/to/oligos.fasta"'
    )
    
    # Mapping parameters
    parser.add_argument(
        '--mapper',
        choices=['pblat', 'minimap2'],
        default='pblat',
        help='Mapping tool (default: pblat). minimap2 is run with -c and the chosen preset'
    )
    parser.add_argument(
        '--strand',
        choices=['both', 'plus', 'minus'],
        default='both',
        help='Count hits on both strands, or only plus or minus strand hits (default: both)'
    )
    parser.add_argument(
        '--minimap2-preset',
        default='sr',
        help='minimap2 -x preset when --mapper minimap2 (default: sr)'
    )
    parser.add_argument(
        '--min-identity',
        type=float,
        default=90.0,
        help='Minimum sequence identity for mapping (default: 90.0)'
    )
    parser.add_argument(
        '--min-length',
        type=int,
        default=100,
        help='Minimum mapping length (default: 100)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Number of threads for mapping (default: 1)'
    )
    
    # Coverage analysis parameters
    parser.add_argument(
        '--min-coverage',
        type=float,
        default=1.0,
        help='Minimum coverage depth threshold (default: 1.0)'
    )
    parser.add_argument(
        '--target-coverage',
        type=float,
        default=10.0,
        help='Target coverage depth for analysis (default: 10.0)'
    )
    
    # Analysis options
    parser.add_argument(
        '--no-statistical-analysis',
        dest='enable_statistical_analysis',
        action='store_false',
        help='Skip statistical significance testing'
    )
    parser.add_argument(
        '--significance-level',
        type=float,
        default=0.05,
        help='Statistical significance threshold (default: 0.05)'
    )
    parser.add_argument(
        '--multiple-comparison-correction',
        choices=['bonferroni', 'holm', 'fdr'],
        default='fdr',
        help='Multiple comparison correction method (default: fdr)'
    )
    
    # Visualization parameters
    parser.add_argument(
        '--plot-format',
        choices=['png', 'pdf', 'svg'],
        default='png',
        help='Output format for plots (default: png)'
    )
    parser.add_argument(
        '--plot-dpi',
        type=int,
        default=300,
        help='Plot resolution in DPI (default: 300)'
    )
    
    # Output control
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='Suppress progress messages'
    )
    parser.add_argument(
        '-l', '--log',
        action='store_true',
        help='Enable detailed logging'
    )


def main(args):
    """Main function for the 'compare' subcommand."""
    setup_logging(args.log, args.quiet)
    
    # Validate inputs
    validate_inputs(args)
    
    # Parse oligo set specifications
    oligo_sets = parse_oligo_sets(args.sets)
    
    logging.info(f"Starting comparative analysis of {len(oligo_sets)} oligo sets...")
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize comparative analyzer
    analyzer = ComparativeAnalyzer(
        reference_file=args.reference,
        output_dir=output_dir,
        min_identity=args.min_identity,
        min_length=args.min_length,
        min_coverage=args.min_coverage,
        target_coverage=args.target_coverage,
        threads=args.threads,
        mapper=args.mapper,
        minimap2_preset=args.minimap2_preset,
        strand=args.strand
    )
    
    # Add each oligo set for analysis
    for name, oligo_file in oligo_sets.items():
        logging.info(f"Adding oligo set '{name}' for analysis...")
        analyzer.add_oligo_set(name, oligo_file)
    
    # Generate comparison matrix
    logging.info("Generating comparison matrix...")
    comparison_matrix = analyzer.generate_comparison_matrix()
    
    # Initialize differential analyzer
    differential_analyzer = None
    if args.enable_statistical_analysis:
        logging.info("Initializing statistical analysis...")
        differential_analyzer = DifferentialAnalyzer(
            significance_level=args.significance_level,
            correction_method=args.multiple_comparison_correction
        )
    
    # Generate visualizations
    logging.info("Generating comparative visualizations...")
    visualizer = ComparativeVisualizer(
        output_dir=output_dir,
        plot_format=args.plot_format,
        dpi=args.plot_dpi
    )
    
    plots = visualizer.generate_all_comparative_plots(analyzer, differential_analyzer)
    
    # Generate comprehensive report
    logging.info("Generating comprehensive comparison report...")
    report_generator = ComparativeReportGenerator(
        analyzer=analyzer,
        output_dir=output_dir,
        differential_analyzer=differential_analyzer
    )
    
    report_file = report_generator.generate_report()
    
    # Export comparison data
    logging.info("Exporting comparison data...")
    exported_files = analyzer.export_comparison_data()
    
    statistics = None
    if differential_analyzer is not None:
        statistics = {'per_reference_metrics': differential_analyzer.compare_quality_metrics(analyzer.oligo_sets)}
        if len(analyzer.oligo_sets) == 2:
            first, second = analyzer.oligo_sets
            statistics['coverage_distribution'] = differential_analyzer.compare_coverage_distributions(first, second)
            statistics['bait_identity'] = differential_analyzer.compare_oligo_identity(first, second)
            statistics['gap_sizes'] = differential_analyzer.analyze_gap_patterns(first, second)
    json_file = write_json({
        'versions': tool_versions(),
        'arguments': arguments_record(args),
        'reference': str(args.reference),
        'parameters': {
            'min_identity': args.min_identity, 'min_length': args.min_length,
            'min_coverage': args.min_coverage, 'target_coverage': args.target_coverage,
            'significance_level': args.significance_level,
            'multiple_comparison_correction': args.multiple_comparison_correction,
        },
        'sets': [{
            'name': r.name, 'file': r.file_path,
            'coverage_stats': r.coverage_stats,
            'gap_analysis': {k: v for k, v in r.gap_analysis.items() if k != 'feature_analysis'},
            'quality_score': r.quality_score.to_dict(),
            'benchmarks': r.benchmark_results,
        } for r in analyzer.oligo_sets],
        'comparison_matrix': comparison_matrix,
        'ranking': analyzer.generate_ranking(),
        'statistics': statistics,
    }, output_dir / "comparison.json")
    exported_files['json'] = str(json_file)
    
    # Print summary
    print_summary(analyzer, comparison_matrix, plots, report_file, exported_files)
    
    logging.info(f"Comparative analysis complete! Results saved to {output_dir}")


def setup_logging(enable_debug: bool, quiet: bool) -> None:
    """Configure logging (kept for direct calls of main)."""
    ensure_logging(verbose=enable_debug, quiet=quiet)

def validate_inputs(args) -> None:
    """Validate input files and parameters."""
    # Check reference file exists
    if not Path(args.reference).exists():
        logging.error(f"Reference file not found: {args.reference}")
        sys.exit(1)
    
    # Check oligo set specifications
    if len(args.sets) < 2:
        logging.error("Need at least 2 oligo sets for comparison")
        sys.exit(1)
    
    # Validate oligo set format and files
    for set_spec in args.sets:
        if ':' not in set_spec:
            logging.error(f"Invalid oligo set specification: {set_spec}")
            logging.error("Use format 'Name:path/to/oligos.fasta'")
            sys.exit(1)
        
        name, file_path = set_spec.split(':', 1)
        if not Path(file_path).exists():
            logging.error(f"Oligo file not found: {file_path}")
            sys.exit(1)
    
    # Check parameter ranges
    if args.min_identity < 0 or args.min_identity > 100:
        logging.error("Minimum identity must be between 0 and 100")
        sys.exit(1)
    
    if args.min_coverage < 0:
        logging.error("Minimum coverage must be >= 0")
        sys.exit(1)
    
    if args.target_coverage < args.min_coverage:
        logging.error("Target coverage must be >= minimum coverage")
        sys.exit(1)
    
    if args.significance_level <= 0 or args.significance_level >= 1:
        logging.error("Significance level must be between 0 and 1")
        sys.exit(1)
    
    if not check_mapper_available(args.mapper):
        logging.error(f"{args.mapper} not found in PATH. Please install {args.mapper}.")
        sys.exit(1)


def parse_oligo_sets(set_specs: List[str]) -> Dict[str, str]:
    """Parse oligo set specifications."""
    oligo_sets = {}
    
    for spec in set_specs:
        name, file_path = spec.split(':', 1)
        
        # Validate name uniqueness
        if name in oligo_sets:
            logging.error(f"Duplicate oligo set name: {name}")
            sys.exit(1)
        
        oligo_sets[name] = file_path
    
    return oligo_sets


def print_summary(analyzer: ComparativeAnalyzer, comparison_matrix, plots: Dict[str, str],
                 report_file: str, exported_files: Dict[str, str]) -> None:
    """Print comparative analysis summary."""
    
    print("\n" + "="*80)
    print("COMPARATIVE ANALYSIS SUMMARY")
    print("="*80)
    
    # Best performer
    best_performer = analyzer.identify_best_performer('quality_score')
    ranking = analyzer.generate_ranking()
    
    print(f"Oligo Sets Compared:    {len(analyzer.oligo_sets)}")
    print(f"Best Performer:         {best_performer.name}")
    print(f"Top Quality Score:      {best_performer.quality_score.overall_score:.2f} (0-1, {best_performer.quality_score.category.value})")
    
    # Coverage range
    coverage_values = [result.coverage_stats.get('coverage_breadth', 0) for result in analyzer.oligo_sets]
    print(f"Coverage Range:         {min(coverage_values):.1f}% - {max(coverage_values):.1f}%")
    
    # Gap range
    gap_values = [result.gap_analysis.get('total_gaps', 0) for result in analyzer.oligo_sets]
    print(f"Gap Count Range:        {min(gap_values):,} - {max(gap_values):,}")
    
    print("\nRanking (Top 3):")
    for i, (name, score) in enumerate(ranking[:3], 1):
        oligo_set = next(result for result in analyzer.oligo_sets if result.name == name)
        print(f"  {i}. {name:<20} Score: {score:.2f}, Category: {oligo_set.quality_score.category.value}")
    
    print("\nResults Summary:")
    print(f"  📊 Interactive Report:   {Path(report_file).name}")
    print(f"  📈 Visualizations:       {len(plots)} plots generated")
    print(f"  📋 Data Exports:         {len(exported_files)} files exported")
    
    print("\nKey Files:")
    print(f"  • {Path(report_file).name}")
    for file_type, file_path in exported_files.items():
        print(f"  • {Path(file_path).name}")
    
    print("="*80)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Comparative analysis of multiple oligo sets'
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)