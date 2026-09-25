#!/usr/bin/env python3

"""
evaluate.py

Comprehensive oligo set coverage evaluation command.
Integrates mapping, coverage analysis, visualization, and reporting
to provide complete assessment of oligo set performance against reference sequences.

Usage:
    baitUtils evaluate -i oligos.fasta -r reference.fasta -o coverage_report/
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List
import tempfile
import shutil

from baitUtils._version import __version__
from baitUtils.logging_utils import ensure_logging
from baitUtils.coverage_stats import CoverageAnalyzer
from baitUtils.coverage_viz import CoverageVisualizer
from baitUtils.json_export import write_json, tool_versions, arguments_record
from baitUtils.mapping_utils import run_mapper, check_mapper_available
from baitUtils.gap_analysis import GapAnalyzer
from baitUtils.reference_analyzer import ReferenceAnalyzer
from baitUtils.quality_scorer import QualityScorer
from baitUtils.report_generator import InteractiveReportGenerator
from baitUtils.interactive_plots import InteractivePlotter
from baitUtils.benchmark import BenchmarkAnalyzer


def add_arguments(parser):
    """Add command-line arguments for the 'evaluate' subcommand."""
    parser.add_argument('--version', action='version', version=f"baitUtils evaluate {__version__}")
    
    # Required arguments
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input oligos/baits FASTA file'
    )
    parser.add_argument(
        '-r', '--reference',
        required=True,
        help='Reference genome/target FASTA file'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory for coverage evaluation results'
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
    
    # Gap analysis parameters
    parser.add_argument(
        '--min-gap-size',
        type=int,
        default=100,
        help='Minimum gap size to report (default: 100)'
    )
    parser.add_argument(
        '--gap-extend',
        type=int,
        default=0,
        help='Extend gaps by N bases on each side for analysis (default: 0)'
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
    
    # Phase 2 features
    parser.add_argument(
        '--no-html-report',
        dest='enable_html_report',
        action='store_false',
        help='Skip the interactive HTML report'
    )
    parser.add_argument(
        '--no-interactive-plots',
        dest='enable_interactive_plots',
        action='store_false',
        help='Skip interactive plots'
    )
    parser.add_argument(
        '--offline-plots',
        action='store_true',
        help='Embed plotly.js in the HTML report instead of loading it from a CDN'
    )
    parser.add_argument(
        '--reference-analysis-window',
        type=int,
        default=1000,
        help='Window size for reference sequence analysis (default: 1000)'
    )
    parser.add_argument(
        '--no-benchmarking',
        dest='enable_benchmarking',
        action='store_false',
        help='Skip benchmarking against design targets'
    )
    
    # Output control
    parser.add_argument(
        '--keep-intermediates',
        action='store_true',
        help='Keep intermediate mapping files'
    )
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
    """Main function for the 'evaluate' subcommand."""
    setup_logging(args.log, args.quiet)
    
    # Validate inputs
    validate_inputs(args)
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        logging.info("Starting comprehensive oligo set evaluation...")
        
        # Step 1: Perform mapping
        logging.info("Step 1/5: Mapping oligos to reference sequences...")
        psl_file = perform_mapping(args, temp_path)
        
        # Step 2: Analyze coverage statistics
        logging.info("Step 2/5: Analyzing coverage statistics...")
        analyzer = CoverageAnalyzer(
            psl_file=psl_file,
            reference_file=args.reference,
            min_coverage=args.min_coverage,
            target_coverage=args.target_coverage,
            min_identity=args.min_identity,
            min_length=args.min_length,
            oligos_file=args.input,
            strand=args.strand
        )
        
        coverage_stats = analyzer.analyze()
        
        # Step 3: Perform gap analysis
        logging.info("Step 3/5: Performing gap analysis...")
        gap_analyzer = GapAnalyzer(
            coverage_data=coverage_stats,
            reference_file=args.reference,
            min_gap_size=args.min_gap_size,
            extend_bp=args.gap_extend,
            coverage_arrays=analyzer.coverage_arrays,
            min_coverage=args.min_coverage
        )
        
        gap_analysis = gap_analyzer.analyze()
        
        # Step 4: Generate visualizations
        logging.info("Step 4/5: Generating visualizations...")
        visualizer = CoverageVisualizer(
            coverage_stats=coverage_stats,
            gap_analysis=gap_analysis,
            output_dir=output_dir,
            plot_format=args.plot_format,
            dpi=args.plot_dpi
        )
        
        visualizer.generate_all_plots()
        
        # Phase 2 features - Extended analysis and reporting
        reference_analysis = None
        quality_score = None
        benchmark_results = None
        
        if args.enable_html_report or args.enable_interactive_plots or args.enable_benchmarking:
            # Step 5: Analyze reference sequences
            logging.info("Step 5/8: Analyzing reference sequences...")
            ref_analyzer = ReferenceAnalyzer(
                reference_file=args.reference,
                coverage_data=coverage_stats,
                window_size=args.reference_analysis_window,
                coverage_arrays=analyzer.coverage_arrays,
                min_coverage=args.min_coverage
            )
            reference_analysis = ref_analyzer.analyze()

            # Step 6: Calculate quality scores
            logging.info("Step 6/8: Calculating quality scores...")
            quality_scorer = QualityScorer(
                coverage_stats=coverage_stats,
                gap_analysis=gap_analysis,
                reference_analysis=reference_analysis
            )
            quality_score = quality_scorer.calculate_quality_score()
            
            # Step 7: Run benchmarking analysis
            if args.enable_benchmarking:
                logging.info("Step 7/9: Running benchmark analysis...")
                benchmark_analyzer = BenchmarkAnalyzer(
                    coverage_stats=coverage_stats,
                    gap_analysis=gap_analysis,
                    reference_analysis=reference_analysis,
                    quality_score=quality_score
                )
                benchmark_results = benchmark_analyzer.run_full_benchmark()
            
            # Step 8: Generate interactive plots
            if args.enable_interactive_plots:
                logging.info("Step 8/9: Generating interactive plots...")
                interactive_plotter = InteractivePlotter(
                    coverage_stats=coverage_stats,
                    gap_analysis=gap_analysis,
                    reference_analysis=reference_analysis,
                    quality_scores=quality_score.to_dict(),
                    output_dir=output_dir
                )
                interactive_plotter.create_all_interactive_plots()
            
            # Step 9: Generate interactive HTML report
            if args.enable_html_report:
                logging.info("Step 9/9: Generating interactive HTML report...")
                report_generator = InteractiveReportGenerator(
                    coverage_stats=coverage_stats,
                    gap_analysis=gap_analysis,
                    reference_analysis=reference_analysis,
                    quality_scores=quality_score.to_dict(),
                    output_dir=output_dir,
                    offline_plots=args.offline_plots
                )
                report_generator.generate_report()
        
        # Step 5/Final: Generate standard reports
        step_num = "5/5" if not (args.enable_html_report or args.enable_interactive_plots or args.enable_benchmarking) else "Final"
        logging.info(f"Step {step_num}: Generating standard reports...")
        generate_reports(coverage_stats, gap_analysis, output_dir, args, benchmark_results,
                         quality_score=quality_score, reference_analysis=reference_analysis)
        
        # Copy intermediate files if requested
        if args.keep_intermediates:
            intermediates_dir = output_dir / "intermediates"
            intermediates_dir.mkdir(exist_ok=True)
            shutil.copy2(psl_file, intermediates_dir / "mapping.psl")
            logging.info(f"Intermediate files saved to {intermediates_dir}")
    
    logging.info(f"Coverage evaluation complete! Results saved to {output_dir}")
    
    # Print enhanced summary to stdout
    print_summary(coverage_stats, gap_analysis, quality_score)


def setup_logging(enable_debug: bool, quiet: bool) -> None:
    """Configure logging (kept for direct calls of main)."""
    ensure_logging(verbose=enable_debug, quiet=quiet)

def validate_inputs(args) -> None:
    """Validate input files and parameters."""
    # Check input files exist
    if not Path(args.input).exists():
        logging.error(f"Input oligos file not found: {args.input}")
        sys.exit(1)
    
    if not Path(args.reference).exists():
        logging.error(f"Reference file not found: {args.reference}")
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
    
    if not check_mapper_available(args.mapper):
        logging.error(f"{args.mapper} not found in PATH. Please install {args.mapper}.")
        sys.exit(1)


def perform_mapping(args, temp_dir: Path) -> Path:
    """Map oligos with the chosen mapper; returns the PSL or PAF path."""
    try:
        psl_file = Path(run_mapper(args.mapper, args.reference, args.input, str(temp_dir), "mapping",
                                   threads=args.threads, min_identity=args.min_identity,
                                   minimap2_preset=args.minimap2_preset))
    except RuntimeError as e:
        logging.error(str(e))
        sys.exit(1)
    
    if not psl_file.exists() or psl_file.stat().st_size == 0:
        logging.error(f"{args.mapper} produced no output")
        sys.exit(1)
    
    logging.info(f"Mapping completed. Results saved to {psl_file}")
    return psl_file


def generate_reports(coverage_stats: Dict, gap_analysis: Dict, output_dir: Path, args,
                     benchmark_results=None, quality_score=None, reference_analysis=None) -> None:
    """Generate text reports and the machine-readable evaluation.json."""
    
    write_json({
        'versions': tool_versions(),
        'arguments': arguments_record(args),
        'inputs': {'oligos': str(args.input), 'reference': str(args.reference)},
        'parameters': {
            'min_identity': args.min_identity, 'min_length': args.min_length,
            'min_coverage': args.min_coverage, 'target_coverage': args.target_coverage,
            'min_gap_size': args.min_gap_size,
        },
        'coverage_stats': coverage_stats,
        'gap_analysis': {k: v for k, v in gap_analysis.items() if k != 'feature_analysis'},
        'quality_score': quality_score.to_dict() if quality_score is not None else None,
        'benchmarks': benchmark_results,
        'reference_summary': (reference_analysis or {}).get('summary'),
    }, output_dir / "evaluation.json")
    
    # Coverage statistics report
    stats_file = output_dir / "coverage_statistics.txt"
    with open(stats_file, 'w') as f:
        f.write("=== OLIGO SET COVERAGE EVALUATION REPORT ===\n\n")
        f.write(f"Input oligos: {args.input}\n")
        f.write(f"Reference: {args.reference}\n")
        f.write(f"Analysis date: {coverage_stats.get('analysis_date', 'Unknown')}\n\n")
        
        f.write("=== OVERALL COVERAGE STATISTICS ===\n")
        f.write(f"Total oligos: {coverage_stats.get('total_oligos', 0):,}\n")
        f.write(f"Mapped oligos: {coverage_stats.get('mapped_oligos', 0):,}\n")
        f.write(f"Mapping efficiency: {coverage_stats.get('mapping_efficiency', 0):.1f}%\n\n")
        
        f.write(f"Reference length: {coverage_stats.get('reference_length', 0):,} bp\n")
        f.write(f"N bases: {coverage_stats.get('n_bases', 0):,} bp\n")
        f.write(f"Covered bases: {coverage_stats.get('covered_bases', 0):,} bp\n")
        f.write(f"Coverage breadth: {coverage_stats.get('coverage_breadth', 0):.1f}%\n")
        f.write(f"Coverage breadth excluding N: {coverage_stats.get('coverage_breadth_non_n', 0):.1f}%\n\n")
        
        f.write(f"Mean coverage depth: {coverage_stats.get('mean_depth', 0):.1f}x\n")
        f.write(f"Median coverage depth: {coverage_stats.get('median_depth', 0):.1f}x\n")
        f.write(f"Coverage std dev: {coverage_stats.get('depth_std', 0):.1f}x\n\n")
        
        f.write("=== COVERAGE DEPTH DISTRIBUTION ===\n")
        depth_dist = coverage_stats.get('depth_distribution', {})
        for threshold, percent in depth_dist.items():
            f.write(f"Bases with >={threshold}x coverage: {percent:.1f}%\n")
        f.write("\n")
    
    # Gap analysis report
    gap_file = output_dir / "gap_analysis.txt"
    with open(gap_file, 'w') as f:
        f.write("=== GAP ANALYSIS REPORT ===\n\n")
        
        f.write(f"Total gaps: {gap_analysis.get('total_gaps', 0):,}\n")
        f.write(f"Total gap length: {gap_analysis.get('total_gap_length', 0):,} bp\n")
        f.write(f"Gap percentage: {gap_analysis.get('gap_percentage', 0):.1f}%\n\n")
        
        f.write(f"Mean gap size: {gap_analysis.get('mean_gap_size', 0):.1f} bp\n")
        f.write(f"Median gap size: {gap_analysis.get('median_gap_size', 0):.1f} bp\n")
        f.write(f"Largest gap: {gap_analysis.get('max_gap_size', 0):,} bp\n\n")
        
        f.write("=== GAP SIZE DISTRIBUTION ===\n")
        size_dist = gap_analysis.get('size_distribution', {})
        for size_range, count in size_dist.items():
            f.write(f"{size_range}: {count} gaps\n")
        f.write("\n")
        
        f.write("=== LARGEST GAPS ===\n")
        largest_gaps = gap_analysis.get('largest_gaps', [])
        for i, gap in enumerate(largest_gaps[:10], 1):
            f.write(f"{i}. {gap['chromosome']}:{gap['start']}-{gap['end']} "
                   f"({gap['length']:,} bp)\n")
    
    # Recommendations report
    rec_file = output_dir / "recommendations.txt"
    with open(rec_file, 'w') as f:
        f.write("=== COVERAGE IMPROVEMENT RECOMMENDATIONS ===\n\n")
        
        # Generate recommendations based on analysis
        recommendations = generate_recommendations(coverage_stats, gap_analysis)
        
        for i, rec in enumerate(recommendations, 1):
            f.write(f"{i}. {rec}\n\n")
    
    # Benchmarking report
    if benchmark_results is not None:
        benchmark_file = output_dir / "benchmark_analysis.txt"
        with open(benchmark_file, 'w') as f:
            # Import here to avoid circular imports
            from baitUtils.benchmark import BenchmarkAnalyzer
            dummy_analyzer = BenchmarkAnalyzer({}, {}, {}, None)
            report_content = dummy_analyzer.generate_benchmark_report(benchmark_results)
            f.write(report_content)
    
    logging.info("Reports generated successfully")


def generate_recommendations(coverage_stats: Dict, gap_analysis: Dict) -> List[str]:
    """Generate improvement recommendations based on analysis."""
    recommendations = []
    
    # Coverage breadth recommendations
    breadth = coverage_stats.get('coverage_breadth', 0)
    if breadth < 80:
        recommendations.append(
            f"Coverage breadth is only {breadth:.1f}%. Consider adding more oligos "
            "to target uncovered regions, particularly the largest gaps."
        )
    
    # Coverage depth recommendations  
    mean_depth = coverage_stats.get('mean_depth', 0)
    depth_std = coverage_stats.get('depth_std', 0)
    if depth_std / mean_depth > 1.0 if mean_depth > 0 else True:
        recommendations.append(
            "Coverage depth is highly variable (high coefficient of variation). "
            "Consider redistributing oligos for more uniform coverage."
        )
    
    # Gap-based recommendations
    total_gaps = gap_analysis.get('total_gaps', 0)
    if total_gaps > 100:
        recommendations.append(
            f"Found {total_gaps} coverage gaps. Use 'baitUtils fill' command "
            "to iteratively select additional oligos for gap closure."
        )
    
    # Large gap recommendations
    max_gap = gap_analysis.get('max_gap_size', 0)
    if max_gap > 10000:
        recommendations.append(
            f"Largest gap is {max_gap:,} bp. Very large gaps may indicate "
            "problematic reference regions requiring specialized oligo design."
        )
    
    # Mapping efficiency recommendations
    mapping_eff = coverage_stats.get('mapping_efficiency', 0)
    if mapping_eff < 80:
        recommendations.append(
            f"Only {mapping_eff:.1f}% of oligos mapped successfully. "
            "Consider reviewing oligo design parameters or reference sequence quality."
        )
    
    if not recommendations:
        recommendations.append(
            "Coverage quality appears good! No major issues identified."
        )
    
    return recommendations


def print_summary(coverage_stats: Dict, gap_analysis: Dict, quality_score=None) -> None:
    """Print a brief summary to stdout."""
    print("\n" + "="*60)
    print("COVERAGE EVALUATION SUMMARY")
    print("="*60)
    
    print(f"Coverage Breadth:     {coverage_stats.get('coverage_breadth', 0):6.1f}%")
    print(f"Mean Coverage Depth:  {coverage_stats.get('mean_depth', 0):6.1f}x")
    print(f"Total Gaps:           {gap_analysis.get('total_gaps', 0):6,}")
    print(f"Largest Gap:          {gap_analysis.get('max_gap_size', 0):6,} bp")
    print(f"Mapping Efficiency:   {coverage_stats.get('mapping_efficiency', 0):6.1f}%")
    
    if quality_score is not None:
        print(f"Overall Quality:      {quality_score.overall_score:6.2f} (0-1)")
        print(f"Quality Category:     {quality_score.category.value:>6s}")
    
    print("="*60)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Comprehensive oligo set coverage evaluation'
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)