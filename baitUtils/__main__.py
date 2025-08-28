# baitUtils/__main__.py

import argparse
import sys
from baitUtils._version import __version__

# Import refactored subcommands
from baitUtils.sequence_statistics import add_arguments as stats_add_args, main as stats_main
from baitUtils.statistical_plots import add_arguments as plot_add_args, main as plot_main
from baitUtils.sequence_mapping import add_arguments as map_add_args, main as map_main
from baitUtils.coverage_evaluation import add_arguments as check_add_args, main as check_main
from baitUtils.gap_filling import add_arguments as fill_add_args, main as fill_main
from baitUtils.evaluate import add_arguments as evaluate_add_args, main as evaluate_main
from baitUtils.compare import add_arguments as compare_add_args, main as compare_main


def main():
    parser = argparse.ArgumentParser(prog="baitUtils")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command")

    # ---------------------
    # stats subcommand
    # ---------------------
    parser_stats = subparsers.add_parser("stats", help="Calculate sequence statistics with filtering")
    stats_add_args(parser_stats)
    parser_stats.set_defaults(func=stats_main)

    # ---------------------
    # plot subcommand
    # ---------------------
    parser_plot = subparsers.add_parser("plot", help="Generate statistical plots from sequence data")
    plot_add_args(parser_plot)
    parser_plot.set_defaults(func=plot_main)

    # ---------------------
    # map subcommand
    # ---------------------
    parser_map = subparsers.add_parser("map", help="Map sequences against reference genomes")
    map_add_args(parser_map)
    parser_map.set_defaults(func=map_main)

    # ---------------------
    # check subcommand
    # ---------------------
    parser_check = subparsers.add_parser("check", help="Evaluate oligo coverage and report gaps")
    check_add_args(parser_check)
    parser_check.set_defaults(func=check_main)

    # ---------------------
    # fill subcommand
    # ---------------------
    parser_fill = subparsers.add_parser("fill", help="Multi-pass gap filling to maximize coverage")
    fill_add_args(parser_fill)
    parser_fill.set_defaults(func=fill_main)

    # ---------------------
    # evaluate subcommand
    # ---------------------
    parser_evaluate = subparsers.add_parser("evaluate", help="Comprehensive oligo set coverage evaluation")
    evaluate_add_args(parser_evaluate)
    parser_evaluate.set_defaults(func=evaluate_main)

    # ---------------------
    # compare subcommand
    # ---------------------
    parser_compare = subparsers.add_parser("compare", help="Comparative analysis of multiple oligo sets")
    compare_add_args(parser_compare)
    parser_compare.set_defaults(func=compare_main)

    # Parse the command line and call the appropriate subcommand
    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()