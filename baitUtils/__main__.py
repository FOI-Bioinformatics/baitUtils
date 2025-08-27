# baitUtils/__main__.py

import argparse
import sys
from baitUtils._version import __version__  # Import version from _version.py
from baitUtils import stats, plot, map

# Import your newly created subcommands
from baitUtils.check import add_arguments as check_add_args, main as check_main
from baitUtils.fill import add_arguments as fill_add_args, main as fill_main
from baitUtils.evaluate import add_arguments as evaluate_add_args, main as evaluate_main
from baitUtils.compare import add_arguments as compare_add_args, main as compare_main


def main():
    parser = argparse.ArgumentParser(prog="baitUtils")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command")

    # ---------------------
    # stats subcommand
    # ---------------------
    parser_stats = subparsers.add_parser("stats", help="Run bait statistics calculations")
    stats.add_arguments(parser_stats)
    parser_stats.set_defaults(func=stats.main)

    # ---------------------
    # plot subcommand
    # ---------------------
    parser_plot = subparsers.add_parser("plot", help="Generate plots from bait statistics")
    plot.add_arguments(parser_plot)
    parser_plot.set_defaults(func=plot.main)

    # ---------------------
    # map subcommand
    # ---------------------
    parser_map = subparsers.add_parser("map", help="Map baits against target and off-target sequences")
    map.add_arguments(parser_map)
    parser_map.set_defaults(func=map.main)

    # ---------------------
    # check subcommand
    # ---------------------
    parser_check = subparsers.add_parser("check", help="Evaluate coverage of oligos")
    check_add_args(parser_check)
    parser_check.set_defaults(func=check_main)

    # ---------------------
    # fill subcommand
    # ---------------------
    parser_fill = subparsers.add_parser("fill", help="Multi-pass selection to fill coverage gaps")
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