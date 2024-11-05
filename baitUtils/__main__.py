# baitUtils/__main__.py

import argparse
import sys
from baitUtils._version import __version__  # Import version from _version.py
from baitUtils import stats, plot, map  # Import the new map module


def main():
    parser = argparse.ArgumentParser(prog="baitUtils")
    parser.add_argument('--version', action='version', version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command")

    # Subcommand for stats
    parser_stats = subparsers.add_parser("stats", help="Run bait statistics calculations")
    stats.add_arguments(parser_stats)
    parser_stats.set_defaults(func=stats.main)

    # Subcommand for plot
    parser_plot = subparsers.add_parser("plot", help="Generate plots from bait statistics")
    plot.add_arguments(parser_plot)
    parser_plot.set_defaults(func=plot.main)

    # Subcommand for map
    parser_map = subparsers.add_parser("map", help="Map baits against target and off-target sequences")
    map.add_arguments(parser_map)
    parser_map.set_defaults(func=map.main)

    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()