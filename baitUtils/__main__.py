# baitUtils/__main__.py

import argparse
import sys
from pathlib import Path

from baitUtils.config import load_config, apply_config
from baitUtils.logging_utils import configure_logging, LEVELS
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
    parser.add_argument("--config", type=Path, default=None,
                        help="JSON file with option defaults per command (see baitUtils.config)")
    parser.add_argument("--log-level", dest="log_level_global", choices=LEVELS, default=None,
                        help="Logging level for any command (default: INFO)")
    subparsers = parser.add_subparsers(dest="command")
    commands = {}

    # ---------------------
    # stats subcommand
    # ---------------------
    parser_stats = subparsers.add_parser("stats", help="Calculate sequence statistics with filtering")
    stats_add_args(parser_stats)
    parser_stats.set_defaults(func=stats_main)
    commands["stats"] = parser_stats

    # ---------------------
    # plot subcommand
    # ---------------------
    parser_plot = subparsers.add_parser("plot", help="Generate statistical plots from sequence data")
    plot_add_args(parser_plot)
    parser_plot.set_defaults(func=plot_main)
    commands["plot"] = parser_plot

    # ---------------------
    # map subcommand
    # ---------------------
    parser_map = subparsers.add_parser("map", help="Map sequences against reference genomes")
    map_add_args(parser_map)
    parser_map.set_defaults(func=map_main)
    commands["map"] = parser_map

    # ---------------------
    # check subcommand
    # ---------------------
    parser_check = subparsers.add_parser("check", help="Evaluate oligo coverage and report gaps")
    check_add_args(parser_check)
    parser_check.set_defaults(func=check_main)
    commands["check"] = parser_check

    # ---------------------
    # fill subcommand
    # ---------------------
    parser_fill = subparsers.add_parser("fill", help="Multi-pass gap filling to maximize coverage")
    fill_add_args(parser_fill)
    parser_fill.set_defaults(func=fill_main)
    commands["fill"] = parser_fill

    # ---------------------
    # evaluate subcommand
    # ---------------------
    parser_evaluate = subparsers.add_parser("evaluate", help="Comprehensive oligo set coverage evaluation")
    evaluate_add_args(parser_evaluate)
    parser_evaluate.set_defaults(func=evaluate_main)
    commands["evaluate"] = parser_evaluate

    # ---------------------
    # compare subcommand
    # ---------------------
    parser_compare = subparsers.add_parser("compare", help="Comparative analysis of multiple oligo sets")
    compare_add_args(parser_compare)
    parser_compare.set_defaults(func=compare_main)
    commands["compare"] = parser_compare

    # Apply the configuration file before parsing the full command line
    argv = sys.argv[1:]
    if "--config" in argv:
        config_path = Path(argv[argv.index("--config") + 1])
        apply_config(commands, load_config(config_path))

    args = parser.parse_args(argv)
    configure_logging(
        level=args.log_level_global or getattr(args, "log_level", None) or "INFO",
        verbose=bool(getattr(args, "log", False) or getattr(args, "verbose", False)),
        quiet=bool(getattr(args, "quiet", False)),
    )
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()