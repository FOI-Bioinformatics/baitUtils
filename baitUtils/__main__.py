# baitUtils/__main__.py

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from baitUtils.config import load_config, apply_config
from baitUtils.cli_utils import logging_kwargs
from baitUtils.logging_utils import configure_logging, LEVELS
from baitUtils._version import __version__

from baitUtils.sequence_statistics import add_arguments as stats_add_args, main as stats_main
from baitUtils.statistical_plots import add_arguments as plot_add_args, main as plot_main
from baitUtils.sequence_mapping import add_arguments as map_add_args, main as map_main
from baitUtils.coverage_evaluation import add_arguments as check_add_args, main as check_main
from baitUtils.gap_filling import add_arguments as fill_add_args, main as fill_main
from baitUtils.evaluate import add_arguments as evaluate_add_args, main as evaluate_main
from baitUtils.compare import add_arguments as compare_add_args, main as compare_main


COMMANDS = (
    ("stats", "Calculate sequence statistics with filtering", stats_add_args, stats_main),
    ("plot", "Generate statistical plots from sequence data", plot_add_args, plot_main),
    ("map", "Map sequences against reference genomes", map_add_args, map_main),
    ("check", "Evaluate oligo coverage and report gaps", check_add_args, check_main),
    ("fill", "Multi-pass gap filling to maximize coverage", fill_add_args, fill_main),
    ("evaluate", "Comprehensive oligo set coverage evaluation", evaluate_add_args, evaluate_main),
    ("compare", "Comparative analysis of multiple oligo sets", compare_add_args, compare_main),
)


def build_parser() -> Tuple[argparse.ArgumentParser, Dict[str, argparse.ArgumentParser]]:
    """The main parser with one subparser per command."""
    parser = argparse.ArgumentParser(prog="baitUtils")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--config", type=Path, default=None,
                        help="JSON file with option defaults per command (see baitUtils.config)")
    parser.add_argument("--log-level", dest="log_level_global", choices=LEVELS, default=None,
                        help="Logging level for any command (default: INFO)")
    subparsers = parser.add_subparsers(dest="command")
    commands = {}
    for name, help_text, add_args, run in COMMANDS:
        sub = subparsers.add_parser(name, help=help_text)
        add_args(sub)
        sub.set_defaults(func=run)
        commands[name] = sub
    return parser, commands


def _extract_config(argv: List[str]) -> Tuple[Optional[Path], List[str]]:
    """
    Find --config FILE (or --config=FILE) anywhere on the command line and
    return the path and the remaining arguments. A --config without a file
    is a usage error.
    """
    pre = argparse.ArgumentParser(prog="baitUtils", add_help=False, allow_abbrev=False)
    pre.add_argument("--config", type=Path, default=None)
    known, rest = pre.parse_known_args(argv)
    return known.config, rest


def main(argv: Optional[List[str]] = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser, commands = build_parser()

    # Apply the configuration file before parsing the full command line so
    # that its values act as defaults and command-line values win
    config_path, argv = _extract_config(argv)
    if config_path is not None:
        apply_config(commands, load_config(config_path))

    args = parser.parse_args(argv)
    args.config = config_path
    configure_logging(**logging_kwargs(args))
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
