"""
config.py

Optional JSON configuration file for command-line defaults.

The file maps command names to option values, using the option's
destination name (for example "min_identity" for --min-identity). Values
under the key "common" apply to every command that has that option.
Options given on the command line override the file.

Example:
    {
      "common": {"threads": 4, "log_level": "DEBUG"},
      "evaluate": {"min_identity": 95.0, "target_coverage": 10},
      "compare": {"multiple_comparison_correction": "holm"}
    }
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, Any

COMMON = "common"


def load_config(path: Path) -> Dict[str, Any]:
    """Read the JSON configuration file."""
    with open(path) as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError("Configuration file must contain a JSON object")
    return data


def apply_config(subparsers: Dict[str, argparse.ArgumentParser], config: Dict[str, Any]) -> None:
    """
    Set defaults on each subparser from the configuration. An option that
    is required on the command line becomes optional when the file sets it.

    Warnings are given for a section that is not a command, for a key in a
    command section that the command does not have, and for a key in the
    common section that no command has. A common key is applied silently
    to the commands that have it.
    """
    dests_per_command = {
        command: {action.dest: action for action in parser._actions}
        for command, parser in subparsers.items()
    }

    unknown_sections = [key for key in config if key != COMMON and key not in subparsers]
    if unknown_sections:
        logging.warning(f"Configuration sections that are not commands: {', '.join(unknown_sections)}")

    common = config.get(COMMON, {}) or {}
    all_dests = set().union(*(set(d) for d in dests_per_command.values())) if dests_per_command else set()
    unused_common = [key for key in common if key not in all_dests]
    if unused_common:
        logging.warning(f"Configuration keys in '{COMMON}' used by no command: {', '.join(unused_common)}")

    for command, parser in subparsers.items():
        dests = dests_per_command[command]
        values = {key: value for key, value in common.items() if key in dests}
        specific = config.get(command, {}) or {}
        unknown = [key for key in specific if key not in dests]
        if unknown:
            logging.warning(f"Configuration keys not used by {command}: {', '.join(unknown)}")
        values.update({key: value for key, value in specific.items() if key in dests})
        if not values:
            continue
        for key in values:
            dests[key].required = False
        parser.set_defaults(**values)
