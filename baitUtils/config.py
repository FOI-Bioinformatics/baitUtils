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
    """
    common = config.get("common", {}) or {}
    for command, parser in subparsers.items():
        values = dict(common)
        values.update(config.get(command, {}) or {})
        if not values:
            continue
        dests = {action.dest: action for action in parser._actions}
        unknown = [key for key in values if key not in dests]
        if unknown:
            logging.warning(f"Configuration keys not used by {command}: {', '.join(unknown)}")
        applied = {key: value for key, value in values.items() if key in dests}
        for key in applied:
            dests[key].required = False
        parser.set_defaults(**applied)
