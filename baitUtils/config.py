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
from typing import Dict, Any, List

from baitUtils.cli_utils import InputError

COMMON = "common"


def load_config(path: Path) -> Dict[str, Any]:
    """Read the JSON configuration file."""
    try:
        with open(path) as fh:
            data = json.load(fh)
    except OSError as e:
        raise InputError(f"cannot read configuration file {path}: {e.strerror or e}")
    except json.JSONDecodeError as e:
        raise InputError(f"configuration file {path} is not valid JSON: {e}")
    if not isinstance(data, dict):
        raise InputError("Configuration file must contain a JSON object")
    for section, values in data.items():
        if values is not None and not isinstance(values, dict):
            raise InputError(f"Configuration section '{section}' must be a JSON object")
    return data


def _is_flag(action: argparse.Action) -> bool:
    """store_true / store_false style options, which take no value."""
    return action.nargs == 0 and isinstance(action.const, bool)


def _convert(action: argparse.Action, value: Any, where: str) -> Any:
    """
    Check and convert one configuration value the way argparse treats the
    same value on the command line: the option's type, its choices, a list
    for options that take several values, and a JSON boolean for flags.
    """
    if _is_flag(action):
        if not isinstance(value, bool):
            raise InputError(f"{where}: expected true or false, got {value!r}")
        return value
    if action.nargs in ("+", "*") or isinstance(action.nargs, int):
        if not isinstance(value, list):
            raise InputError(f"{where}: expected a list, got {value!r}")
        if action.nargs == "+" and not value:
            raise InputError(f"{where}: expected at least one value")
        return [_convert_one(action, item, where) for item in value]
    if isinstance(value, (list, dict)):
        raise InputError(f"{where}: expected a single value, got {value!r}")
    return _convert_one(action, value, where)


def _convert_one(action: argparse.Action, value: Any, where: str) -> Any:
    if isinstance(value, bool):
        raise InputError(f"{where}: expected a value, got {value!r}")
    if value is not None and action.type is not None:
        try:
            value = action.type(value if isinstance(value, str) else str(value))
        except (argparse.ArgumentTypeError, TypeError, ValueError) as e:
            raise InputError(f"{where}: invalid value {value!r} ({e})")
    if action.choices is not None and value not in action.choices:
        allowed = ", ".join(repr(choice) for choice in action.choices)
        raise InputError(f"{where}: invalid choice {value!r} (choose from {allowed})")
    return value


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

    errors: List[str] = []
    for command, parser in subparsers.items():
        dests = dests_per_command[command]
        values = {key: (COMMON, value) for key, value in common.items() if key in dests}
        specific = config.get(command, {}) or {}
        unknown = [key for key in specific if key not in dests]
        if unknown:
            logging.warning(f"Configuration keys not used by {command}: {', '.join(unknown)}")
        values.update({key: (command, value) for key, value in specific.items() if key in dests})
        if not values:
            continue
        converted = {}
        for key, (section, value) in values.items():
            try:
                where = f"'{section}.{key}'" if section == COMMON else f"'{command}.{key}'"
                converted[key] = _convert(dests[key], value, where)
            except InputError as e:
                errors.append(str(e))
        for key in converted:
            dests[key].required = False
        parser.set_defaults(**converted)
    if errors:
        raise InputError("invalid configuration: " + "; ".join(dict.fromkeys(errors)))
