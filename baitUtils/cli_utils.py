"""
cli_utils.py

Helpers shared by the command-line parsers.

Every command spells its options in kebab-case (``--min-coverage``). The
spellings used before 0.6 (``--min_coverage``, ``--minIdentity``, ``-q`` for
the reference of ``map``) are still accepted as aliases but are not shown in
the help text, so that existing scripts keep working while the documented
interface stays uniform.
"""

import argparse
from typing import Iterable

from baitUtils.logging_utils import LEVELS


def register_alias(parser: argparse.ArgumentParser, action: argparse.Action, *aliases: str) -> None:
    """
    Accept each alias as another spelling of an existing option. The alias
    is resolved to the same action but does not appear in the help text.
    """
    for alias in aliases:
        if alias in parser._option_string_actions:
            raise ValueError(f"option {alias} is already defined")
        parser._option_string_actions[alias] = action


def add_argument(parser: argparse.ArgumentParser, *flags: str,
                 aliases: Iterable[str] = (), **kwargs) -> argparse.Action:
    """parser.add_argument with hidden alias spellings."""
    action = parser.add_argument(*flags, **kwargs)
    register_alias(parser, action, *aliases)
    return action


def add_logging_arguments(parser: argparse.ArgumentParser,
                          verbose_aliases: Iterable[str] = ()) -> None:
    """
    The logging options every command has: -v/--verbose (DEBUG), --quiet
    (WARNING) and --log-level. The dispatcher reads them; command modules
    call ensure_logging with the same values when run directly.
    """
    add_argument(parser, '-v', '--verbose', action='store_true',
                 help='Debug-level logging', aliases=verbose_aliases)
    add_argument(parser, '--quiet', action='store_true',
                 help='Log warnings and errors only')
    add_argument(parser, '--log-level', dest='log_level', choices=LEVELS, default=None,
                 help='Logging level (default: INFO)', aliases=('--log_level',))


def logging_kwargs(args: argparse.Namespace) -> dict:
    """Keyword arguments for ensure_logging / configure_logging from parsed arguments."""
    return {
        'level': getattr(args, 'log_level_global', None) or getattr(args, 'log_level', None) or 'INFO',
        'verbose': bool(getattr(args, 'verbose', False)),
        'quiet': bool(getattr(args, 'quiet', False)),
    }


def percent(value: str) -> float:
    """argparse type for a percentage between 0 and 100."""
    number = float(value)
    if not 0 <= number <= 100:
        raise argparse.ArgumentTypeError(f"{value} is not between 0 and 100")
    return number


def positive_int(value: str) -> int:
    """argparse type for an integer of at least 1."""
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError(f"{value} must be at least 1")
    return number


def non_negative_int(value: str) -> int:
    """argparse type for an integer of at least 0."""
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError(f"{value} must not be negative")
    return number
