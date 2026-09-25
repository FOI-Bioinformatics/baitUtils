"""
logging_utils.py

One logging configuration for all commands. Level resolution, highest
priority first: verbose flag (DEBUG), quiet flag (WARNING), explicit level.
"""

import logging
from typing import Optional

LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
LEVELS = ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL")


_configured = False


def configure_logging(level: Optional[str] = "INFO", verbose: bool = False, quiet: bool = False) -> None:
    """Configure the root logger; safe to call more than once."""
    global _configured
    _configured = True
    if verbose:
        resolved = logging.DEBUG
    elif quiet:
        resolved = logging.WARNING
    else:
        resolved = getattr(logging, str(level or "INFO").upper(), logging.INFO)
    logging.basicConfig(level=resolved, format=LOG_FORMAT, datefmt=DATE_FORMAT, force=True)


def ensure_logging(level: Optional[str] = "INFO", verbose: bool = False, quiet: bool = False) -> None:
    """
    Configure logging unless the command-line dispatcher already did. Used
    by command modules so that they also log sensibly when called directly.
    """
    if not _configured:
        configure_logging(level=level, verbose=verbose, quiet=quiet)
