"""
bedtools_support.py

Dependency check for the commands that rely on pybedtools and the bedtools
binary (check and fill). Other commands do not need bedtools.
"""

import logging
import shutil
import sys

INSTALL_HINT = (
    "The check and fill commands require pybedtools and the bedtools binary. "
    "Install with: pip install pybedtools; conda install -c bioconda bedtools"
)


def has_bedtools() -> bool:
    """Return True when both pybedtools and the bedtools binary are available."""
    try:
        import pybedtools  # noqa: F401
    except ImportError:
        return False
    return shutil.which("bedtools") is not None


def require_bedtools() -> None:
    """Exit with a clear message when bedtools support is missing."""
    if not has_bedtools():
        logging.error(INSTALL_HINT)
        sys.exit(INSTALL_HINT)
