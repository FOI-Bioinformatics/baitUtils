"""
json_export.py

Helpers for writing analysis results as JSON with numpy and dataclass values
converted to plain Python types.
"""

import dataclasses
import json
import shutil
import subprocess
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd


def to_jsonable(value: Any) -> Any:
    """Recursively convert numpy, pandas, dataclass, NamedTuple and Enum values."""
    if isinstance(value, dict):
        return {str(k): to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)) and not hasattr(value, '_asdict'):
        return [to_jsonable(v) for v in value]
    if hasattr(value, '_asdict'):
        return to_jsonable(value._asdict())
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return to_jsonable(dataclasses.asdict(value))
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, pd.DataFrame):
        return to_jsonable(value.to_dict(orient='records'))
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if np.isnan(value) else float(value)
    if isinstance(value, float) and value != value:
        return None
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def _first_line(cmd, timeout: float = 10.0) -> Optional[str]:
    """First non-empty output line of a command, or None when unavailable."""
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired):
        return None
    for line in result.stdout.splitlines():
        if line.strip():
            return line.strip()
    return None


def tool_versions() -> Dict[str, Optional[str]]:
    """Versions of baitUtils, Python, key libraries and external tools on PATH."""
    import platform
    import numpy
    import pandas
    import scipy
    import Bio
    from baitUtils._version import __version__

    versions: Dict[str, Optional[str]] = {
        "baitUtils": __version__,
        "python": platform.python_version(),
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
        "scipy": scipy.__version__,
        "biopython": Bio.__version__,
    }
    try:
        import RNA
        versions["viennarna"] = getattr(RNA, "__version__", "present")
    except ImportError:
        versions["viennarna"] = None
    versions["pblat"] = _first_line(["pblat"]) if shutil.which("pblat") else None
    versions["minimap2"] = _first_line(["minimap2", "--version"]) if shutil.which("minimap2") else None
    versions["bedtools"] = _first_line(["bedtools", "--version"]) if shutil.which("bedtools") else None
    return versions


def arguments_record(args: Any) -> Dict[str, Any]:
    """Plain dictionary of parsed command-line arguments without callables."""
    return {key: to_jsonable(value) for key, value in vars(args).items() if not callable(value)}


def write_json(data: Any, path: Path) -> Path:
    """Write data as indented JSON and return the path."""
    path = Path(path)
    with open(path, 'w') as fh:
        json.dump(to_jsonable(data), fh, indent=2)
    return path
