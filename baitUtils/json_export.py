"""
json_export.py

Helpers for writing analysis results as JSON with numpy and dataclass values
converted to plain Python types.
"""

import dataclasses
import json
from enum import Enum
from pathlib import Path
from typing import Any

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


def write_json(data: Any, path: Path) -> Path:
    """Write data as indented JSON and return the path."""
    path = Path(path)
    with open(path, 'w') as fh:
        json.dump(to_jsonable(data), fh, indent=2)
    return path
