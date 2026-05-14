import os
from typing import Any, Union

import orjson

from swanlab.sdk.internal.pkg import safe


def json_loads(s: Union[str, bytes, bytearray]) -> Any:
    """Parse a JSON string using orjson (fast path)."""
    return orjson.loads(s)


def validate_path(base_dir: str, file_path: str) -> str | None:
    """
    Resolve *file_path* relative to *base_dir* and guard against directory traversal.

    Returns the absolute path if it lies inside *base_dir*, otherwise ``None``.
    """
    if not file_path:
        return None
    abs_base = os.path.abspath(base_dir)
    abs_path = os.path.abspath(os.path.join(base_dir, file_path))
    if not abs_path.startswith(abs_base + os.sep) and abs_path != abs_base:
        return None
    return abs_path


def proto_items_to_dict(items) -> dict:
    """
    Convert a sequence of wandb protobuf config items to a plain dict.

    Each item is expected to have ``key`` / ``nested_key`` and ``value_json``.
    JSON decode errors are silently swallowed per-item via ``safe.block``.
    """
    mapping = {}
    for item in items:
        key = item.key or "/".join(item.nested_key)
        if not key or key.startswith("_"):
            continue
        with safe.block(message=f"Could not decode json for key '{key}'"):
            mapping[key] = json_loads(item.value_json)
    return mapping
