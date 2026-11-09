"""
@author: Nexisato
@file: locate.py
@time: 2026/11/9
@description: 定位 swanlab-core：环境变量优先，其次包内 swanlab/bin。
"""

import os
from pathlib import Path

from swanlab.sdk.internal.core_go.config import BINARY_NAME
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage

__all__ = ["CORE_BINARY_ENV", "locate_core_binary"]

CORE_BINARY_ENV = "SWANLAB_CORE_BINARY"

_IS_WINDOWS = os.name == "nt"


def _binary_name() -> str:
    return f"{BINARY_NAME}.exe" if _IS_WINDOWS else BINARY_NAME


def _packaged_binary() -> Path:
    # swanlab/sdk/internal/core_go/locate.py -> swanlab/bin/
    return Path(__file__).resolve().parents[3] / "bin" / _binary_name()


def _missing(path: Path, *, source: str) -> CoreGoError:
    return CoreGoError(
        Stage.SPAWN,
        ReasonCode.BINARY_NOT_FOUND,
        f"swanlab-core binary not found: {path} ({source}); "
        f"set {CORE_BINARY_ENV} or install a platform wheel that bundles swanlab-core",
    )


def locate_core_binary() -> Path:
    """返回可用的 core 二进制路径；找不到时抛出 CoreGoError。"""
    override = os.environ.get(CORE_BINARY_ENV, "").strip()
    if override:
        path = Path(override).expanduser()
        if not path.is_file():
            raise _missing(path, source=f"from {CORE_BINARY_ENV}")
        return path.resolve()

    packaged = _packaged_binary()
    if packaged.is_file():
        return packaged
    raise _missing(packaged, source="bundled binary")
