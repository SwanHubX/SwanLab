"""
@author: cunyue
@file: require.py
@time: 2026/8/11
@description: swanlab.require 方法，过渡态用于动态选择 core/probe 后端实现

设计参考 wandb.require，属于中间过渡态 API：
- 默认使用 Python 内置实现（CorePython / ProbePython）
- 在 swanlab.init 之前调用 swanlab.require(...) 切换后端实现：
  - "core" / "probe"：切换到新后端（Go core / Rust probe）
  - "core_python" / "probe_python"：显式切回 Python 内置实现
- 也可通过环境变量 SWANLAB_REQUIRE=core,probe 在导入时生效（应用逻辑见 helper.impl）
- 待新后端成为默认实现后，此 API 将退化为 no-op 或废弃
"""

from typing import List

from swanlab.sdk.cmd.guard import with_cmd_lock, without_run
from swanlab.sdk.internal.pkg import console, helper

__all__ = ["require"]


@with_cmd_lock
@without_run("require")
def require(*requirements: str) -> List[str]:
    """Enable a transitional SwanLab runtime backend (process-global).

    Must be called before ``swanlab.init()``. Currently supported requirements:

    - ``"core"``: switch to the SwanLab Go core runtime (swanlab-core).
    - ``"probe"``: switch to the SwanLab Rust probe runtime.
    - ``"core_python"``: explicitly use the built-in Python core.
    - ``"probe_python"``: explicitly use the built-in Python probe.

    Equivalent to setting the ``SWANLAB_REQUIRE`` environment variable (comma-separated).
    Idempotent; unknown tokens raise ``ValueError``.

    :param requirements: One or more requirement tokens, e.g. ``"core"``, ``"probe"``.
    :return: The list of activated requirement tokens.
    :raises ValueError: If a requirement token is unknown.
    :raises RuntimeError: If called while a run is active.

    Examples:

        >>> import swanlab
        >>> swanlab.require("core")
        >>> swanlab.init()
    """
    activated: List[str] = []
    for raw in requirements:
        token = raw.strip().lower()
        label = helper.apply_requirement(token)
        console.debug(f"require('{token}') -> {label}")
        activated.append(token)
    return activated
