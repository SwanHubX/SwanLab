"""
@author: cunyue
@file: require.py
@time: 2026/8/11
@description: swanlab.require 方法，过渡态用于动态选择 core/probe 后端实现

设计参考 wandb.require，属于中间过渡态 API：
- 默认使用 Go core（CoreGo，二进制不可用时自动回退 CorePython）和 Python probe（ProbePython）
- 在 swanlab.init 之前调用 swanlab.require(...) 切换后端实现：
  - "core" / "probe"：使用新后端（Go core / Rust probe）
  - "core_python" / "probe_python"：显式切回 Python 内置实现
- 也可通过环境变量 SWANLAB_REQUIRE=core_python 等在导入时生效（应用逻辑见 helper.impl）
- 待全后端切换稳定后，此 API 将退化为 no-op 或废弃
"""

from typing import List

from swanlab.sdk.cmd.guard import with_cmd_lock, without_run
from swanlab.sdk.internal.pkg import console, helper

__all__ = ["require"]


@with_cmd_lock
@without_run("require")
def require(*requirements: str) -> List[str]:
    """Select a SwanLab runtime backend (process-global).

    Must be called before ``swanlab.init()``. By default, SwanLab uses the Go core
    (falling back to the Python core if the binary is unavailable). Supported requirements:

    - ``"core"``: use the SwanLab Go core runtime (default).
    - ``"core_python"``: switch back to the built-in Python core.
    - ``"probe"``: switch to the SwanLab Rust probe runtime.
    - ``"probe_python"``: explicitly use the built-in Python probe.

    Equivalent to setting the ``SWANLAB_REQUIRE`` environment variable (comma-separated).
    Calling ``require("core")`` while Go core is already active emits a warning.

    :param requirements: One or more requirement tokens, e.g. ``"core_python"``, ``"core"``.
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
        was_go = helper.get_core_impl() == "go"
        label = helper.apply_requirement(token)
        if token == "core" and was_go:
            console.warning("`swanlab.require('core')` is a no-op as Go core is now the default behavior.")
        console.debug(f"require('{token}') -> {label}")
        activated.append(token)
    return activated
