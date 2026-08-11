"""
@author: cunyue
@file: impl.py
@time: 2026/5/17
@description: SwanLab SDK component implementation helpers

NOTE: this is a transition module.
"""

import os
from typing import Callable, Dict, Literal, Optional, Tuple

from swanlab.sdk.protocol import CoreEnum, ProbeEnum

require_env = os.environ.get("SWANLAB_REQUIRE", "")
# probe/core 实现选择，进程级全局状态，默认 Python 实现
# 过渡期通过 swanlab.require("core"|"probe") 或 SWANLAB_REQUIRE 环境变量切换
_probe_impl: ProbeEnum = ProbeEnum.PROBE_PYTHON
_core_impl: CoreEnum = CoreEnum.CORE_PYTHON


def set_probe_impl(impl: ProbeEnum) -> None:
    """设置当前 probe 实现类型（进程级，需在 swanlab.init 之前调用）。"""
    global _probe_impl
    _probe_impl = impl


def set_core_impl(impl: CoreEnum) -> None:
    """设置当前 core 实现类型（进程级，需在 swanlab.init 之前调用）。"""
    global _core_impl
    _core_impl = impl


def get_probe_impl() -> Literal["python", "rust"]:
    """获取当前 probe 实现类型。"""
    if _probe_impl == ProbeEnum.PROBE_PYTHON:
        return "python"
    elif _probe_impl == ProbeEnum.PROBE:
        return "rust"
    else:
        raise ValueError(f"Unknown probe impl: {_probe_impl}")


def get_core_impl() -> Literal["python", "go"]:
    """获取当前 core 实现类型。"""
    if _core_impl == CoreEnum.CORE_PYTHON:
        return "python"
    elif _core_impl == CoreEnum.CORE:
        return "go"
    else:
        raise ValueError(f"Unknown core impl: {_core_impl}")


# requirement token -> (人类可读说明, 设置器)
# 过渡态：通过 swanlab.require(...) 或 SWANLAB_REQUIRE 环境变量切换后端实现
REQUIREMENTS: Dict[str, Tuple[str, Callable[[], None]]] = {
    "core": ("SwanLab Go core (swanlab-core)", lambda: set_core_impl(CoreEnum.CORE)),
    "probe": ("SwanLab Rust probe", lambda: set_probe_impl(ProbeEnum.PROBE)),
    "core_python": ("SwanLab Python core (built-in)", lambda: set_core_impl(CoreEnum.CORE_PYTHON)),
    "probe_python": ("SwanLab Python probe (built-in)", lambda: set_probe_impl(ProbeEnum.PROBE_PYTHON)),
}


def _try_apply(token: str) -> Optional[str]:
    """已知 token 则执行设置器并返回说明；未知返回 None（静默跳过）。"""
    entry = REQUIREMENTS.get(token)
    if entry is None:
        return None
    label, setter = entry
    setter()
    return label


# 导入时读取 SWANLAB_REQUIRE 环境变量并静默应用，逗号分隔
for _t in (_x.strip().lower() for _x in require_env.split(",") if _x.strip()):
    _try_apply(_t)


def apply_requirement(token: str) -> str:
    """应用单个 requirement token，返回人类可读说明；未知 token 抛 ValueError。"""
    label = _try_apply(token)
    if label is None:
        raise ValueError(f"Unknown requirement: {token!r}. Available: {sorted(REQUIREMENTS)}")
    return label
