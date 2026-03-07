"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 12:45
@description: SwanLab 运行时上下文，存储Key等存粹的动态信息
"""

from contextvars import ContextVar
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

from swanlab.sdk.internal.settings import Settings

__all__ = ["RunConfig", "RunContext", "set_context", "clear_context", "get_context", "has_context"]


# 运行配置，包含当前运行上下文的必要状态
class RunConfig:
    settings: Settings


# 指标状态，实验运行过程中不断更新
@dataclass
class RunMetrics:
    global_step: int = 0
    _metric_steps: Dict[str, int] = field(default_factory=dict)
    summary: Dict[str, Any] = field(default_factory=dict)

    def update(self, metrics: Dict[str, Any], explicit_step: Optional[int] = None):
        if explicit_step is None:
            self.global_step += 1
        for key, value in metrics.items():
            current_step = explicit_step if explicit_step is not None else self._metric_steps.get(key, 0) + 1
            self._metric_steps[key] = current_step
            self.summary[key] = value


# 上下文宿主，纯粹的数据容器
class RunContext:
    def __init__(self, config: RunConfig):
        self.config: RunConfig = config
        self.metrics: RunMetrics = RunMetrics()


# ContextVar 现在只存这个轻量级的数据宿主
_current_ctx: ContextVar[Optional[RunContext]] = ContextVar("swanlab_run_ctx", default=None)


def set_context(ctx: RunContext):
    """
    设置SwanLab运行上下文，这将在实验开始前/开始时调用。
    """
    _current_ctx.set(ctx)


def clear_context():
    """
    清空SwanLab运行上下文，这将在每个实验结束后调用。
    """
    _current_ctx.set(None)


def get_context() -> RunContext:
    """
    获取SwanLab运行上下文。

    若上下文未初始化则抛出RuntimeError。
    """
    ctx = _current_ctx.get()
    if ctx is None:
        raise RuntimeError("SwanLab Context is not initialized.")
    return ctx


def has_context() -> bool:
    """
    检查SwanLab运行上下文是否已初始化。
    """
    return _current_ctx.get() is not None
