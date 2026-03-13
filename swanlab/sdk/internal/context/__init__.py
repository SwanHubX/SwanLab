"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 12:45
@description: SwanLab 运行时上下文，存储Key等存粹的动态信息
"""

from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Generator, Optional

from swanlab.sdk.internal.settings import Settings

from .callbacker import CallbackManager, callbacker
from .metrics import MediaMetric, RunMetrics, ScalarMetric
from .transformer import TransformMediaType, TransformType

__all__ = [
    # value
    "use_context",
    "callbacker",
    # class
    "RunContext",
    "RunConfig",
    "RunMetrics",
    "TransformType",
    "TransformMediaType",
    "CallbackManager",
    "ScalarMetric",
    "MediaMetric",
]


# 运行配置，包含当前运行上下文的必要状态
@dataclass(frozen=True)
class RunConfig:
    run_dir: Path
    settings: Settings


# 上下文宿主
class RunContext:
    def __init__(self, config: RunConfig):
        self.config: RunConfig = config
        self.callbacker: CallbackManager = CallbackManager()
        # 使用 callbacker.registered_callbacks 作为初始回调集合
        self.callbacker.merge_callbacks(callbacker.registered_callbacks)
        self.metrics: RunMetrics = RunMetrics()

    @cached_property
    def run_dir(self) -> Path:
        return self.config.run_dir

    @cached_property
    def media_dir(self) -> Path:
        return self.config.run_dir / "media"

    @cached_property
    def files_dir(self) -> Path:
        return self.config.run_dir / "files"

    @cached_property
    def metadata_file(self) -> Path:
        return self.files_dir / "swanlab-metadata.json"

    @cached_property
    def config_file(self) -> Path:
        return self.files_dir / "config.yaml"

    @cached_property
    def requirements_file(self) -> Path:
        return self.files_dir / "requirements.txt"

    @cached_property
    def conda_file(self) -> Path:
        return self.files_dir / "conda.yaml"

    @cached_property
    def debug_dir(self) -> Path:
        return self.config.run_dir / "debug"

    @cached_property
    def run_file(self) -> Path:
        run_id = self.config.settings.run.id
        assert run_id, "Run ID is not set."
        return self.config.run_dir / f"run-{run_id}.swanlab"


# ContextVar 现在只存这个轻量级的数据宿主
_current_ctx: ContextVar[Optional[RunContext]] = ContextVar("swanlab_run_ctx", default=None)


def has_context() -> bool:
    """
    检查SwanLab运行上下文是否已初始化。
    """
    return _current_ctx.get() is not None


def get_context() -> RunContext:
    """
    获取SwanLab运行上下文。

    若上下文未初始化则抛出RuntimeError。
    """
    ctx = _current_ctx.get()
    if ctx is None:
        raise RuntimeError("SwanLab Context is not initialized.")
    return ctx


@contextmanager
def use_context(ctx: RunContext) -> Generator[RunContext, None, None]:
    """
    临时使用 SwanLab 运行上下文的上下文管理器。

    前置条件：仅允许在当前不存在上下文时使用。
    退出行为：无论是否发生异常，离开 with 块时都会自动清空上下文。
    """
    # 1. 严格检查：如果已经存在上下文，直接拦截报错
    if _current_ctx.get() is not None:
        raise RuntimeError("SwanLab Context is already active. Cannot nest or overwrite temp contexts.")

    # 2. 前置操作：设置上下文
    _current_ctx.set(ctx)

    try:
        # 3. 交出执行权，并把 ctx yield 出去，方便外部直接用 `as` 接收
        yield ctx
    finally:
        # 4. 最终清理（回退）：无论业务代码报什么错，绝对保证上下文被清空，不会污染全局 ContextVar
        _current_ctx.set(None)
