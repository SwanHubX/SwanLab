"""
@author: cunyue
@file: context.py
@time: 2026/3/10 18:31
@description: SwanLab 运行时上下文，包含：
1. 实验配置
2. 运行时指标状态
"""

import sys
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Dict, Generator, Literal, Optional, Union

from swanlab.sdk.internal.context.transformer import TransformType
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.settings import Settings

from .callbacker import CallbackManager, callbacker


# 运行配置，包含当前运行上下文的必要状态
@dataclass(frozen=True)
class RunConfig:
    run_dir: Path
    settings: Settings


# 开启 slots 可以减少内存占用，提高性能，但是仅适用于 Python 3.10 及以上版本
# 因此我们动态判断版本：如果是 3.10 及以上，开启 slots；否则传空字典
DATACLASS_KWARGS = {"slots": True} if sys.version_info >= (3, 10) else {}


# 使用解包的方式传入参数
@dataclass(**DATACLASS_KWARGS)
class RunScalarMetric:
    # 必须是普通实例变量，初始值为 None
    _type: TransformType
    # 指标对应的单实验图表索引
    _chart: str
    # 指标对应的图表名称，默认为列名称
    _chart_name: Optional[str] = None
    # 指标名称，默认为列名称
    _name: Optional[str] = None
    # 指标颜色
    _color: Optional[str] = None
    # x轴，可以是其他的标量，也可以是系统值"_step"或"_relative_time"
    _x_axis: Optional[Union[str, Literal["_step", "_relative_time"]]] = "_step"
    # 是否为系统指标
    _system: bool = False
    # 指标最新值
    latest: Optional[Union[float, int]] = None
    # 指标的最大值
    max: Optional[Union[float, int]] = None
    # 指标的最小值
    min: Optional[Union[float, int]] = None

    def update(self, _type: "TransformType", value: Union[float, int]) -> bool:
        """
        更新指标值
        :param _type: 传入的指标类型
        :param value: 指标值
        :return: 是否成功更新
        """
        # 首次写入逻辑：如果当前没有类型，则锁定为传入的类型
        if self._type is None:
            self._type = _type

        # 后续写入逻辑：校验类型是否一致
        elif self._type != _type:
            console.error(f"Invalid type for scalar metric: expected {self._type}, got {_type}")
            return False

        # 更新值
        self.latest = value
        if self.max is None or value > self.max:
            self.max = value
        if self.min is None or value < self.min:
            self.min = value
        return True


@dataclass(**DATACLASS_KWARGS)
class RunMediaMetric:
    # 媒体类型
    _type: TransformType
    # 媒体存储路径
    path: Path

    def update(self, _type: TransformType, path: Path) -> bool:
        """
        更新媒体信息
        :param _type: 媒体类型
        :param path: 媒体存储路径
        :return: 是否成功更新
        """
        if _type != self._type:
            console.error(f"Invalid type for media metric: expected {_type}, got {_type}")
            return False
        self.path = path
        return True


# 指标状态，实验运行过程中不断更新
@dataclass
class RunMetrics:
    global_step: int = 0
    _metrics: Dict[str, Union[RunScalarMetric, RunMediaMetric]] = field(default_factory=dict)


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
    def backup_file(self) -> Path:
        return self.config.run_dir / "backup.swanlab"


# ContextVar 现在只存这个轻量级的数据宿主
_current_ctx: ContextVar[Optional[RunContext]] = ContextVar("swanlab_run_ctx", default=None)


def has_context() -> bool:
    """
    检查SwanLab运行上下文是否已初始化。
    """
    return _current_ctx.get() is not None


def set_context(ctx: RunContext):
    """
    设置SwanLab运行上下文，这将在实验开始前/开始时调用。
    如果上下文已存在，则报错
    """
    if has_context():
        raise RuntimeError("SwanLab Context is already active. Cannot set another context.")
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


@contextmanager
def use_temp_context(ctx: RunContext) -> Generator[RunContext, None, None]:
    """
    临时使用 SwanLab 运行上下文的上下文管理器。

    前置条件：仅允许在当前不存在上下文时使用。
    退出行为：无论是否发生异常，离开 with 块时都会自动清空上下文。
    """
    # 1. 严格检查：如果已经存在上下文，直接拦截报错
    if has_context():
        raise RuntimeError("SwanLab Context is already active. Cannot nest or overwrite temp contexts.")

    # 2. 前置操作：设置上下文
    set_context(ctx)

    try:
        # 3. 交出执行权，并把 ctx yield 出去，方便外部直接用 `as` 接收
        yield ctx
    finally:
        # 4. 最终清理（回退）：无论业务代码报什么错，绝对保证上下文被清空，不会污染全局 ContextVar
        clear_context()
