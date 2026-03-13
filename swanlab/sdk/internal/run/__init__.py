"""
@author: cunyue
@file: __init__.py
@time: 2026/3/12
@description: SwanLab SDK 运行模块，涉及：
1. 数据处理 (基于 Event-Bus 事件驱动架构)
2. 运行、实验上下文管理
3. 触发异步微批处理落盘与回调
"""

import threading
from functools import cached_property, wraps
from pathlib import Path
from typing import Any, Literal, Mapping, Optional, Union, cast, get_args

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.sdk.internal.bus import RunEmitter
from swanlab.sdk.internal.bus.events import MetricDefineEvent, MetricLogEvent, RunFinishEvent, RunStartEvent
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.core_python import CorePython
from swanlab.sdk.internal.pkg import console, log
from swanlab.sdk.typings.run import FinishType
from swanlab.sdk.typings.run.data import ScalarXAxisType

from . import utils_fmt as fmt
from .consumer import BackgroundConsumer
from .data.transforms import Text
from .record_builder import RecordBuilder

__all__ = ["SwanLabRun", "has_run", "get_run", "set_run", "clear_run"]


class SwanLabRun:
    """
    The SwanLabRun class is used for logging during a single run.
    There should be only one instance of the SwanLabRun class for each experiment.
    """

    def __init__(self, ctx: RunContext):
        self._ctx = ctx
        self._state: Union[FinishType, Literal["running"]] = "running"
        # 外部API锁，防止并发调用
        self._api_lock = threading.RLock()

        # 事件发射器：唯一的队列写入入口，可注入给内部系统组件
        self._emitter = RunEmitter(maxsize=100_000)

        # Core：Record 落盘与后端交互的统一入口
        self._core = CorePython(ctx)

        # 记录构建器：负责将事件转换为 Record
        self._builder = RecordBuilder(ctx)

        # 后台消费者：从 emitter.queue 消费事件并落盘
        self._consumer = BackgroundConsumer(ctx, self._emitter.queue, self._builder, self._core)

        # 触发启动事件
        ts = Timestamp()
        ts.GetCurrentTime()
        self._emitter.emit(RunStartEvent(timestamp=ts))

        # TODO: 硬件监控与metadata采集
        # TODO: Config 事件

        # 启动后台消费者
        self._consumer.start()

        # 设置全局运行实例
        set_run(self)

        # 绑定日志文件
        if self._ctx.config.settings.mode != "disabled":
            log.bindfile(self._ctx.debug_dir)

    # ----------------------------------
    # 属性 (Properties)
    # ----------------------------------

    @cached_property
    def id(self) -> str:
        assert self._ctx.config.settings.run.id is not None, "Run id is not set."
        return self._ctx.config.settings.run.id

    @cached_property
    def run_dir(self) -> Path:
        assert self._ctx.run_dir is not None, "Run dir is not set."
        return self._ctx.run_dir

    @cached_property
    def project_url(self) -> Optional[str]:
        settings = self._ctx.config.settings
        if settings.mode != "cloud":
            return None
        return f"{settings.web_host}/@{settings.project.workspace}/{settings.project.name}"

    @cached_property
    def url(self) -> str:
        settings = self._ctx.config.settings
        if settings.mode != "cloud":
            raise RuntimeError("Run URL is only available in cloud mode.")
        return f"{self.project_url}/runs/{settings.run.id}"

    # ----------------------------------
    # 私有 API：内部辅助函数
    # ----------------------------------
    @staticmethod
    def _with_lock(func):
        """线程安全装饰器，自动获取和释放外部API锁"""

        @wraps(func)
        def wrapper(self, *args, **kwargs):
            with self._api_lock:
                return func(self, *args, **kwargs)

        return wrapper

    # ----------------------------------
    # 公开 API：只负责验证输入并发事件
    # ----------------------------------
    @_with_lock
    def log(self, data: Mapping[str, Any], step: Optional[int] = None):
        """记录一组日志（可能触发隐式列创建）"""
        if self._state != "running":
            console.error("Run has already finished or is not active, cannot call log() again.")
            return
        if not (this_data := fmt.safe_validate_log_data(data)):
            console.error(f"Log data must be a dict, but got {type(data).__name__}. SwanLab will ignore this log.")
            return
        if step is not None:
            if not isinstance(step, int):
                console.error(
                    f"Step must be an integer or None, but got {type(step).__name__}. SwanLab will ignore this log."
                )
                return
            if step < 0:
                console.error(f"Step must be non-negative, but got {step}. SwanLab will ignore this log.")
                return

        next_step = self._ctx.metrics.next_step(step)

        ts = Timestamp()
        ts.GetCurrentTime()

        # 展平字典并在内部进行合规性验证和截断
        flatten_data = fmt.flatten_dict(this_data)

        # 推送日志事件
        self._emitter.emit(MetricLogEvent(data=flatten_data, step=next_step, timestamp=ts))

    @_with_lock
    def log_text(self, key: str, data: Union[str, Text], caption: Optional[str] = None, step: Optional[int] = None):
        """
        A syntactic sugar for logging text data.
        :param key: The key for the text data.
        :param data: The text data itself or a Text object.
        :param caption: Optional caption for the text data.
        :param step: Optional step for the text data.
        """
        if not isinstance(data, Text):
            data = Text(data, caption=caption)
        self.log({key: data}, step=step)

    @_with_lock
    def define_scalar(
        self,
        key: str,
        name: Optional[str] = None,
        color: Optional[str] = None,
        x_axis: Optional[ScalarXAxisType] = None,
        chart_name: Optional[str] = None,
    ):
        """
        Explicitly define a scalar column.
        :param key: The key for the scalar column.
        :param name: Optional name for the scalar column.
        :param color: Optional color for the scalar column.
        :param x_axis: Optional x-axis for the scalar column.
        :param chart_name: Optional name for the chart.
        """
        if not (this_key := fmt.safe_validate_key(key)):
            return console.error(
                f"Invalid key for define scalar: {key}, please use valid characters (alphanumeric, '.', '-', '/') and avoid special characters."
            )

        original_name = name
        if name and not (name := fmt.safe_validate_name(name)):
            return console.error(f"Invalid name for define scalar: {original_name}, must be a string.")

        original_color = color
        if color and not (color := fmt.safe_validate_color(color)):
            return console.error(f"Invalid color for define scalar: {original_color}, must be a hex color code.")

        if (this_x_axis := fmt.safe_validate_x_axis(x_axis)) is None:
            return console.error(f"Invalid x_axis for define scalar: {x_axis}, must be a valid ScalarXAxisType.")

        original_chart_name = chart_name
        if chart_name and not (chart_name := fmt.safe_validate_chart_name(chart_name)):
            return console.error(f"Invalid chart_name for define scalar: {original_chart_name}, must be a string.")

        self._emitter.emit(
            MetricDefineEvent(
                key=this_key,
                name=name,
                color=color,
                system=False,
                x_axis=this_x_axis,
                chart_name=chart_name,
                chart=None,
            )
        )

    @_with_lock
    def finish(self, state: FinishType = "success", error: Optional[str] = None):
        """安全关闭当前 Run，等待所有日志落盘"""
        if self._state != "running":
            console.error("Run has already finished or is not active, cannot call finish() again.")
            return
        self._state = state
        state = state.lower()  # type: ignore
        if not (this_state := fmt.safe_validate_state(cast(FinishType, state))):
            console.error(f"Invalid state: {state}, allowed values are {get_args(FinishType)}")
            return

        if state == "crashed" and error is None:
            console.warning("Crashed reason has been set to 'unknown' due to missing error message.")
            error = "unknown"

        console.debug("SwanLab Run is finishing, waiting for logs to flush...")
        # 线程退出
        ts = Timestamp()
        ts.GetCurrentTime()
        self._emitter.emit(RunFinishEvent(state=this_state, error=error, timestamp=ts))
        # 阻塞主线程，等待后台队列消费完毕
        self._consumer.join()
        # 清理全局运行实例
        clear_run()
        console.debug(f"Run finished with state: {state}")
        # 释放全局logger
        if self._ctx.config.settings.mode != "disabled":
            log.reset()


_current_run: Optional[SwanLabRun] = None


def has_run() -> bool:
    return _current_run is not None


def get_run() -> SwanLabRun:
    if _current_run is None:
        raise RuntimeError("No active SwanLabRun. Call swanlab.init() first.")
    return _current_run


def set_run(run: SwanLabRun) -> None:
    global _current_run
    _current_run = run


def clear_run() -> None:
    global _current_run
    _current_run = None
