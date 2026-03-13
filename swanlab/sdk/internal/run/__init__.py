"""
@author: cunyue
@file: __init__.py
@time: 2026/3/12
@description: SwanLab SDK 运行模块，涉及：
1. 数据处理 (基于 Event-Bus 事件驱动架构)
2. 运行、实验上下文管理
3. 触发异步微批处理落盘与回调
"""

import queue
from functools import cached_property
from pathlib import Path
from typing import Any, Mapping, Optional, Union, cast, get_args

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run import FinishType
from swanlab.sdk.typings.run.data import ScalarXAxisType

from . import utils_fmt as F
from .callbackers import CloudCallback, LocalCallback, OfflineCallback
from .consumer import BackgroundConsumer
from .data.transforms import Text
from .events import DefineEvent, EventPayload, FinishEvent, LogEvent
from .record import RecordBuilder

__all__ = ["SwanLabRun", "CloudCallback", "LocalCallback", "OfflineCallback"]


class SwanLabRun:
    def __init__(self, ctx: RunContext):
        self._ctx = ctx

        # 事件总线：最大积压 10 万个事件（防 OOM 兜底）
        self._queue: queue.Queue[EventPayload] = queue.Queue(maxsize=100_000)

        # 记录构建器：负责将事件转换为 Record
        self._builder = RecordBuilder(ctx)
        # 后台消费者：负责将 Record 刷盘到文件系统
        self._consumer = BackgroundConsumer(self._queue, self._builder, ctx.callbacker, ctx.metrics)
        self._consumer.start()
        # TODO: 触发启动事件
        # TODO: 硬件监控与metadata采集
        # TODO: Config 事件

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

    # ----------------------------------
    # 公开 API：只负责验证输入并发事件
    # ----------------------------------

    def log(self, data: Mapping[str, Any], step: Optional[int] = None):
        """记录一组日志（可能触发隐式列创建）"""
        if not (this_data := F.safe_validate_log_data(data)):
            console.error(f"Log data must be a dict, but got {type(data).__name__}. SwanLab will ignore it.")
            return
        if not (this_step := F.safe_validate_step(step)):
            console.error(f"Step must be an integer or None, but got {type(step).__name__}. SwanLab will ignore it.")
            return

        next_step = self._ctx.metrics.next_step(this_step)

        ts = Timestamp()
        ts.GetCurrentTime()

        # 展平字典并在内部进行合规性验证和截断
        flatten_data = F.flatten_dict(this_data)

        # 推送日志事件
        self._queue.put(LogEvent(data=flatten_data, step=next_step, timestamp=ts), block=True)

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
        if not (this_key := F.safe_validate_key(key)):
            return console.error(
                f"Invalid key for define scalar: {key}, please use valid characters (alphanumeric, '.', '-', '/') and avoid special characters."
            )

        original_name = name
        if name and not (name := F.safe_validate_name(name)):
            return console.error(f"Invalid name for define scalar: {original_name}, must be a string.")

        original_color = color
        if color and not (color := F.safe_validate_color(color)):
            return console.error(f"Invalid color for define scalar: {original_color}, must be a hex color code.")

        if (this_x_axis := F.safe_validate_x_axis(x_axis)) is None:
            return console.error(f"Invalid x_axis for define scalar: {x_axis}, must be a valid ScalarXAxisType.")

        original_chart_name = chart_name
        if chart_name and not (chart_name := F.safe_validate_chart_name(chart_name)):
            return console.error(f"Invalid chart_name for define scalar: {original_chart_name}, must be a string.")

        self._define_scalar(
            key=this_key,
            name=name,
            color=color,
            x_axis=this_x_axis,
            system=False,
            chart_name=chart_name,
        )

    def finish(self, state: FinishType = "success", error: Optional[str] = None):
        """安全关闭当前 Run，等待所有日志落盘"""
        state = state.lower()  # type: ignore
        if not (this_state := F.safe_validate_state(cast(FinishType, state))):
            console.error(f"Invalid state: {state}, allowed values are {get_args(FinishType)}")
            return

        if state == "crashed" and error is None:
            console.warning("Crashed reason has been set to 'unknown' due to missing error message.")
            error = "unknown"

        console.debug("SwanLab Run is finishing, waiting for logs to flush...")
        # 线程退出
        ts = Timestamp()
        ts.GetCurrentTime()
        self._queue.put(FinishEvent(state=this_state, error=error, timestamp=ts))
        # 阻塞主线程，等待后台队列消费完毕
        self._consumer.join()
        console.debug(f"Run finished with state: {state}")

    def _define_scalar(
        self,
        key: str,
        name: Optional[str],
        color: Optional[str],
        x_axis: ScalarXAxisType = "_step",
        system: bool = False,
        chart_name: Optional[str] = None,
        chart_index: Optional[str] = None,
    ):
        """
        定义一个标量指标
        内部方法，供系统指标使用
        :param key: 指标键
        :param name: 指标名称，默认为列名称
        :param color: 指标颜色
        :param x_axis: x轴，可以是其他的标量，也可以是系统值"_step"或"_relative_time"
        :param system: 是否为系统指标
        :param chart_name: 图表名称，默认为列名称
        :param chart_index: 图表索引
        """
        self._queue.put(
            DefineEvent(
                key=key,
                name=name,
                color=color,
                system=system,
                x_axis=x_axis,
                chart_name=chart_name,
                chart=chart_index,
            ),
            block=True,
        )
