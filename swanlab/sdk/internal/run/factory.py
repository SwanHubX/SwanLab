"""
@author: cunyue
@file: factory.py
@time: 2026/4/14
@description: 运行时组件工厂，根据 mode 返回真实组件或 Null 组件
"""

from queue import Queue
from typing import Optional

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.sdk.internal.bus import RunEmitter
from swanlab.sdk.internal.bus.emitter import EmitterProtocol, RunQueue
from swanlab.sdk.internal.bus.events import CondaEvent, EventPayload, MetadataEvent, RequirementsEvent
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import fs
from swanlab.sdk.internal.run import system

from .config import Config, create_run_config, create_unbound_run_config
from .consumer import BackgroundConsumer, ConsumerProtocol
from .record_builder import RecordBuilder


class NullEmitter(EmitterProtocol):
    """空发射器，emit 为 no-op"""

    def __init__(self):
        self._queue = Queue()

    def emit(self, event: EventPayload) -> None:
        pass

    @property
    def queue(self) -> RunQueue:
        return self._queue


def factory_emitter(ctx: RunContext) -> EmitterProtocol:
    """emitter对象工厂

    :param ctx: 运行上下文，包含配置信息和运行时状态
    """
    if ctx.config.settings.mode == "disabled":
        return NullEmitter()
    return RunEmitter()


class NullConsumer(ConsumerProtocol):
    """空消费者，start/join 为 no-op"""

    def start(self) -> None:
        pass

    def stop(self) -> None:
        pass

    def join(self) -> None:
        pass


def factory_consumer(
    ctx: RunContext,
    emitter: EmitterProtocol,
    builder: RecordBuilder,
) -> ConsumerProtocol:
    """consumer对象工厂

    :param ctx: 运行上下文
    :param emitter: 事件发射器
    :param builder: 构建器
    """
    if ctx.config.settings.mode == "disabled":
        return NullConsumer(ctx, emitter.queue, builder)
    return BackgroundConsumer(ctx, emitter.queue, builder)


def factory_config(ctx: RunContext, emitter: EmitterProtocol) -> Config:
    """
    config 对象工厂函数

    :param ctx: 运行上下文
    :param emitter: 事件发射器，绑定 config 时用于发出 ConfigEvent
    """

    if ctx.config.settings.mode == "disabled":
        return create_unbound_run_config()
    return create_run_config(ctx.config_file, emitter.emit)


def factory_monitor(ctx: RunContext, emitter: EmitterProtocol) -> Optional["system.Monitor"]:
    """
    系统信息采集函数
    :param ctx:
    :param emitter:
    :return:
    """
    if ctx.config.settings.mode == "disabled":
        return None
    this_monitor: Optional["system.Monitor"] = None
    sys_info, monitor = system.new(ctx)
    ts = Timestamp()
    ts.GetCurrentTime()
    if sys_info.metadata:
        fs.safe_write(ctx.metadata_file, sys_info.metadata.model_dump_json())
        emitter.emit(MetadataEvent(timestamp=ts))
    if sys_info.requirements:
        fs.safe_write(ctx.requirements_file, sys_info.requirements)
        emitter.emit(RequirementsEvent(timestamp=ts))
    if sys_info.conda:
        fs.safe_write(ctx.conda_file, sys_info.conda)
        emitter.emit(CondaEvent(timestamp=ts))
    if monitor is not None and monitor.start(ctx, emitter):
        this_monitor = monitor
    return this_monitor
