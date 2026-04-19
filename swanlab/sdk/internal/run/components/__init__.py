"""
@author: cunyue
@file: __init__.py
@time: 2026/4/19 18:59
@description: Run 组件工厂，提供 Run 所需组件的创建和初始化
"""

from typing import Tuple

from swanlab.sdk.internal.bus import RunEmitter
from swanlab.sdk.internal.bus.emitter import EmitterProtocol
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.run.components.asynctask import AsyncTaskManager
from swanlab.sdk.internal.run.components.builder import RecordBuilder
from swanlab.sdk.internal.run.components.config import Config, create_run_config, create_unbound_run_config
from swanlab.sdk.internal.run.components.consumer import BackgroundConsumer, ConsumerProtocol
from swanlab.sdk.internal.run.components.null import NullConsumer, NullEmitter

__all__ = [
    "new",
    "Config",
    "ConsumerProtocol",
    "EmitterProtocol",
    "AsyncTaskManager",
    "RecordBuilder",
    "BackgroundConsumer",
    "NullConsumer",
    "NullEmitter",
]


def new(ctx: RunContext) -> Tuple[AsyncTaskManager, ConsumerProtocol, EmitterProtocol, Config]:
    """
    创建 Run 运行所需必要组件
    :param ctx: 运行上下文
    :return: 任务管理器，消费者，系统监控，配置
    """
    a = AsyncTaskManager()
    b = RecordBuilder(ctx)
    e = _factory_emitter(ctx)
    c = _factory_consumer(ctx, e, b)
    return a, c, e, _factory_config(ctx, e)


def _factory_emitter(ctx: RunContext) -> EmitterProtocol:
    """emitter对象工厂

    :param ctx: 运行上下文，包含配置信息和运行时状态
    """
    if ctx.config.settings.mode == "disabled":
        return NullEmitter()
    return RunEmitter()


def _factory_consumer(
    ctx: RunContext,
    e: EmitterProtocol,
    b: RecordBuilder,
) -> ConsumerProtocol:
    """consumer对象工厂

    :param ctx: 运行上下文
    :param e: 事件发射器
    :param b: 构建器
    """
    if ctx.config.settings.mode == "disabled":
        return NullConsumer(ctx, e.queue, b)
    return BackgroundConsumer(ctx, e.queue, b)


def _factory_config(ctx: RunContext, e: EmitterProtocol) -> Config:
    """
    config 对象工厂函数

    :param ctx: 运行上下文
    :param e: 事件发射器，绑定 config 时用于发出 ConfigEvent
    """

    if ctx.config.settings.mode == "disabled":
        return create_unbound_run_config()
    return create_run_config(ctx.config_file, e.emit)
