"""
@author: cunyue
@file: __init__.py
@time: 2026/4/19 18:59
@description: Run 组件工厂与生命周期管理
"""

from __future__ import annotations

from typing import Optional

from swanlab.sdk.internal.bus import EmitterProtocol, RunEmitter
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console, fork
from swanlab.sdk.internal.run import system
from swanlab.sdk.internal.run.components.asynctask import AsyncTaskManager
from swanlab.sdk.internal.run.components.builder import RecordBuilder
from swanlab.sdk.internal.run.components.config import (
    Config,
    create_run_config,
    create_unbound_run_config,
    deactivate_run_config,
)
from swanlab.sdk.internal.run.components.consumer import BackgroundConsumer, ConsumerProtocol
from swanlab.sdk.internal.run.components.null import NullConsumer, NullEmitter, NullTerminalProxy
from swanlab.sdk.internal.run.components.terminal import TerminalProxy, TerminalProxyProtocol

__all__ = [
    "Components",
    "Config",
    "ConsumerProtocol",
    "TerminalProxyProtocol",
    "AsyncTaskManager",
    "RecordBuilder",
    "BackgroundConsumer",
    "NullConsumer",
    "NullEmitter",
    "NullTerminalProxy",
]


class Components:
    """Run 运行时组件容器，统一管理所有组件的创建与生命周期。

    职责：
    - 创建所有运行时组件（工厂模式）
    - 编排组件启动顺序
    - 编排组件停止顺序

    Run 只持有 Components 引用，通过属性访问需要的组件。
    """

    def __init__(self, ctx: RunContext) -> None:
        self._ctx = ctx
        self._init_pid = fork.current_pid()

        # 核心组件
        self._asynctask = AsyncTaskManager()
        self._builder = RecordBuilder(ctx)
        self._emitter = _factory_emitter(ctx)
        self._consumer = _factory_consumer(ctx, self._emitter, self._builder)
        self._config = _factory_config(ctx, self._emitter)

        # 附加组件（延迟创建，在 start() 中初始化）
        self._monitor: Optional[system.Monitor] = None
        self._terminal: TerminalProxyProtocol = _factory_terminal(ctx, self._emitter, self._init_pid)

    # ----------------------------------
    # 组件访问
    # ----------------------------------

    @property
    def asynctask(self) -> AsyncTaskManager:
        return self._asynctask

    @property
    def consumer(self) -> ConsumerProtocol:
        return self._consumer

    @property
    def emitter(self) -> EmitterProtocol:
        return self._emitter

    @property
    def config(self) -> Config:
        return self._config

    @property
    def monitor(self) -> Optional[system.Monitor]:
        return self._monitor

    @property
    def terminal(self) -> TerminalProxyProtocol:
        return self._terminal

    # ----------------------------------
    # 生命周期：启动
    # ----------------------------------

    def start(self) -> None:
        """启动所有组件，按依赖顺序编排。"""
        ctx = self._ctx

        # 1. 启动硬件监控（依赖 emitter）
        self._monitor = system.create_monitor(ctx, self._emitter)

        # 2. 启动终端代理（依赖 emitter）
        self._terminal.install()

        # 3. 启动消费者线程（消费 emitter 队列）
        self._consumer.start()

    # ----------------------------------
    # 生命周期：停止
    # ----------------------------------

    def stop(self, async_log_timeout: float | None = None) -> None:
        """停止所有组件，按反向依赖顺序编排。

        :param async_log_timeout: async_log 任务等待超时
        """
        # 1. 等待异步任务完成
        console.debug("Waiting for async_log tasks to complete...")
        self._asynctask.shutdown(timeout=async_log_timeout)

        # 2. 停止硬件监控
        if self._monitor is not None:
            console.debug("Stopping hardware monitor...")
            self._monitor.stop()

        # 3. 停止终端代理（可能发射最后的 ConsoleEvent）
        console.debug("Stopping terminal proxy...")
        self._terminal.uninstall()

        # 4. 解绑 config
        deactivate_run_config()

        # 5. 停止消费者线程（消费剩余事件包括最后的 ConsoleEvent）
        console.debug("SwanLab Run is finishing, waiting for logs to flush...")
        self._consumer.stop()
        self._consumer.join()


# ==========================================
# 工厂函数
# ==========================================


def _factory_emitter(ctx: RunContext) -> EmitterProtocol:
    if ctx.config.settings.mode == "disabled":
        return NullEmitter()
    return RunEmitter()


def _factory_consumer(
    ctx: RunContext,
    e: EmitterProtocol,
    b: RecordBuilder,
) -> ConsumerProtocol:
    if ctx.config.settings.mode == "disabled":
        return NullConsumer(ctx, e.queue, b)
    return BackgroundConsumer(ctx, e.queue, b)


def _factory_config(ctx: RunContext, e: EmitterProtocol) -> Config:
    if ctx.config.settings.mode == "disabled":
        return create_unbound_run_config()
    return create_run_config(ctx.config_file, e.emit)


def _factory_terminal(ctx: RunContext, e: EmitterProtocol, init_pid: int) -> TerminalProxyProtocol:
    settings = ctx.config.settings
    if settings.mode == "disabled" or settings.console.proxy_type == "none":
        return NullTerminalProxy()
    return TerminalProxy(
        emitter=e,
        proxy_type=settings.console.proxy_type,
        max_log_length=settings.console.max_log_length,
        init_pid=init_pid,
    )
