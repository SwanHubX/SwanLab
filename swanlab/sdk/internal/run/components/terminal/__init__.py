"""
@author: cunyue
@file: __init__.py
@time: 2026/4/19 20:26
@description: 终端代理组件 — 拦截 stdout/stderr，经终端模拟器处理，发射 ConsoleEvent

生命周期：
    install() → [运行中] → uninstall()

架构：
    主线程 → StreamCapture → queue → emulator 线程 → emulator.write() → _flush()
    只发射 committed 行（光标已通过 \\n 离开），跳过 \\r 覆盖的中间状态。
"""

from __future__ import annotations

import queue
import threading
from typing import Literal, Protocol

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.system.v1.console_pb2 import StreamType
from swanlab.sdk.internal.bus import ConsoleEvent, EmitterProtocol
from swanlab.sdk.internal.pkg import safe

from .capture import StreamCapture
from .emulator import TerminalEmulator

__all__ = ["TerminalProxyProtocol", "TerminalProxy"]


class TerminalProxyProtocol(Protocol):
    """终端代理协议"""

    def install(self) -> None: ...

    def uninstall(self) -> None: ...


class TerminalProxy(TerminalProxyProtocol):
    """终端代理组件。

    拦截 stdout/stderr 输出，通过终端模拟器处理 ANSI 序列和 \\r 覆盖，
    提取完整行后发射 ConsoleEvent 到事件总线。

    stdout 和 stderr 各使用独立的 TerminalEmulator，确保 stream 字段准确。
    """

    def __init__(
        self,
        emitter: EmitterProtocol,
        proxy_type: Literal["all", "stdout", "stderr", "none"],
        max_log_length: int,
        init_pid: int,
    ) -> None:
        """
        :param emitter: 事件总线发射器
        :param proxy_type: 捕获哪些流
        :param max_log_length: 每行最大字符数
        :param init_pid: 创建时 PID，用于 fork 检测
        """
        self._emitter = emitter
        self._proxy_type = proxy_type
        self._max_log_length = max_log_length
        self._init_pid = init_pid

        # 双 Emulator
        self._stdout_emulator = TerminalEmulator()
        self._stderr_emulator = TerminalEmulator()

        # 双队列
        self._stdout_queue: queue.Queue[str] = queue.Queue()
        self._stderr_queue: queue.Queue[str] = queue.Queue()

        # 线程控制
        self._stopped = threading.Event()
        self._installed = False

        # Capture 实例
        self._stdout_capture: StreamCapture | None = None
        self._stderr_capture: StreamCapture | None = None

        # 后台线程
        self._worker_thread: threading.Thread | None = None

    # ----------------------------------
    # 生命周期
    # ----------------------------------

    def install(self) -> None:
        """启动终端代理：安装 capture，启动后台线程。"""
        if self._installed or self._proxy_type == "none":
            return

        # 1. 安装 StreamCapture
        if self._proxy_type in ("all", "stdout"):
            self._stdout_capture = StreamCapture(
                stream_name="stdout",
                on_write=self._on_stdout_write,
                init_pid=self._init_pid,
            )
            self._stdout_capture.install()

        if self._proxy_type in ("all", "stderr"):
            self._stderr_capture = StreamCapture(
                stream_name="stderr",
                on_write=self._on_stderr_write,
                init_pid=self._init_pid,
            )
            self._stderr_capture.install()

        # 2. 启动工作线程（emulator 写入 + flush）
        self._worker_thread = threading.Thread(
            target=self._worker_loop,
            name="SwanLab·TerminalProxy",
            daemon=True,
        )
        self._worker_thread.start()

        self._installed = True

    def uninstall(self) -> None:
        """停止终端代理：卸载 capture，等待线程排空，最终 flush。"""
        if not self._installed:
            return

        # 1. 卸载 capture（停止拦截 write）
        if self._stdout_capture is not None:
            self._stdout_capture.uninstall()
            self._stdout_capture = None
        if self._stderr_capture is not None:
            self._stderr_capture.uninstall()
            self._stderr_capture = None

        # 2. 通知线程停止
        self._stopped.set()

        # 3. 等待工作线程排空队列
        if self._worker_thread is not None:
            self._worker_thread.join(timeout=5)
            self._worker_thread = None

        # 4. 强制提交 pending 行
        self._stdout_emulator.finalize()
        self._stderr_emulator.finalize()

        # 5. 最终 flush
        self._flush()

        self._installed = False

    # ----------------------------------
    # Write 回调（由 StreamCapture 在主线程调用）
    # ----------------------------------

    def _on_stdout_write(self, text: str, stream_type: StreamType) -> None:
        self._stdout_queue.put(text)

    def _on_stderr_write(self, text: str, stream_type: StreamType) -> None:
        self._stderr_queue.put(text)

    # ----------------------------------
    # 工作线程
    # ----------------------------------

    def _worker_loop(self) -> None:
        """后台线程：从队列读取数据，写入模拟器，立即 flush。"""
        while True:
            stdout_items = self._drain_queue(self._stdout_queue)
            if stdout_items:
                self._stdout_emulator.write("".join(stdout_items))

            stderr_items = self._drain_queue(self._stderr_queue)
            if stderr_items:
                self._stderr_emulator.write("".join(stderr_items))

            # 有数据写入就 flush
            if stdout_items or stderr_items:
                self._flush()

            # 退出条件：已停止且队列排空
            if self._stopped.is_set() and self._stdout_queue.empty() and self._stderr_queue.empty():
                return

            if not stdout_items and not stderr_items:
                self._stopped.wait(0.1)

    @staticmethod
    def _drain_queue(q: queue.Queue[str]) -> list[str]:
        """批量排空队列。"""
        items: list[str] = []
        try:
            items.append(q.get_nowait())
            while True:
                items.append(q.get_nowait())
        except queue.Empty:
            pass
        return items

    # ----------------------------------
    # Flush
    # ----------------------------------

    def _flush(self) -> None:
        """读取双模拟器 diff，发射 ConsoleEvent。

        只发射新增行（is_new_line=True），跳过 \\r 覆盖的中间状态。
        """
        with safe.block(message="Terminal proxy flush error", write_to_tty=False):
            # stdout diff
            for line_text, is_new_line in self._stdout_emulator.read():
                if not is_new_line:
                    continue
                self._emit_line(line_text, StreamType.STREAM_TYPE_STDOUT)

            # stderr diff
            for line_text, is_new_line in self._stderr_emulator.read():
                if not is_new_line:
                    continue
                self._emit_line(line_text, StreamType.STREAM_TYPE_STDERR)

    # ----------------------------------
    # ConsoleEvent 发射
    # ----------------------------------

    def _emit_line(self, line: str, stream: StreamType) -> None:
        """发射一行 ConsoleEvent，应用 max_log_length 截断。"""
        # 截断
        if len(line) > self._max_log_length:
            line = line[: self._max_log_length]

        ts = Timestamp()
        ts.GetCurrentTime()
        self._emitter.emit(ConsoleEvent(line=line, stream=stream, timestamp=ts))
