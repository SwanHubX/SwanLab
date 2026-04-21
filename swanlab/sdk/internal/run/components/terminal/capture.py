"""
@author: cunyue
@file: capture.py
@time: 2026/4/20
@description: stdout/stderr write 拦截 — 按需启动/停止，重入保护，fork 安全

设计要点：
- Passthrough-first: 先调用原始 write，再通知回调
- 按需启动: install() 时替换 write，uninstall() 时恢复，不影响全局
- 重入保护: RLock + _is_writing + ContextVar 防止回调内 write 死循环
- Fork 安全: PID guard 静默跳过 + fork.register() 重置锁
- 异常安全: 回调异常静默吞掉，绝不让 callback 错误传播到用户 print()
"""

from __future__ import annotations

import contextvars
import sys
import threading
from typing import Callable, Literal

from swanlab.proto.swanlab.system.v1.console_pb2 import StreamType
from swanlab.sdk.internal.pkg import fork

# ----------------------------------
# 模块级重入保护状态
# ----------------------------------

_module_rlock = threading.RLock()
_is_writing: bool = False
_in_callback: contextvars.ContextVar[bool] = contextvars.ContextVar(
    "_swanlab_console_in_callback",
    default=False,
)


class StreamCapture:
    """拦截 sys.stdout 或 sys.stderr 的 write 方法。

    install() 时替换流对象的 write 方法为包装器：
    1. 调用原始 write (passthrough)
    2. 调用 on_write 回调

    uninstall() 时恢复原始 write 方法。
    """

    def __init__(
        self,
        stream_name: Literal["stdout", "stderr"],
        on_write: Callable[[str, StreamType], None],
        init_pid: int,
    ) -> None:
        """
        :param stream_name: "stdout" 或 "stderr"
        :param on_write: 回调函数 callback(written_text, stream_type)
        :param init_pid: 构造时的 PID，用于 fork 检测
        """
        self._stream_name = stream_name
        self._on_write = on_write
        self._init_pid = init_pid
        self._original_write: Callable | None = None
        self._installed = False
        self._fork_reset_callback: Callable[[], None] | None = None

        # 预计算 stream_type
        self._stream_type = StreamType.STREAM_TYPE_STDOUT if stream_name == "stdout" else StreamType.STREAM_TYPE_STDERR

    def install(self) -> None:
        """替换 sys.stdout.write 或 sys.stderr.write。"""
        if self._installed:
            return
        stream = getattr(sys, self._stream_name)
        self._original_write = stream.write
        stream.write = self._make_wrapper()  # type: ignore[method-assign]
        self._installed = True
        # 注册 fork 清理回调
        self._fork_reset_callback = _reset_locks_in_child
        fork.register(self._fork_reset_callback)

    def uninstall(self) -> None:
        """恢复原始 write 方法。"""
        if not self._installed:
            return
        stream = getattr(sys, self._stream_name)
        stream.write = self._original_write  # type: ignore[method-assign]
        self._original_write = None
        self._installed = False
        # 注销 fork 回调
        if self._fork_reset_callback is not None:
            fork.unregister(self._fork_reset_callback)
            self._fork_reset_callback = None

    def _make_wrapper(self):
        """创建 write() 包装器。"""
        original_write = self._original_write
        on_write = self._on_write
        init_pid = self._init_pid
        stream_type = self._stream_type
        assert original_write is not None

        def write_wrapper(data) -> int:
            # 1. Passthrough: 始终先调用原始 write
            n = original_write(data)

            # 2. Fork 安全: fork 子进程静默跳过
            if fork.is_forked(init_pid):
                return n

            # 3. 重入保护
            global _is_writing
            with _module_rlock:
                if _is_writing or _in_callback.get():
                    return n
                _is_writing = True
                _in_callback.set(True)

            try:
                # 4. 解码 bytes → str
                if isinstance(data, bytes):
                    text = data[:n].decode("utf-8", errors="replace")
                else:
                    text = data[:n]
                # 5. 调用回调
                on_write(text, stream_type)
            except Exception:
                # 静默吞掉异常，绝不影响用户的 print()
                pass
            finally:
                with _module_rlock:
                    _is_writing = False
                    _in_callback.set(False)

            return n

        return write_wrapper


def _reset_locks_in_child() -> None:
    """Fork 回调：在子进程中重置锁状态。"""
    global _module_rlock, _is_writing
    _module_rlock = threading.RLock()
    _is_writing = False
