"""
@author: cunyue
@file: capture.py
@time: 2026/4/20
@description: stdout/stderr write 拦截 — 按需启动/停止，重入保护，fork 安全

设计要点：
- Passthrough-first: 先调用原始 write，再通知回调
- 按需启动: install() 时替换 write，uninstall() 时恢复，不影响全局
- 重入保护: ContextVar 防止同线程回调内 write 死循环，不同线程可并发
- Fork 安全: 每次 write 动态检查 PID，install 后 fork 的子进程自动跳过
- 异常安全: safe.block 捕获回调异常，写日志文件但不打终端，绝不影响用户 print()
"""

from __future__ import annotations

import contextvars
import sys
from typing import Callable, Literal, cast

from swanlab.sdk.internal.pkg import fork, safe

# ----------------------------------
# 模块级重入保护状态
# ----------------------------------

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
        on_write: Callable[[str], None],
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

    def install(self) -> None:
        """替换 sys.stdout.write 或 sys.stderr.write。"""
        if self._installed:
            return
        stream = getattr(sys, self._stream_name)
        self._original_write = stream.write
        stream.write = self._make_wrapper()  # type: ignore[method-assign]
        self._installed = True

    def uninstall(self) -> None:
        """恢复原始 write 方法。"""
        if not self._installed:
            return
        stream = getattr(sys, self._stream_name)
        stream.write = self._original_write  # type: ignore[method-assign]
        self._original_write = None
        self._installed = False

    def _make_wrapper(self):
        """创建 write() 包装器。"""
        original_write = cast(Callable[[str], int], self._original_write)
        assert original_write is not None
        on_write = self._on_write
        init_pid = self._init_pid

        def write_wrapper(data) -> int:
            # 1. Passthrough: 始终先调用原始 write
            n = original_write(data)

            # 2. Fork 子进程动态检查（install 后 fork 的子进程自动跳过）
            if fork.is_forked(init_pid):
                return n

            # 3. 重入保护：仅防同线程递归，不同线程可并发
            if _in_callback.get():
                return n
            _in_callback.set(True)

            try:
                with safe.block(message="StreamCapture on_write error", write_to_tty=False):
                    # 4. 解码 bytes → str
                    if isinstance(data, bytes):
                        text = data[:n].decode("utf-8", errors="replace")
                    else:
                        text = data[:n]
                    # 5. 调用回调
                    on_write(text)
            finally:
                _in_callback.set(False)

            return n

        return write_wrapper
