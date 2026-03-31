"""上传进度追踪器，在 `swanlab.save()` 上传文件时展示 rich spinner 动画。"""

import threading

from rich.status import Status


class _SaveProgress:
    """后台上传 spinner，实时展示文件上传计数进度
    线程安全，支持多次 ``save()`` 调用累加总量，全部完成后自动消失。
    """

    def __init__(self, total: int) -> None:
        self._total: int = 0
        self._completed: int = 0
        self._lock = threading.Lock()
        self._status = Status("", spinner="dots")
        self._add_unlocked(total)

    # ------------------------------------------------------------------
    # 公开接口
    # ------------------------------------------------------------------

    def add(self, n: int) -> None:
        """增加待上传文件数（多次 ``save()`` 调用时累加）。"""
        with self._lock:
            self._add_unlocked(n)

    def done(self) -> None:
        """标记一个文件上传完成。当所有文件完成时自动停止 spinner。"""
        with self._lock:
            self._completed += 1
            self._refresh_display()
            if self._completed >= self._total and self._total > 0:
                self._status.stop()

    def stop(self) -> None:
        """强制停止 spinner（在 ``close()`` 时调用）。"""
        self._status.stop()

    # ------------------------------------------------------------------
    # 内部方法
    # ------------------------------------------------------------------

    def _add_unlocked(self, n: int) -> None:
        """必须在持有 ``_lock`` 时调用。"""
        self._total += n
        self._refresh_display()
        if self._total > self._completed and not self._status._live.is_started:
            self._status.start()
        # else: 所有文件已完成，无需启动

    def _refresh_display(self) -> None:
        self._status.update(f"Uploading files ({self._completed}/{self._total})...")
