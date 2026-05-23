"""
@file: tracker.py
@description: Thread-safe upload progress tracker.
"""

from __future__ import annotations

import threading
import time
from dataclasses import dataclass, field
from typing import Optional

from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState, FileProgress, OperationStats


@dataclass
class _RateTracker:
    """EMA 速率追踪器：在采样窗口内累积 delta，到达窗口边界时计算瞬时速率并混合。"""

    rate: float = 0.0
    last_advance_time: Optional[float] = None
    pending_delta: int = 0

    # EMA 平滑因子，值越大速率对近期吞吐越敏感
    _ALPHA: float = 0.3
    # 采样窗口（秒），窗口内的 delta 累积后一次性刷入，避免逐字节抖动
    _SAMPLE_INTERVAL: float = 0.5

    def advance(self, delta: int) -> None:
        """累积 delta 并在采样窗口边界时更新 EMA 速率。"""
        if delta < 0:
            return
        now = time.monotonic()
        self.pending_delta += delta

        if self.last_advance_time is None:
            self.last_advance_time = now
            return
        elapsed = now - self.last_advance_time
        if elapsed < self._SAMPLE_INTERVAL:
            return
        instant_rate = self.pending_delta / elapsed
        if self.rate <= 0:
            self.rate = instant_rate
        else:
            self.rate = self._ALPHA * instant_rate + (1 - self._ALPHA) * self.rate
        self.last_advance_time = now
        self.pending_delta = 0

    def decayed_rate(self) -> float:
        """返回经指数衰减后的当前速率。"""
        if self.rate <= 0 or self.last_advance_time is None:
            return 0.0
        now = time.monotonic()
        elapsed = now - self.last_advance_time
        if elapsed < self._SAMPLE_INTERVAL:
            return self.rate
        intervals = elapsed / self._SAMPLE_INTERVAL
        decayed = self.rate * ((1 - self._ALPHA) ** intervals)
        return decayed


@dataclass
class _TrackedFileProgress:
    progress: FileProgress
    sequence: int = 0
    rate_tracker: _RateTracker = field(default_factory=_RateTracker)
    parts: dict[str, int] = field(default_factory=dict)


class UploadTracker:
    """Thread-safe aggregate progress tracker with optional file display rows."""

    _RATE_SAMPLE_INTERVAL = _RateTracker._SAMPLE_INTERVAL

    def __init__(self, uploading_files_limit: int = 5) -> None:
        self._lock = threading.Lock()
        self._uploading_files_limit = max(uploading_files_limit, 0)
        self._state: CoreState = CoreState.CORE_STATE_NOT_STARTED
        self._total_number: int = 0
        self._uploaded_number: int = 0
        self._total_size: int = 0
        self._uploaded_size: int = 0
        self._total_records: int = 0
        self._uploaded_records: int = 0
        self._global_rate = _RateTracker()
        self._files: dict[str, _TrackedFileProgress] = {}
        self._registered_keys: set[str] = set()
        self._file_sequence: int = 0

    def set_state(self, state: CoreState) -> None:
        with self._lock:
            self._state = state

    def add_total(self, n: int) -> None:
        amount = self._normalize_amount(n)
        if amount == 0:
            return
        with self._lock:
            self._total_number += amount

    def advance(self, delta: int) -> None:
        amount = self._normalize_amount(delta)
        if amount == 0:
            return
        with self._lock:
            self._uploaded_number += amount

    def set_total_records(self, total: int) -> None:
        """Set the known final record total. A zero value means the total is unknown."""
        amount = self._normalize_amount(total)
        with self._lock:
            self._total_records = amount

    def add_total_records(self, delta: int) -> None:
        amount = self._normalize_amount(delta)
        if amount == 0:
            return
        with self._lock:
            self._total_records += amount

    def advance_records(self, delta: int) -> None:
        amount = self._normalize_amount(delta)
        if amount == 0:
            return
        with self._lock:
            self._uploaded_records += amount

    def add_file(self, key: str, path: str, total: int) -> None:
        amount = self._normalize_amount(total)
        if amount == 0:
            return
        with self._lock:
            if key in self._registered_keys:
                return
            now = time.monotonic()
            self._registered_keys.add(key)
            self._file_sequence += 1
            self._files[key] = _TrackedFileProgress(
                progress=FileProgress(path=path, total=amount),
                sequence=self._file_sequence,
                rate_tracker=_RateTracker(last_advance_time=now),
            )
            self._total_number += 1
            self._total_size += amount
            if self._global_rate.last_advance_time is None:
                self._global_rate.last_advance_time = now

    def update_file_progress(self, file_key: str, upload_key: str, current_bytes: int) -> None:
        current_bytes = self._normalize_amount(current_bytes)
        with self._lock:
            file_progress = self._files.get(file_key)
            if file_progress is None:
                return

            # Compute delta for this specific upload/part key
            prev_bytes = file_progress.parts.get(upload_key, 0)
            delta = current_bytes - prev_bytes
            if delta == 0:
                return
            if delta < 0 and current_bytes > 0:
                return

            # Update part progress
            file_progress.parts[upload_key] = current_bytes

            # Update file's total uploaded bytes
            file_progress.progress.uploaded += delta

            # Update rate for this file
            file_progress.rate_tracker.advance(delta)

            # Update global uploaded bytes
            self._uploaded_size += delta

            # Update global rate
            self._global_rate.advance(delta)

    def finish_file(self, key: str) -> None:
        with self._lock:
            file_progress = self._files.pop(key, None)
            if file_progress is not None:
                remaining = file_progress.progress.total - file_progress.progress.uploaded
                if remaining > 0:
                    # 将已完成/跳过的文件视为全部上传，使总进度条可达 100%，
                    # 即使后端接受了部分上传（如去重跳过）也能正确收敛。
                    self._uploaded_size += remaining
                self._uploaded_number += 1

    def snapshot(self) -> OperationStats:
        with self._lock:
            stats = OperationStats(
                state=self._state,
                total_number=self._total_number,
                uploaded_number=self._uploaded_number,
                total_size=self._total_size,
                uploaded_size=self._uploaded_size,
                rate=self._global_rate.decayed_rate(),
                total_records=self._total_records,
                uploaded_records=self._uploaded_records,
            )
            for file_progress in self._visible_files():
                prog = FileProgress()
                prog.CopyFrom(file_progress.progress)
                prog.rate = file_progress.rate_tracker.decayed_rate()
                stats.files.add().CopyFrom(prog)
            return stats

    def _visible_files(self) -> list[_TrackedFileProgress]:
        """返回前 N 个未完成文件，已开始上传的优先显示。"""
        unfinished = [file for file in self._files.values() if file.progress.uploaded < file.progress.total]
        unfinished.sort(key=lambda file: (file.progress.uploaded <= 0, file.sequence))
        return unfinished[: self._uploading_files_limit]

    @staticmethod
    def _normalize_amount(amount: int) -> int:
        if amount <= 0:
            return 0
        return int(amount)


__all__ = ["UploadTracker"]
