"""Upload progress display for run finish.
@file: progress.py
@description: Rich Live poll/render progress printer for upload progress display.
"""

from __future__ import annotations

import contextlib
import os
import time
from collections.abc import Callable, Iterator
from typing import Optional, TypeVar

from rich.bar import Bar
from rich.live import Live
from rich.spinner import Spinner
from rich.table import Table
from rich.text import Text

from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState, OperationStats
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.run.fmt import fmt_bytes, fmt_items, fmt_rate

_T = TypeVar("_T")

_BAR_WIDTH = 24
_MAX_FILES = 5
_POLL_INTERVAL = 0.25


# ── Unit helpers ──────────────────────────────────────────────


def _resolve_unit(stats: OperationStats, unit: str) -> str:
    """自动推断显示单位：有文件上传时用 bytes，否则用 records。"""
    if unit != "auto":
        return unit
    if stats.total_number > 0 or stats.total_size > 0:
        return "bytes"
    if stats.total_records > 0 or stats.uploaded_records > 0:
        return "records"
    return "records"


def _resolve_progress(stats: OperationStats, unit: str) -> tuple[int, int]:
    """Return (total, done) for an already resolved unit."""
    if unit == "bytes":
        return stats.total_size, stats.uploaded_size
    return stats.total_records, stats.uploaded_records


def _should_show(stats: OperationStats, unit: str) -> bool:
    unit = _resolve_unit(stats, unit)
    total, done = _resolve_progress(stats, unit)
    return total > 0 or done > 0


# ── Rich renderable ──────────────────────────────────────────


def _build_renderable(stats: OperationStats, unit: str) -> Table:
    unit = _resolve_unit(stats, unit)
    table = Table.grid(padding=(0, 2))
    if unit == "bytes" and (stats.total_number > 0 or stats.total_size > 0):
        _render_files(table, stats)
    else:
        _render_records(table, stats, unit)
    return table


def _render_files(table: Table, stats: OperationStats) -> None:
    pct = int(stats.uploaded_size / max(stats.total_size, 1) * 100)

    left = Spinner(
        "dots",
        Text.assemble(
            ("Syncing files ", "bold"),
            (f"({stats.uploaded_number}/{stats.total_number}):", ""),
        ),
    )
    bar = Bar(size=max(stats.total_size, 1), begin=0, end=stats.uploaded_size, width=_BAR_WIDTH)
    right = Text.assemble(
        (f"{fmt_bytes(stats.uploaded_size)} / {fmt_bytes(stats.total_size)} ", "bold"),
        (f"({pct}%)", "bold"),
    )
    _append_rate(right, stats.rate)
    table.add_row(left, bar, right)

    for fp in stats.files[:_MAX_FILES]:
        fpct = int(fp.uploaded / max(fp.total, 1) * 100)
        left = Spinner("dots", Text.assemble(("  Uploading ", "dim"), (fp.path, "bold")))
        bar = Bar(size=max(fp.total, 1), begin=0, end=fp.uploaded, width=12)
        right = Text.assemble(
            (f"{fmt_bytes(fp.uploaded)} / {fmt_bytes(fp.total)} ", ""),
            (f"({fpct}%)", ""),
        )
        _append_rate(right, fp.rate)
        table.add_row(left, bar, right)


def _render_records(table: Table, stats: OperationStats, unit: str) -> None:
    total, done = _resolve_progress(stats, unit)
    if total > 0:
        pct = int(done / total * 100)
        left = Spinner("dots", Text(f"Uploading {unit}: ", "bold"))
        bar = Bar(size=total, begin=0, end=done, width=_BAR_WIDTH)
        right = Text.assemble(
            (f"{fmt_items(done, unit)} / {fmt_items(total, unit)} ", "bold"),
            (f"({pct}%)", "bold"),
        )
        table.add_row(left, bar, right)
    else:
        left = Spinner("dots", Text(f"Uploaded {fmt_items(done, unit)}", "bold"))
        table.add_row(left, Text(""), Text(""))


def _append_rate(text: Text, rate: float) -> None:
    if rate > 0:
        text.append(Text(f", avg {fmt_rate(rate, 'bytes')}", "dim"))


# ── Summary printer ──────────────────────────────────────────


def _print_summary(stats: OperationStats, unit: str, finished: bool) -> None:
    unit = _resolve_unit(stats, unit)
    is_files = unit == "bytes" and (stats.total_number > 0 or stats.total_size > 0)

    if is_files:
        n, nt = stats.uploaded_number, stats.total_number
        done, total = stats.uploaded_size, stats.total_size
        if finished and done >= total:
            msg = f"Upload complete: {n} file(s) ({fmt_bytes(done)})"
        elif not finished:
            msg = f"Upload timeout: {n}/{nt} files ({fmt_bytes(done)}/{fmt_bytes(total)})"
        else:
            msg = f"Upload finished: {n}/{nt} files ({fmt_bytes(done)}/{fmt_bytes(total)})"
    else:
        total, done = _resolve_progress(stats, unit)
        if finished and done >= total:
            msg = f"Upload complete: {fmt_items(done, unit)}"
        elif not finished:
            msg = f"Upload timeout: {fmt_items(done, unit)}/{fmt_items(total, unit)}"
        else:
            msg = f"Upload finished: {fmt_items(done, unit)}/{fmt_items(total, unit)}"

    console.info(msg)


# ── Live display ─────────────────────────────────────────────


@contextlib.contextmanager
def progress_display() -> Iterator[ProgressDisplay]:
    """Context manager providing a Rich Live progress area for upload stats."""
    display = ProgressDisplay()
    try:
        yield display
    finally:
        display.close()


class ProgressDisplay:
    """管理 Rich Live 进度渲染区域。

    仅在 TTY 终端激活。summary 输出也收敛在 display 内部，
    避免调用方依赖渲染状态字段。
    """

    def __init__(self) -> None:
        self._is_tty = console.c.is_terminal and hasattr(os, "isatty") and os.isatty(1)
        self._live: Optional[Live] = None
        self._has_shown = False

    def update(self, stats: OperationStats, unit: str = "auto") -> None:
        if not _should_show(stats, unit):
            return
        self._has_shown = True
        if self._is_tty and self._live is None:
            self._live = Live(
                console=console.c,
                transient=True,
                refresh_per_second=4,
                vertical_overflow="visible",
            )
            self._live.start()
        if self._live is not None:
            self._live.update(_build_renderable(stats, unit), refresh=True)

    def close(self) -> None:
        if self._live is not None:
            self._live.stop()

    def print_summary(self, stats: OperationStats, unit: str, finished: bool) -> None:
        if self._has_shown:
            _print_summary(stats, unit, finished)


# ── Progress waiter ──────────────────────────────────────────


def run_with_progress(
    stats_fn: Callable[[], OperationStats],
    blocking_fn: Callable[[], _T],
    timeout: Optional[float] = None,
    unit: str = "auto",
) -> _T:
    """Poll *stats_fn* while waiting for upload completion, then call *blocking_fn*.

    Displays a Rich Live progress area that updates every poll interval.
    The loop exits when the tracker reports ``CORE_STATE_FINISHED`` or
    when *timeout* elapses.  *blocking_fn* is called once after the
    display closes (e.g. to join the transport and confirm with backend).
    """
    deadline = time.monotonic() + timeout if timeout is not None else None
    finished = False

    with progress_display() as display:
        while True:
            stats = stats_fn()
            display.update(stats, unit)
            if stats.state == CoreState.CORE_STATE_FINISHED:
                finished = True
                break
            if deadline is not None and time.monotonic() >= deadline:
                break
            time.sleep(_POLL_INTERVAL)

    display.print_summary(stats_fn(), unit, finished)

    return blocking_fn()
