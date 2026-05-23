"""Unit tests for progress display."""

from types import SimpleNamespace
from typing import Optional

from rich.text import Text

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import GetOperationStatsResponse
from swanlab.proto.swanlab.operation.v1.operation_pb2 import (
    CoreState,
    FileProgress,
    OperationStats,
)
from swanlab.sdk.internal.run import progress as progress_module
from swanlab.sdk.internal.run.fmt import (
    fmt_bytes as _format_bytes,
)
from swanlab.sdk.internal.run.fmt import (
    fmt_items as _format_items,
)
from swanlab.sdk.internal.run.fmt import (
    fmt_rate as _format_rate,
)
from swanlab.sdk.internal.run.progress import (
    _build_renderable,
    progress_display,
    run_with_progress,
)

# ============================================================
# Format helpers
# ============================================================


class TestFormatBytes:
    def test_zero(self):
        assert _format_bytes(0) == "0 B"

    def test_bytes(self):
        assert _format_bytes(42) == "42 B"

    def test_kilobytes(self):
        assert _format_bytes(1024) == "1.0 KB"

    def test_megabytes(self):
        assert _format_bytes(2 * 1024 * 1024) == "2.0 MB"

    def test_gigabytes(self):
        assert _format_bytes(1536 * 1024 * 1024) == "1.5 GB"


class TestFormatRate:
    def test_bytes_rate(self):
        assert _format_rate(1024.0, "bytes") == "1.0 KB/s"

    def test_records_rate(self):
        assert _format_rate(32.0, "records") == "32.0 records/s"

    def test_zero_rate(self):
        assert _format_rate(0.0, "records") == "0.0 records/s"


class TestFormatItems:
    def test_bytes(self):
        assert _format_items(1024, "bytes") == "1.0 KB"

    def test_records(self):
        assert _format_items(120, "records") == "120 records"


# ============================================================
# Build renderable
# ============================================================


def _make_stats(
    total_number: int = 1,
    uploaded_number: int = 0,
    total_size: int = 200 * 1024 * 1024,
    uploaded_size: int = 100 * 1024 * 1024,
    rate: float = 2.0 * 1024 * 1024,
    state=CoreState.CORE_STATE_RUNNING,
    files: Optional[list] = None,
    total_records: int = 0,
    uploaded_records: int = 0,
) -> OperationStats:
    return OperationStats(
        state=state,  # pyright: ignore[reportArgumentType]
        total_number=total_number,
        uploaded_number=uploaded_number,
        total_size=total_size,
        uploaded_size=uploaded_size,
        rate=rate,
        files=files or [],
        total_records=total_records,
        uploaded_records=uploaded_records,
    )


def _render_plain(stats: OperationStats, unit: str = "bytes") -> str:
    """Build renderable and return its plain text for assertion."""
    from rich.console import Console

    c = Console(width=120, force_terminal=True, no_color=True, legacy_windows=False)
    r = _build_renderable(stats, unit)
    with c.capture() as cap:
        c.print(r, end="")
    return cap.get()


class TestBuildRenderable:
    def test_bytes_with_total(self):
        stats = _make_stats()
        text = _render_plain(stats, "bytes")
        assert "100.0 MB" in text
        assert "200.0 MB" in text

    def test_bytes_no_total(self):
        stats = _make_stats(total_number=0, total_size=0, uploaded_size=50 * 1024 * 1024, rate=0)
        text = _render_plain(stats, "bytes")
        assert "50.0 MB" in text

    def test_records_with_total_shows_bar(self):
        stats = _make_stats(total_records=200, uploaded_records=120, rate=32.0)
        text = _render_plain(stats, "records")
        assert "120 records / 200 records" in text
        assert "█" in text

    def test_records_no_total_no_bar(self):
        stats = _make_stats(total_records=0, uploaded_records=120, rate=32.0)
        text = _render_plain(stats, "records")
        assert "Uploaded 120 records" in text
        assert "█" not in text

    def test_auto_files_take_precedence_over_records(self):
        stats = _make_stats(total_number=1, total_size=100, uploaded_size=50, total_records=200, uploaded_records=120)
        text = _render_plain(stats, "auto")
        assert "50 B" in text
        assert "100 B" in text
        assert "records" not in text

    def test_uploading_files_shown(self):
        fp = FileProgress(
            path="/path/to/checkpoint.pt",
            total=200 * 1024 * 1024,
            uploaded=80 * 1024 * 1024,
            rate=2.0 * 1024 * 1024,
        )
        stats = _make_stats(total_number=1, files=[fp])
        text = _render_plain(stats, "bytes")
        assert "/path/to/checkpoint.pt" in text
        assert "80.0 MB" in text

    def test_uploading_file_count_uses_completed_files(self):
        fp = FileProgress(
            path="/path/to/checkpoint.pt",
            total=200 * 1024 * 1024,
            uploaded=80 * 1024 * 1024,
            rate=2.0 * 1024 * 1024,
        )
        stats = _make_stats(total_number=3, uploaded_number=1, files=[fp])

        text = _render_plain(stats, "bytes")

        assert "(1/3):" in text

    def test_uploading_files_rows_do_not_have_spacer(self):
        fp = FileProgress(
            path="/path/to/checkpoint.pt",
            total=200 * 1024 * 1024,
            uploaded=80 * 1024 * 1024,
            rate=2.0 * 1024 * 1024,
        )
        stats = _make_stats(total_number=1, files=[fp])

        text = _render_plain(stats, "bytes")

        assert all(line.strip() for line in text.splitlines())

    def test_progress_bars_use_default_terminal_colors(self, monkeypatch):
        calls = []

        def fake_bar(*args, **kwargs):
            calls.append(kwargs)
            return Text("bar")

        fp = FileProgress(
            path="/path/to/checkpoint.pt",
            total=200 * 1024 * 1024,
            uploaded=80 * 1024 * 1024,
            rate=2.0 * 1024 * 1024,
        )
        stats = _make_stats(total_number=1, files=[fp])
        monkeypatch.setattr(progress_module, "Bar", fake_bar)

        _build_renderable(stats, "bytes")

        assert calls
        assert all("color" not in kwargs and "bgcolor" not in kwargs for kwargs in calls)

    def test_uploading_files_rate_is_labeled_as_average_with_dynamic_units(self):
        kb_file = FileProgress(
            path="/path/to/small.txt",
            total=200 * 1024,
            uploaded=80 * 1024,
            rate=512.0 * 1024,
        )
        mb_file = FileProgress(
            path="/path/to/checkpoint.pt",
            total=200 * 1024 * 1024,
            uploaded=80 * 1024 * 1024,
            rate=2.0 * 1024 * 1024,
        )
        stats = _make_stats(total_number=2, files=[kb_file, mb_file], rate=2.5 * 1024 * 1024)

        text = _render_plain(stats, "bytes")

        assert "avg 2.5 MB/s" in text
        assert "avg 512.0 KB/s" in text
        assert "avg 2.0 MB/s" in text

    def test_uploading_files_limit(self):
        files = [FileProgress(path=f"/file{i}.pt", total=100, uploaded=0, rate=0) for i in range(10)]
        stats = _make_stats(total_number=10, files=files)
        text = _render_plain(stats, "bytes")
        # default limit is 5
        assert "/file4.pt" in text
        assert "/file5.pt" not in text


# ============================================================
# progress_display context manager
# ============================================================


class TestProgressDisplay:
    def test_context_manager_closes_live(self, monkeypatch):
        closed = []

        class FakeLive:
            def __init__(self, *a, **kw):
                pass

            def start(self):
                pass

            def update(self, *a, **kw):
                pass

            def stop(self):
                closed.append(True)

        monkeypatch.setattr(progress_module, "Live", FakeLive)
        monkeypatch.setattr(progress_module.console, "c", SimpleNamespace(is_terminal=True))
        monkeypatch.setattr(progress_module.os, "isatty", lambda _: True)

        with progress_display() as display:
            display.update(_make_stats(total_number=1, total_size=100, uploaded_size=50))

        assert closed  # stop() was called


# ============================================================
# run_with_progress integration
# ============================================================


class _FakeProvider:
    """Stats provider that returns pre-set stats and optionally transitions to FINISHED."""

    def __init__(self, stats_sequence: Optional[list] = None):
        self._iter = iter(stats_sequence or [])

    def get_operation_stats(self) -> GetOperationStatsResponse:
        try:
            return GetOperationStatsResponse(success=True, message="OK", stats=next(self._iter))
        except StopIteration:
            return GetOperationStatsResponse(
                success=True,
                message="OK",
                stats=OperationStats(state=CoreState.CORE_STATE_FINISHED),
            )


class TestRunWithProgress:
    def test_immediate_finished(self):
        provider = _FakeProvider(
            [
                OperationStats(state=CoreState.CORE_STATE_FINISHED, total_number=0, uploaded_number=0, rate=0),
            ]
        )
        called = []
        run_with_progress(provider.get_operation_stats, lambda: called.append(True), timeout=5.0)
        assert called == [True]

    def test_timeout(self):
        running_stats = OperationStats(
            state=CoreState.CORE_STATE_RUNNING,
            total_number=100,
            uploaded_number=10,
            rate=1.0,
            total_size=100,
            uploaded_size=10,
        )

        class AlwaysRunningProvider:
            def get_operation_stats(self) -> GetOperationStatsResponse:
                return GetOperationStatsResponse(success=True, message="OK", stats=running_stats)

        called = []
        run_with_progress(
            AlwaysRunningProvider().get_operation_stats,
            lambda: called.append(True),
            timeout=0.3,
            unit="bytes",
        )
        assert called == [True]

    def test_transitions_to_finished(self):
        provider = _FakeProvider(
            [
                OperationStats(
                    state=CoreState.CORE_STATE_RUNNING,
                    total_number=100,
                    uploaded_number=50,
                    rate=10.0,
                    total_records=100,
                    uploaded_records=50,
                ),
                OperationStats(
                    state=CoreState.CORE_STATE_RUNNING,
                    total_number=100,
                    uploaded_number=80,
                    rate=10.0,
                    total_records=100,
                    uploaded_records=80,
                ),
                OperationStats(
                    state=CoreState.CORE_STATE_FINISHED,
                    total_number=100,
                    uploaded_number=100,
                    rate=0,
                    total_records=100,
                    uploaded_records=100,
                ),
            ]
        )
        called = []
        run_with_progress(provider.get_operation_stats, lambda: called.append(True), timeout=5.0, unit="records")
        assert called == [True]

    def test_blocking_fn_receives_return_value(self):
        provider = _FakeProvider(
            [
                OperationStats(state=CoreState.CORE_STATE_FINISHED),
            ]
        )

        class CustomResult:
            pass

        expected = CustomResult()
        result = run_with_progress(provider.get_operation_stats, lambda: expected, timeout=5.0)
        assert result is expected

    def test_live_updates_refresh_immediately(self, monkeypatch):
        class FakeLive:
            instances = []

            def __init__(self, *args, **kwargs):
                self.updates = []
                FakeLive.instances.append(self)

            def start(self):
                pass

            def update(self, renderable, *, refresh=False):
                self.updates.append(refresh)

            def stop(self):
                pass

        provider = _FakeProvider(
            [
                # 第一次：初始检查，发现需要显示进度
                OperationStats(
                    state=CoreState.CORE_STATE_RUNNING,
                    total_number=1,
                    uploaded_number=1,
                    total_size=100,
                    uploaded_size=100,
                ),
                # 第二次：此时立刻被设置为 CORE_STATE_FINISHED ，刷新一次，退出循环
                OperationStats(
                    state=CoreState.CORE_STATE_FINISHED,
                    total_number=1,
                    uploaded_number=1,
                    total_size=100,
                    uploaded_size=100,
                ),
                # 第三次：获取最终状态，确认显示最终状态时不需要刷新了
                OperationStats(
                    state=CoreState.CORE_STATE_FINISHED,
                    total_number=1,
                    uploaded_number=1,
                    total_size=100,
                    uploaded_size=100,
                ),
            ]
        )
        monkeypatch.setattr(progress_module, "Live", FakeLive)
        monkeypatch.setattr(progress_module.console, "c", SimpleNamespace(is_terminal=True))
        monkeypatch.setattr(progress_module.os, "isatty", lambda _: True)
        monkeypatch.setattr(progress_module, "_print_summary", lambda *args: None)

        result = run_with_progress(provider.get_operation_stats, lambda: "ok", timeout=5.0)

        assert result == "ok"
        assert FakeLive.instances[0].updates == [True]
