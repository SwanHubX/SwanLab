import threading
import time
from concurrent.futures import thread as futures_thread
from pathlib import Path
from types import SimpleNamespace
from typing import cast
from unittest.mock import MagicMock

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.probe_python.context import ProbeConfig, ProbeContext
from swanlab.sdk.internal.probe_python.monitor import Monitor
from swanlab.sdk.internal.probe_python.protocol import CollectorProtocol
from swanlab.sdk.internal.probe_python.typings import SystemScalar, SystemShim
from swanlab.sdk.protocol import CoreProtocol


class _Collector(CollectorProtocol):
    def collect(self):
        return [("cpu.util", 0.5)]


def _start_monitor(monkeypatch, tmp_path, collector):
    """挂载单个采集器启动 Monitor，返回 (monitor, trace 替身, upsert_scalars 替身)。"""
    shim = SystemShim(slug="linux", enable_cpu=True)
    scalar = SystemScalar(key="cpu.util", chart_name="CPU")
    upsert_scalars = MagicMock()
    core = cast(CoreProtocol, SimpleNamespace(upsert_columns=MagicMock(), upsert_scalars=upsert_scalars))
    ctx = ProbeContext(
        config=ProbeConfig(
            run_id="run-id",
            run_dir=tmp_path,
            global_system_step=0,
            hardware=True,
            runtime=False,
            requirements=False,
            conda=False,
            git=False,
            swanlab=False,
            monitor=True,
            monitor_interval=3600,
            monitor_disk_dir=Path(tmp_path),
        )
    )

    monkeypatch.setattr("swanlab.sdk.internal.probe_python.monitor.CPU.new", lambda _shim: (collector, [scalar]))
    monkeypatch.setattr("swanlab.sdk.internal.probe_python.monitor.Memory.new", lambda _shim: None)
    monkeypatch.setattr("swanlab.sdk.internal.probe_python.monitor.Apple.new", lambda _shim: None)
    monkeypatch.setattr("swanlab.sdk.internal.pkg.timer.Timer.start", lambda self: self)

    trace_mock = MagicMock()
    monkeypatch.setattr(console, "trace", trace_mock)

    monitor = Monitor(shim, core).start(ctx)
    assert monitor is not None
    return monitor, trace_mock, upsert_scalars


def test_monitor_collects_without_executor_error_during_interpreter_shutdown(monkeypatch, tmp_path):
    shim = SystemShim(slug="linux", enable_cpu=True)
    monitor, trace_mock, upsert_scalars = _start_monitor(monkeypatch, tmp_path, _Collector(shim))

    monkeypatch.setattr(futures_thread, "_shutdown", True)
    monitor._timer._execute_once()  # type: ignore[union-attr]
    monkeypatch.setattr(futures_thread, "_shutdown", False)
    monitor.stop()

    trace_mock.assert_not_called()
    upsert_scalars.assert_called_once()


def test_monitor_collect_timeout_degrades(monkeypatch, tmp_path):
    """采集任务超过 COLLECT_TIMEOUT 时按失败处理，不阻塞也不产生标量。"""
    release = threading.Event()

    class _SlowCollector(CollectorProtocol):
        def collect(self):
            release.wait(5)
            return [("cpu.util", 0.5)]

    monkeypatch.setattr("swanlab.sdk.internal.probe_python.monitor.COLLECT_TIMEOUT", 0.05)
    shim = SystemShim(slug="linux", enable_cpu=True)
    monitor, trace_mock, upsert_scalars = _start_monitor(monkeypatch, tmp_path, _SlowCollector(shim))

    try:
        monitor._timer._execute_once()  # type: ignore[union-attr]
    finally:
        release.set()
        monitor.stop()

    trace_mock.assert_called()
    assert upsert_scalars.call_args[0][0] == []


def test_monitor_stop_is_bounded_when_collector_hangs(monkeypatch, tmp_path):
    """采集任务卡死时 stop() 不能被无限拖住。"""
    release = threading.Event()
    shim = SystemShim(slug="linux", enable_cpu=True)
    monitor, _, _ = _start_monitor(monkeypatch, tmp_path, _Collector(shim))

    executor = monitor._executor
    assert executor is not None
    executor.submit(release.wait, 5)  # 模拟卡在外部命令上的采集任务
    monkeypatch.setattr("swanlab.sdk.internal.probe_python.monitor.STOP_TIMEOUT", 0.05)

    started = time.monotonic()
    try:
        monitor.stop()
    finally:
        release.set()

    assert time.monotonic() - started < 1


def test_monitor_does_not_upsert_after_stop_timeout(monkeypatch, tmp_path):
    """collector 卡住导致 stop() 超时返回后，collector 再返回也不会继续调 Core。"""
    started = threading.Event()
    release = threading.Event()

    class _BlockedCollector(CollectorProtocol):
        def collect(self):
            started.set()
            release.wait(5)
            return [("cpu.util", 0.5)]

    shim = SystemShim(slug="linux", enable_cpu=True)
    monitor, _, upsert_scalars = _start_monitor(monkeypatch, tmp_path, _BlockedCollector(shim))
    monkeypatch.setattr("swanlab.sdk.internal.probe_python.monitor.STOP_TIMEOUT", 0.05)

    assert monitor._timer is not None
    worker = threading.Thread(target=monitor._timer._execute_once, daemon=True)
    worker.start()
    try:
        assert started.wait(2), "collector did not start in time"
        stopped_at = time.monotonic()
        monitor.stop()
        assert time.monotonic() - stopped_at < 1
    finally:
        release.set()
        worker.join(2)

    upsert_scalars.assert_not_called()


def test_monitor_collect_timeout_uses_shared_round_deadline(monkeypatch, tmp_path):
    """单轮采集超时对所有采集器共享同一 deadline，而不是每个采集器各等一轮。"""
    release = threading.Event()

    class _SlowCollector(CollectorProtocol):
        def collect(self):
            release.wait(5)
            return [("cpu.util", 0.5)]

    shim = SystemShim(slug="linux", enable_cpu=True, enable_memory=True)
    cpu_scalar = SystemScalar(key="cpu.util", chart_name="CPU")
    mem_scalar = SystemScalar(key="mem.used", chart_name="Mem")
    core = cast(CoreProtocol, SimpleNamespace(upsert_columns=MagicMock(), upsert_scalars=MagicMock()))
    ctx = ProbeContext(
        config=ProbeConfig(
            run_id="run-id",
            run_dir=tmp_path,
            global_system_step=0,
            hardware=True,
            runtime=False,
            requirements=False,
            conda=False,
            git=False,
            swanlab=False,
            monitor=True,
            monitor_interval=3600,
            monitor_disk_dir=Path(tmp_path),
        )
    )

    monkeypatch.setattr(
        "swanlab.sdk.internal.probe_python.monitor.CPU.new", lambda s: (_SlowCollector(s), [cpu_scalar])
    )
    monkeypatch.setattr(
        "swanlab.sdk.internal.probe_python.monitor.Memory.new", lambda s: (_SlowCollector(s), [mem_scalar])
    )
    monkeypatch.setattr("swanlab.sdk.internal.probe_python.monitor.Apple.new", lambda _shim: None)
    monkeypatch.setattr("swanlab.sdk.internal.pkg.timer.Timer.start", lambda self: self)
    monkeypatch.setattr("swanlab.sdk.internal.probe_python.monitor.COLLECT_TIMEOUT", 0.12)
    monkeypatch.setattr(console, "trace", MagicMock())

    monitor = Monitor(shim, core).start(ctx)
    assert monitor is not None
    started = time.monotonic()
    try:
        monitor._timer._execute_once()  # type: ignore[union-attr]
        elapsed = time.monotonic() - started
    finally:
        release.set()
        monitor.stop()

    assert elapsed < 0.21
