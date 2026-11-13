"""CollectorProtocol.collect 的失败降噪测试。

契约：
- 指标采集失败不再每条 ERROR 刷屏：同一指标只告警一次；
- 同一轮新失败的指标聚合成一条 warning；
- 后续失败只写 debug，指标恢复后继续正常返回。
"""

from typing import List

from swanlab.sdk.internal.pkg import console as pkg_console
from swanlab.sdk.internal.probe_python.protocol import CollectorProtocol


def _collector(*handlers) -> CollectorProtocol:
    collector = CollectorProtocol(None)  # type: ignore[arg-type]
    collector._handlers.extend(handlers)
    return collector


def test_collect_returns_handler_values():
    collector = _collector(("a", lambda: 1), ("b", lambda: 2))

    assert collector.collect() == [("a", 1), ("b", 2)]


def test_first_failure_warns_once_and_aggregates(monkeypatch):
    warnings: List[tuple] = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))

    def unsupported():
        raise RuntimeError("NVMLError_NotSupported: Not Supported")

    collector = _collector(
        ("gpu.0.power", unsupported),
        ("gpu.1.power", unsupported),
        ("gpu.0.pct", lambda: 50),
    )

    assert collector.collect() == [("gpu.0.pct", 50)]
    assert len(warnings) == 1
    message = warnings[0][0]
    assert "gpu.0.power" in message
    assert "gpu.1.power" in message
    assert "NVMLError_NotSupported" in message

    # 后续轮次不再重复告警
    assert collector.collect() == [("gpu.0.pct", 50)]
    assert len(warnings) == 1


def test_metric_recovers_without_extra_warning(monkeypatch):
    warnings: List[tuple] = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    state = {"fail": True}

    def flaky():
        if state["fail"]:
            raise RuntimeError("boom")
        return 42

    collector = _collector(("gpu.0.power", flaky))

    assert collector.collect() == []
    assert len(warnings) == 1

    state["fail"] = False
    assert collector.collect() == [("gpu.0.power", 42)]
    assert len(warnings) == 1


def test_repeated_failure_is_debug_only(monkeypatch):
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: None)
    debugs: List[tuple] = []
    monkeypatch.setattr(pkg_console, "debug", lambda *args, **kwargs: debugs.append((args, kwargs)))

    def unsupported():
        raise RuntimeError("Not Supported")

    collector = _collector(("gpu.0.power", unsupported))

    collector.collect()
    collector.collect()

    assert len(debugs) == 1
    assert debugs[0][1]["write_to_tty"] is False
