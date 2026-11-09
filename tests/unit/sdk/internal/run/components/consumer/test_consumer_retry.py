"""BackgroundConsumer 上行失败分类：可重试回塞、终态丢弃、有界退避。

关键契约：
- 可重试错误必须回塞原批次，退避到期后重新提交，记录不丢；
- 终态错误只告警一次并终结该批，不回塞；
- 退避通过 queue.get(timeout=...) 表达，不 sleep 轮询。
"""

import queue
from pathlib import Path
from types import SimpleNamespace
from typing import cast
from unittest.mock import MagicMock

import pytest

from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScalarRecord
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console as pkg_console
from swanlab.sdk.internal.run.components.consumer import BackgroundConsumer
from swanlab.sdk.protocol import CoreProtocolError


class FakeClock:
    def __init__(self) -> None:
        self.now = 100.0

    def __call__(self) -> float:
        return self.now

    def advance(self, seconds: float) -> None:
        self.now += seconds


def _make_consumer(
    tmp_path: Path = Path("."),
    *,
    retry_base_delay: float = 0.5,
    retry_max_delay: float = 2.0,
    flush_timeout: float = 0.5,
):
    core = MagicMock()
    ctx = SimpleNamespace(core=core, callbacker=MagicMock(), media_dir=tmp_path)
    clock = FakeClock()
    consumer = BackgroundConsumer(
        cast(RunContext, cast(object, ctx)),
        queue.Queue(),
        flush_timeout=flush_timeout,
        retry_base_delay=retry_base_delay,
        retry_max_delay=retry_max_delay,
        clock=clock,
    )
    return consumer, core, clock


def _scalar(key: str = "loss") -> ScalarRecord:
    return ScalarRecord(key=key, step=1)


def test_retryable_error_restores_batch_and_schedules_backoff():
    consumer, core, clock = _make_consumer()
    record = _scalar()
    consumer._scalar_batch.append(record)
    core.upsert_scalars.side_effect = [CoreProtocolError("busy", retryable=True), None]

    consumer._flush()

    assert consumer._scalar_batch == [record]
    assert consumer._retry_at == clock.now + 0.5
    assert consumer._retry_due() is False

    clock.advance(0.5)
    assert consumer._retry_due() is True
    consumer._flush()

    assert consumer._scalar_batch == []
    assert core.upsert_scalars.call_count == 2
    assert core.upsert_scalars.call_args_list[0].args[0] == [record]
    assert core.upsert_scalars.call_args_list[1].args[0] == [record]
    assert consumer._retry_at is None
    assert consumer._retry_attempt == 0


def test_retryable_error_logs_no_warning(monkeypatch):
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    consumer, core, _ = _make_consumer()
    consumer._scalar_batch.append(_scalar())
    core.upsert_scalars.side_effect = CoreProtocolError("busy", retryable=True)

    consumer._flush()

    assert warnings == []


def test_terminal_error_drops_batch_with_one_warning(monkeypatch, tmp_path: Path):
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    consumer, core, clock = _make_consumer(tmp_path)
    consumer._scalar_batch.append(_scalar())
    core.upsert_scalars.side_effect = CoreProtocolError("bad record", retryable=False)

    consumer._flush()

    assert consumer._scalar_batch == []
    assert consumer._retry_at is None
    assert len(warnings) == 1

    consumer._flush()  # 空批次不再产生日志
    assert len(warnings) == 1


def test_unknown_exception_is_retried():
    consumer, core, _ = _make_consumer()
    record = _scalar()
    consumer._scalar_batch.append(record)
    core.upsert_scalars.side_effect = RuntimeError("socket closed")

    consumer._flush()

    assert consumer._scalar_batch == [record]
    assert consumer._retry_at is not None


def test_terminal_group_does_not_block_other_groups(monkeypatch):
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: None)
    consumer, core, _ = _make_consumer()
    consumer._scalar_batch.append(_scalar())
    consumer._log_batch.append(object())  # type: ignore[arg-type]
    core.upsert_scalars.side_effect = CoreProtocolError("bad record", retryable=False)

    consumer._flush()

    core.upsert_logs.assert_called_once()
    assert consumer._scalar_batch == []
    assert consumer._log_batch == []


def test_backoff_is_bounded():
    consumer, core, clock = _make_consumer(retry_base_delay=1.0, retry_max_delay=2.0)
    consumer._scalar_batch.append(_scalar())
    core.upsert_scalars.side_effect = CoreProtocolError("busy", retryable=True)

    delays = []
    for _ in range(6):
        consumer._flush()
        delays.append(consumer._retry_at - clock.now)
        clock.advance(consumer._retry_at - clock.now)

    assert delays[:3] == [1.0, 2.0, 2.0]
    assert max(delays) <= 2.0
    assert consumer._scalar_batch != []


@pytest.mark.parametrize(
    "retry_in, flush_timeout, expected",
    [
        (None, 0.5, 0.5),
        (10.0, 0.5, 0.5),
        (0.2, 0.5, 0.2),
        (0.0, 0.5, 0.0),
    ],
)
def test_queue_timeout_is_min_of_flush_timeout_and_remaining(retry_in, flush_timeout, expected):
    consumer, _, clock = _make_consumer(flush_timeout=flush_timeout)
    consumer._retry_at = None if retry_in is None else clock.now + retry_in

    assert consumer._queue_timeout() == pytest.approx(expected)
