"""BackgroundConsumer 测试共享夹具与 FakeClock。"""

import queue
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Dict, Optional, Tuple, cast
from unittest.mock import MagicMock

import pytest

from swanlab.proto.swanlab.save.v1.save_pb2 import SavePolicy, SaveRecord
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.run.components.consumer import BackgroundConsumer


class FakeClock:
    """可控的单调时钟测试替身。"""

    def __init__(self) -> None:
        self.now = 100.0

    def __call__(self) -> float:
        return self.now

    def advance(self, seconds: float) -> None:
        self.now += seconds


def _build_consumer_context(work_dir: Path, core: MagicMock) -> RunContext:
    ctx = SimpleNamespace(
        core=core,
        callbacker=MagicMock(),
        metrics=MagicMock(),
        media_dir=work_dir,
        run_dir=work_dir,
        config=SimpleNamespace(
            settings=SimpleNamespace(
                core=SimpleNamespace(section_rule=0, skip_store=False, media_max_size=50 * 1024**2)
            )
        ),
    )
    return cast(RunContext, cast(object, ctx))


@pytest.fixture
def make_consumer(tmp_path: Path):
    """创建基础 BackgroundConsumer 实例的工厂夹具。"""

    def _create(work_dir: Optional[Path] = None) -> BackgroundConsumer:
        return BackgroundConsumer(_build_consumer_context(work_dir or tmp_path, MagicMock()), queue.Queue())

    return _create


@pytest.fixture
def make_retry_consumer(tmp_path: Path):
    """创建带 FakeClock 的 BackgroundConsumer 及其 mock core 的工厂夹具。"""

    def _create(
        work_dir: Optional[Path] = None,
        *,
        event_queue: Optional[queue.Queue] = None,
        retry_base_delay: float = 0.5,
        retry_max_delay: float = 2.0,
        flush_timeout: float = 0.5,
        batch_size: int = 100,
        rpc_batch_size: Optional[int] = None,
        stop_timeout: Optional[float] = 0.0,
        sleep: Optional[Callable[[float], None]] = None,
        max_buffered_records: int = 50_000,
        max_buffered_bytes: int = 128 << 20,
    ) -> Tuple[BackgroundConsumer, MagicMock, FakeClock]:
        core = MagicMock()
        clock = FakeClock()
        kwargs: Dict[str, Any] = {
            "flush_timeout": flush_timeout,
            "batch_size": batch_size,
            "stop_timeout": stop_timeout,
            "retry_base_delay": retry_base_delay,
            "retry_max_delay": retry_max_delay,
            "clock": clock,
            "max_buffered_records": max_buffered_records,
            "max_buffered_bytes": max_buffered_bytes,
        }
        if rpc_batch_size is not None:
            kwargs["rpc_batch_size"] = rpc_batch_size
        if sleep is not None:
            kwargs["sleep"] = sleep
        consumer = BackgroundConsumer(
            _build_consumer_context(work_dir or tmp_path, core),
            queue.Queue() if event_queue is None else event_queue,
            **kwargs,
        )
        return consumer, core, clock

    return _create


@pytest.fixture
def make_save_consumer(tmp_path: Path):
    """创建带 mock builder 的 BackgroundConsumer 工厂夹具（用于 save 批次测试）。"""

    def _create(batch_size: int = 100) -> Tuple[BackgroundConsumer, MagicMock, MagicMock]:
        core = MagicMock()
        builder = MagicMock()
        builder.build_save.return_value = SaveRecord(
            name="model.pt",
            source_path=str(tmp_path / "model.pt"),
            policy=SavePolicy.SAVE_POLICY_NOW,
        )
        consumer = BackgroundConsumer(_build_consumer_context(tmp_path, core), queue.Queue(), batch_size=batch_size)
        # builder 由 consumer 内部创建，测试以 mock 替换以拦截 build_save 结果
        consumer._builder = builder
        return consumer, core, builder

    return _create
