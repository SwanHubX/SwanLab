"""consumer 出批上限测试：单次 RPC 不能超过 core 的 record_batch。"""

import queue
from types import SimpleNamespace
from typing import cast
from unittest.mock import MagicMock

import pytest

from swanlab.sdk.internal.bus import EmitterProtocol
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.run.components import _factory_consumer
from swanlab.sdk.internal.run.components.consumer import BackgroundConsumer


def _make_consumer(record_batch: int) -> BackgroundConsumer:
    ctx = SimpleNamespace(
        core=MagicMock(),
        callbacker=MagicMock(),
        config=SimpleNamespace(
            settings=SimpleNamespace(
                mode="local",
                core=SimpleNamespace(record_batch=record_batch, skip_store=False, media_max_size=50 * 1024**2),
            )
        ),
    )
    emitter = SimpleNamespace(queue=queue.Queue())
    consumer = _factory_consumer(cast(RunContext, ctx), cast(EmitterProtocol, emitter))
    assert isinstance(consumer, BackgroundConsumer)
    return consumer


@pytest.mark.parametrize(
    ("record_batch", "batch_size", "rpc_batch_size"),
    [(7, 7, 7), (100, 100, 100), (10_000, 100, 2_000)],
)
def test_batch_size_capped_by_record_batch(record_batch, batch_size, rpc_batch_size):
    consumer = _make_consumer(record_batch)

    assert consumer._batch_size == batch_size
    assert consumer._rpc_batch_size == rpc_batch_size
    assert consumer._stop_timeout is None
