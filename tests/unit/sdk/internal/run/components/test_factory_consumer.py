"""consumer 出批上限测试：单次 RPC 不能超过 core 的 record_batch。"""

import queue
from types import SimpleNamespace
from typing import cast
from unittest.mock import MagicMock

from swanlab.sdk.internal.bus import EmitterProtocol
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.run.components import _factory_consumer
from swanlab.sdk.internal.run.components.consumer import BackgroundConsumer


def _make_consumer(record_batch: int) -> BackgroundConsumer:
    ctx = SimpleNamespace(
        core=MagicMock(),
        callbacker=MagicMock(),
        config=SimpleNamespace(settings=SimpleNamespace(mode="local", core=SimpleNamespace(record_batch=record_batch))),
    )
    emitter = SimpleNamespace(queue=queue.Queue())
    consumer = _factory_consumer(cast(RunContext, ctx), cast(EmitterProtocol, emitter))
    assert isinstance(consumer, BackgroundConsumer)
    return consumer


def test_batch_size_capped_by_record_batch():
    c7 = _make_consumer(7)
    assert c7._batch_size == 7
    assert c7._rpc_batch_size == 7

    c100 = _make_consumer(100)
    assert c100._batch_size == 100
    assert c100._rpc_batch_size == 100

    c10k = _make_consumer(10_000)
    assert c10k._batch_size == 100
    assert c10k._rpc_batch_size == 2_000
    assert c10k._stop_timeout is None
