"""
@description: 验证 _factory_consumer 把 CoreSettings.consumer_batch 透传给 BackgroundConsumer
"""

import queue
from pathlib import Path
from types import SimpleNamespace
from typing import cast
from unittest.mock import MagicMock

from swanlab.sdk.internal.bus.emitter import RunEmitter
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.run.components import _factory_consumer
from swanlab.sdk.internal.run.components.consumer import BackgroundConsumer
from swanlab.sdk.internal.run.components.null import NullConsumer
from swanlab.sdk.internal.settings import Settings


def _make_ctx(tmp_path: Path, *, mode: str = "online", consumer_batch: int = 100):
    """构造一个最小可用的 RunContext 替身。

    通过真实 Settings 拿到真实 CoreSettings，再用 model_copy 覆盖 consumer_batch，
    保证测的是真实字段链路 ctx.config.settings.core.consumer_batch。
    """
    settings = Settings()
    core = settings.core.model_copy(update={"consumer_batch": consumer_batch})
    settings_ns = SimpleNamespace(mode=mode, core=core)
    config_ns = SimpleNamespace(settings=settings_ns)
    ctx = SimpleNamespace(
        config=config_ns,
        core=MagicMock(),
        run_dir=tmp_path,
        metrics=MagicMock(),
        callbacker=MagicMock(),
    )
    return cast(RunContext, cast(object, ctx))


def _make_emitter() -> RunEmitter:
    e = MagicMock(spec=RunEmitter)
    e.queue = queue.Queue()
    return e


def test_factory_passes_consumer_batch_to_background_consumer(tmp_path: Path):
    ctx = _make_ctx(tmp_path, mode="online", consumer_batch=1234)
    consumer = _factory_consumer(ctx, _make_emitter(), MagicMock())

    assert isinstance(consumer, BackgroundConsumer)
    assert consumer._batch_size == 1234


def test_factory_default_consumer_batch_is_100(tmp_path: Path):
    ctx = _make_ctx(tmp_path, mode="online")  # 默认 100
    consumer = _factory_consumer(ctx, _make_emitter(), MagicMock())

    assert isinstance(consumer, BackgroundConsumer)
    assert consumer._batch_size == 100


def test_factory_returns_null_consumer_when_disabled(tmp_path: Path):
    ctx = _make_ctx(tmp_path, mode="disabled", consumer_batch=999)
    consumer = _factory_consumer(ctx, _make_emitter(), MagicMock())

    assert isinstance(consumer, NullConsumer)
