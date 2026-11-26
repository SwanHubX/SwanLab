"""
@author: mx
@file: test_consumer_score.py
@time: 2026/11/25
@description: 后台消费者的评测分数事件分流
"""

import math
import queue
from pathlib import Path
from types import SimpleNamespace
from typing import cast
from unittest.mock import MagicMock

import numpy as np
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.sdk.internal.bus.events import ScoreLogEvent
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.run.components.consumer import BackgroundConsumer


def _make_consumer(tmp_path: Path) -> BackgroundConsumer:
    ctx = SimpleNamespace(
        core=MagicMock(),
        callbacker=MagicMock(),
        media_dir=tmp_path,
        config=SimpleNamespace(
            settings=SimpleNamespace(core=SimpleNamespace(section_rule=0, media_max_size=50 * 1024**2))
        ),
    )
    run_ctx = cast(RunContext, cast(object, ctx))
    return BackgroundConsumer(run_ctx, queue.Queue())


def _score_event(data: dict, dataset: str = "MMLU", seconds: int = 1) -> ScoreLogEvent:
    return ScoreLogEvent(dataset=dataset, data=data, timestamp=Timestamp(seconds=seconds))


class TestScoreRecords:
    def test_event_builds_score_records(self, tmp_path: Path):
        consumer = _make_consumer(tmp_path)
        consumer._handle_event(
            _score_event(
                {
                    "accuracy": 0.72,
                    "as_true": True,
                    "as_nan": float("nan"),
                    "as_inf": float("inf"),
                    "as_np": np.float32(0.5),
                }
            )
        )

        records = {record.key: record for record in consumer._score_batch}
        assert records["accuracy"].dataset == "MMLU"
        assert records["accuracy"].value == 0.72
        assert records["as_true"].value == 1.0
        assert math.isnan(records["as_nan"].value)
        assert math.isnan(records["as_inf"].value)
        assert records["as_np"].value == 0.5
        assert records["accuracy"].timestamp.ToMilliseconds() == 1000

    def test_same_key_timestamp_increments(self, tmp_path: Path):
        consumer = _make_consumer(tmp_path)
        consumer._handle_event(_score_event({"accuracy": 0.72}))
        consumer._handle_event(_score_event({"accuracy": 0.85}))

        first, second = consumer._score_batch
        assert second.timestamp.ToMilliseconds() == first.timestamp.ToMilliseconds() + 1

    def test_different_datasets_share_key_name(self, tmp_path: Path):
        consumer = _make_consumer(tmp_path)
        consumer._handle_event(_score_event({"accuracy": 0.72}, dataset="MMLU"))
        consumer._handle_event(_score_event({"accuracy": 0.85}, dataset="GSM8K"))

        first, second = consumer._score_batch
        assert first.timestamp.ToMilliseconds() == second.timestamp.ToMilliseconds()

    def test_invalid_key_and_value_ignored(self, tmp_path: Path):
        consumer = _make_consumer(tmp_path)
        consumer._handle_event(_score_event({"": 1.0, "k" * 513: 1.0, "bad": object(), "ok": 1.0}))

        assert [record.key for record in consumer._score_batch] == ["ok"]


class TestScoreFlush:
    def test_flush_calls_upsert_scores(self, tmp_path: Path):
        consumer = _make_consumer(tmp_path)
        core = cast(MagicMock, consumer._core)
        consumer._handle_event(_score_event({"accuracy": 0.72}))

        consumer._flush()

        core.upsert_scores.assert_called_once()
        assert consumer._score_batch == []

    def test_flush_failure_restores_score_batch(self, tmp_path: Path):
        consumer = _make_consumer(tmp_path)
        core = cast(MagicMock, consumer._core)
        core.upsert_scores.side_effect = RuntimeError("core failed")
        consumer._handle_event(_score_event({"accuracy": 0.72}))

        consumer._flush()

        assert len(consumer._score_batch) == 1
        assert consumer._score_batch[0].key == "accuracy"
