import queue
from pathlib import Path
from types import SimpleNamespace
from typing import cast
from unittest.mock import MagicMock

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import UpdateType
from swanlab.proto.swanlab.save.v1.save_pb2 import SavePolicy, SaveRecord
from swanlab.proto.swanlab.system.v1.console_pb2 import StreamType
from swanlab.sdk.internal.bus.events import ConfigEvent, ConsoleEvent, FileSaveEvent
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.run.components.consumer import BackgroundConsumer


def _make_consumer(tmp_path: Path, batch_size: int = 100):
    core = MagicMock()
    ctx = SimpleNamespace(
        core=core,
        run_dir=tmp_path,
        metrics=MagicMock(),
        callbacker=MagicMock(),
    )
    builder = MagicMock()
    builder.build_save.return_value = SaveRecord(
        name="model.pt",
        source_path=str(tmp_path / "model.pt"),
        policy=SavePolicy.SAVE_POLICY_NOW,
    )
    return BackgroundConsumer(cast(RunContext, ctx), queue.Queue(), builder, batch_size=batch_size), core, builder


def _make_save_event(name: str, tmp_path: Path, policy: str) -> FileSaveEvent:
    return FileSaveEvent(
        source_path=str(tmp_path / name),
        name=name,
        policy=policy,
    )


# ── live ──


def test_live_save_appends_to_batch(tmp_path: Path):
    consumer, _, builder = _make_consumer(tmp_path)
    event = _make_save_event("model.pt", tmp_path, "live")

    consumer._handle_event(event)

    builder.build_save.assert_called_once_with(event)
    assert list(consumer._save_batch) == [builder.build_save.return_value]


# ── now ──


def test_now_save_appends_to_batch(tmp_path: Path):
    consumer, _, builder = _make_consumer(tmp_path)
    event = _make_save_event("model.pt", tmp_path, "now")

    consumer._handle_event(event)

    builder.build_save.assert_called_once_with(event)
    assert list(consumer._save_batch) == [builder.build_save.return_value]


# ── end ──


def test_end_save_appends_to_batch(tmp_path: Path):
    consumer, _, builder = _make_consumer(tmp_path)
    event = _make_save_event("model.pt", tmp_path, "end")

    consumer._handle_event(event)

    builder.build_save.assert_called_once_with(event)
    assert list(consumer._save_batch) == [builder.build_save.return_value]


# ── multiple events ──


def test_mixed_policies_batch_together(tmp_path: Path):
    consumer, _, builder = _make_consumer(tmp_path)
    events = [
        _make_save_event("a.pt", tmp_path, "now"),
        _make_save_event("b.pt", tmp_path, "live"),
    ]

    for e in events:
        consumer._handle_event(e)

    assert builder.build_save.call_count == 2
    assert len(consumer._save_batch) == 2


# ── flush ──


def test_flush_pushes_save_batch_to_core(tmp_path: Path):
    consumer, core, builder = _make_consumer(tmp_path)
    record = builder.build_save.return_value
    consumer._save_batch.append(record)

    consumer._flush()

    core.upsert_saves.assert_called_once_with([record])
    assert consumer._save_batch == []


def test_flush_does_not_call_upsert_saves_when_save_batch_empty(tmp_path: Path):
    consumer, core, _ = _make_consumer(tmp_path)

    consumer._flush()

    core.upsert_saves.assert_not_called()


def test_save_batch_flushes_alongside_other_batch_types(tmp_path: Path):
    consumer, core, builder = _make_consumer(tmp_path)
    save_record = builder.build_save.return_value
    consumer._save_batch.append(save_record)
    consumer._console_batch.append(
        builder.build_console(ConsoleEvent(line="hello", stream=StreamType.STREAM_TYPE_STDOUT, timestamp=Timestamp()))
    )

    consumer._flush()

    core.upsert_saves.assert_called_once_with([save_record])
    core.upsert_consoles.assert_called_once()
    assert consumer._save_batch == []
    assert consumer._console_batch == []


def test_save_batch_auto_flushes_when_full(tmp_path: Path):
    consumer, core, _ = _make_consumer(tmp_path, batch_size=3)
    for name in ("a.pt", "b.pt", "c.pt"):
        consumer._handle_event(_make_save_event(name, tmp_path, "now"))

    assert consumer._batch_full
    consumer._flush()
    core.upsert_saves.assert_called_once()
    assert consumer._save_batch == []


# ── flush failure recovery ──


def test_failed_upsert_saves_restores_batch(tmp_path: Path):
    consumer, core, builder = _make_consumer(tmp_path)
    consumer._save_batch.append(builder.build_save.return_value)
    core.upsert_saves.side_effect = RuntimeError("save failed")

    consumer._flush()

    assert len(consumer._save_batch) == 1


def test_failed_upsert_saves_does_not_double_restore_successful_batches(tmp_path: Path):
    """upsert_saves 失败时，已成功的 config 批次不会被恢复。"""
    consumer, core, builder = _make_consumer(tmp_path)
    consumer._save_batch.append(builder.build_save.return_value)
    consumer._config_batch.append(
        builder.build_config(ConfigEvent(update=UpdateType.UPDATE_TYPE_PATCH, timestamp=Timestamp()))
    )
    core.upsert_saves.side_effect = RuntimeError("save failed")

    consumer._flush()

    assert consumer._config_batch == []
    assert len(consumer._save_batch) == 1
