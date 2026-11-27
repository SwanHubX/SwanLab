from pathlib import Path

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogLevel
from swanlab.sdk.internal.bus.events import FileSaveEvent, LogEvent
from swanlab.sdk.protocol import CoreProtocolError


def _make_save_event(name: str, tmp_path: Path, policy: str) -> FileSaveEvent:
    return FileSaveEvent(
        source_path=str(tmp_path / name),
        name=name,
        policy=policy,
    )


# ── live ──


def test_live_save_appends_to_batch(tmp_path: Path, make_save_consumer):
    consumer, _, builder = make_save_consumer()
    event = _make_save_event("model.pt", tmp_path, "live")

    consumer._handle_event(event)

    builder.build_save.assert_called_once_with(event)
    assert list(consumer._save_batch) == [builder.build_save.return_value]


# ── now ──


def test_now_save_appends_to_batch(tmp_path: Path, make_save_consumer):
    consumer, _, builder = make_save_consumer()
    event = _make_save_event("model.pt", tmp_path, "now")

    consumer._handle_event(event)

    builder.build_save.assert_called_once_with(event)
    assert list(consumer._save_batch) == [builder.build_save.return_value]


# ── end ──


def test_end_save_appends_to_batch(tmp_path: Path, make_save_consumer):
    consumer, _, builder = make_save_consumer()
    event = _make_save_event("model.pt", tmp_path, "end")

    consumer._handle_event(event)

    builder.build_save.assert_called_once_with(event)
    assert list(consumer._save_batch) == [builder.build_save.return_value]


# ── multiple events ──


def test_mixed_policies_batch_together(tmp_path: Path, make_save_consumer):
    consumer, _, builder = make_save_consumer()
    events = [
        _make_save_event("a.pt", tmp_path, "now"),
        _make_save_event("b.pt", tmp_path, "live"),
    ]

    for e in events:
        consumer._handle_event(e)

    assert builder.build_save.call_count == 2
    assert len(consumer._save_batch) == 2


# ── flush ──


def test_flush_pushes_save_batch_to_core(make_save_consumer):
    consumer, core, builder = make_save_consumer()
    record = builder.build_save.return_value
    consumer._save_batch.append(record)

    consumer._flush()

    core.upsert_saves.assert_called_once_with([record])
    assert consumer._save_batch == []


def test_flush_does_not_call_upsert_saves_when_save_batch_empty(make_save_consumer):
    consumer, core, _ = make_save_consumer()

    consumer._flush()

    core.upsert_saves.assert_not_called()


def test_save_batch_flushes_alongside_other_batch_types(make_save_consumer):
    consumer, core, builder = make_save_consumer()
    save_record = builder.build_save.return_value
    consumer._save_batch.append(save_record)
    consumer._log_batch.append(
        builder.build_log(LogEvent(line="hello", level=LogLevel.LOG_LEVEL_INFO, timestamp=Timestamp()))
    )

    consumer._flush()

    core.upsert_saves.assert_called_once_with([save_record])
    core.upsert_logs.assert_called_once()
    assert consumer._save_batch == []
    assert consumer._log_batch == []


def test_save_batch_auto_flushes_when_full(tmp_path: Path, make_save_consumer):
    consumer, core, _ = make_save_consumer(batch_size=3)
    for name in ("a.pt", "b.pt", "c.pt"):
        consumer._handle_event(_make_save_event(name, tmp_path, "now"))

    assert consumer._batch_full
    consumer._flush()
    core.upsert_saves.assert_called_once()
    assert consumer._save_batch == []


# ── flush failure recovery ──


def test_failed_upsert_saves_restores_batch(make_save_consumer):
    consumer, core, builder = make_save_consumer()
    consumer._save_batch.append(builder.build_save.return_value)
    core.upsert_saves.side_effect = CoreProtocolError("save failed", retryable=True)

    consumer._flush()

    assert len(consumer._save_batch) == 1


def test_failed_upsert_saves_does_not_double_restore_successful_batches(make_save_consumer):
    """upsert_saves 失败时，已成功的 log 批次不会被恢复。"""
    consumer, core, builder = make_save_consumer()
    consumer._save_batch.append(builder.build_save.return_value)
    consumer._log_batch.append(
        builder.build_log(LogEvent(line="hello", level=LogLevel.LOG_LEVEL_INFO, timestamp=Timestamp()))
    )
    core.upsert_saves.side_effect = CoreProtocolError("save failed", retryable=True)

    consumer._flush()

    assert consumer._log_batch == []
    assert len(consumer._save_batch) == 1
