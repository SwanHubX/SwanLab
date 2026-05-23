import threading

import pytest

from swanlab.proto.swanlab.operation.v1 import operation_pb2
from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState
from swanlab.sdk.internal.core_python.transport import tracker as tracker_module
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker

_MIB = 1024 * 1024


def test_tracker_defaults_to_not_started_empty_snapshot():
    tracker = UploadTracker()

    stats = tracker.snapshot()

    assert stats.state == CoreState.CORE_STATE_NOT_STARTED
    assert stats.total_number == 0
    assert stats.uploaded_number == 0
    assert stats.total_size == 0
    assert stats.uploaded_size == 0
    assert stats.total_records == 0
    assert stats.uploaded_records == 0
    assert stats.rate == 0.0
    assert list(stats.files) == []


def test_tracker_tracks_aggregate_progress_and_state():
    tracker = UploadTracker()

    tracker.set_state(CoreState.CORE_STATE_RUNNING)
    tracker.add_total(100)
    tracker.advance(30)

    stats = tracker.snapshot()
    assert stats.state == CoreState.CORE_STATE_RUNNING
    assert stats.total_number == 100
    assert stats.uploaded_number == 30


def test_tracker_file_progress_is_optional_and_limited():
    tracker = UploadTracker(uploading_files_limit=2)

    for index in range(4):
        tracker.add_file(f"file-{index}", f"/tmp/file-{index}.bin", 100)
    tracker.update_file_progress("file-1", "part-0", 40)

    stats = tracker.snapshot()
    files = list(stats.files)
    assert stats.total_number == 4
    assert stats.total_size == 400
    assert stats.uploaded_size == 40
    assert [file.path for file in files] == ["/tmp/file-1.bin", "/tmp/file-0.bin"]
    assert files[0].total == 100
    assert files[0].uploaded == 40


def test_tracker_snapshot_uses_proto_file_progress_message():
    tracker = UploadTracker()

    tracker.add_file("file", "/tmp/file.bin", 100)

    proto_file_progress_type = getattr(operation_pb2, "FileProgress", None)
    assert proto_file_progress_type is not None
    stats = tracker.snapshot()
    files = list(stats.files)
    assert len(files) == 1
    assert isinstance(files[0], proto_file_progress_type)


def test_tracker_finish_file_removes_file_display_but_keeps_aggregate_progress():
    tracker = UploadTracker()

    tracker.add_file("ckpt", "/tmp/checkpoint.pt", 200)
    tracker.update_file_progress("ckpt", "part-0", 200)
    tracker.finish_file("ckpt")

    stats = tracker.snapshot()
    assert stats.total_number == 1
    assert stats.uploaded_number == 1
    assert stats.total_size == 200
    assert stats.uploaded_size == 200
    assert list(stats.files) == []


def test_tracker_repeated_file_key_is_not_registered_twice():
    tracker = UploadTracker()

    tracker.add_file("same", "/tmp/first.bin", 100)
    tracker.add_file("same", "/tmp/second.bin", 100)

    stats = tracker.snapshot()
    assert stats.total_size == 100
    assert [file.path for file in stats.files] == ["/tmp/first.bin"]


def test_tracker_aggregate_progress_is_thread_safe():
    tracker = UploadTracker()
    tracker.add_total(1000)

    def advance_many():
        for _ in range(100):
            tracker.advance(1)

    threads = [threading.Thread(target=advance_many) for _ in range(10)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    stats = tracker.snapshot()
    assert stats.uploaded_number == 1000


def test_tracker_allows_record_progress_before_total_is_known():
    tracker = UploadTracker()

    tracker.advance_records(3)

    stats = tracker.snapshot()
    assert stats.total_records == 0
    assert stats.uploaded_records == 3

    tracker.set_total_records(10)

    stats = tracker.snapshot()
    assert stats.total_records == 10
    assert stats.uploaded_records == 3


def test_tracker_file_rate_waits_for_sample_window_before_display(monkeypatch):
    clock = {"now": 0.0}
    monkeypatch.setattr(tracker_module.time, "monotonic", lambda: clock["now"])
    tracker = UploadTracker()

    tracker.add_file("model", "/tmp/model.bin", 10 * _MIB)

    clock["now"] = 0.01
    tracker.update_file_progress("model", "model", _MIB)
    stats = tracker.snapshot()
    assert stats.rate == 0.0
    assert stats.files[0].rate == 0.0

    clock["now"] = 0.02
    tracker.update_file_progress("model", "model", 2 * _MIB)
    stats = tracker.snapshot()
    assert stats.rate == 0.0
    assert stats.files[0].rate == 0.0

    clock["now"] = 1.0
    tracker.update_file_progress("model", "model", 3 * _MIB)
    stats = tracker.snapshot()
    assert stats.rate == pytest.approx(3 * _MIB)
    assert stats.files[0].rate == pytest.approx(3 * _MIB)


def test_tracker_file_rate_applies_ema_between_sample_windows(monkeypatch):
    clock = {"now": 0.0}
    monkeypatch.setattr(tracker_module.time, "monotonic", lambda: clock["now"])
    tracker = UploadTracker()

    tracker.add_file("model", "/tmp/model.bin", 10 * _MIB)

    clock["now"] = 1.0
    tracker.update_file_progress("model", "model", 4 * _MIB)
    stats = tracker.snapshot()
    assert stats.rate == pytest.approx(4 * _MIB)
    assert stats.files[0].rate == pytest.approx(4 * _MIB)

    clock["now"] = 2.0
    tracker.update_file_progress("model", "model", 5 * _MIB)
    stats = tracker.snapshot()
    expected_rate = (0.3 * 1 * _MIB) + (0.7 * 4 * _MIB)
    assert stats.rate == pytest.approx(expected_rate)
    assert stats.files[0].rate == pytest.approx(expected_rate)


def test_tracker_rate_decays_over_time_when_inactive(monkeypatch):
    clock = {"now": 0.0}
    monkeypatch.setattr(tracker_module.time, "monotonic", lambda: clock["now"])
    tracker = UploadTracker()

    tracker.add_file("model", "/tmp/model.bin", 10 * _MIB)

    clock["now"] = 1.0
    tracker.update_file_progress("model", "model", 4 * _MIB)
    stats = tracker.snapshot()
    assert stats.rate == pytest.approx(4 * _MIB)
    assert stats.files[0].rate == pytest.approx(4 * _MIB)

    # Move forward by one sample interval without update.
    clock["now"] = 1.0 + UploadTracker._RATE_SAMPLE_INTERVAL
    stats = tracker.snapshot()
    expected_decayed = 4 * _MIB * 0.7
    assert stats.rate == pytest.approx(expected_decayed)
    assert stats.files[0].rate == pytest.approx(expected_decayed)

    # Move forward by 10 sample intervals since last update.
    clock["now"] = 1.0 + (10 * UploadTracker._RATE_SAMPLE_INTERVAL)
    stats = tracker.snapshot()
    expected_decayed_more = 4 * _MIB * (0.7**10)
    assert stats.rate == pytest.approx(expected_decayed_more)
    assert stats.files[0].rate == pytest.approx(expected_decayed_more)
