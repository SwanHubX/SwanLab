from io import BytesIO
from pathlib import Path
from unittest.mock import ANY, MagicMock, patch

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem, MediaRecord, MediaValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
from swanlab.sdk.internal.core_python.utils import ProgressFileWrapper


class _FakeResponse:
    def __init__(self, etag: str = "etag") -> None:
        self.headers = {"ETag": etag}

    def raise_for_status(self) -> None:
        pass


class _FakeApiErrorResponse:
    def __init__(self, status_code: int) -> None:
        self.status_code = status_code
        self.request = MagicMock()
        self.url = "https://api.test/upload"


class _FakeSession:
    def __init__(self) -> None:
        self.puts: list[dict] = []

    def put(self, url, *, data, headers):
        payload = b""
        if hasattr(data, "read"):
            for chunk in iter(lambda: data.read(2), b""):
                payload += chunk
        self.puts.append({"url": url, "data": data, "payload": payload, "headers": headers})
        return _FakeResponse(etag=f"etag-{url}")


def _make_sender(
    tmp_path: Path,
    save_batch: int = 100,
    save_split: int = 100 * 1024 * 1024,
    save_part: int = 32 * 1024 * 1024,
) -> HttpRecordSender:
    ctx = CoreContext(
        config=CoreConfig(
            run_id="test-run-id",
            run_dir=tmp_path,
            section_rule=0,
            record_batch=10000,
            record_interval=5.0,
            save_split=save_split,
            save_size=50 * 1024 * 1024 * 1024,
            save_part=save_part,
            save_batch=save_batch,
        )
    )
    ctx.set_online_params(
        username="alice",
        project="demo",
        project_id="project-id",
        experiment_id="experiment-id",
    )
    return HttpRecordSender(ctx=ctx)


def _make_save_record(source: Path, name: str = "checkpoints/model.txt") -> Record:
    return Record(save=SaveRecord(name=name, source_path=str(source), target_path=str(source)))


def _make_media_record(filename: str) -> Record:
    timestamp = Timestamp()
    timestamp.GetCurrentTime()
    return Record(
        media=MediaRecord(
            key="examples/image",
            step=1,
            type=ColumnType.COLUMN_TYPE_IMAGE,
            timestamp=timestamp,
            value=MediaValue(items=[MediaItem(filename=filename)]),
        )
    )


def test_upload_save_delegates_small_file_s3_put_to_upload_resources(tmp_path: Path):
    source = tmp_path / "model.txt"
    source.write_text("weights", encoding="utf-8")
    sender = _make_sender(tmp_path)
    session = MagicMock()

    def _fake_upload(session, *, resources, tracker=None):
        return {r["source_path"] for r in resources}

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.prepare_save_files",
            return_value={"urls": ["https://s3.test/model"]},
        ) as mock_prepare,
        patch("swanlab.sdk.internal.core_python.transport.sender.complete_save_files") as mock_complete,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=session),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.upload_saves",
            side_effect=_fake_upload,
        ) as mock_upload_saves,
    ):
        sender.upload_save([_make_save_record(source)])

    mock_prepare.assert_called_once()
    mock_upload_saves.assert_called_once()
    assert mock_upload_saves.call_args.args == (session,)
    resource = mock_upload_saves.call_args.kwargs["resources"][0]
    assert resource == {
        "url": "https://s3.test/model",
        "source_path": str(source),
        "content_type": "text/plain",
        "size": source.stat().st_size,
        "tracker_key": ANY,
    }
    assert resource["tracker_key"].startswith(f"checkpoints/model.txt:{source.stat().st_size}:")
    mock_complete.assert_called_once_with(
        "experiment-id", files=[{"path": "checkpoints/model.txt", "state": "UPLOADED"}]
    )


def test_upload_save_skips_complete_when_s3_upload_fails(tmp_path: Path):
    source = tmp_path / "failed.bin"
    source.write_bytes(b"failed")
    sender = _make_sender(tmp_path)

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.prepare_save_files",
            return_value={"urls": ["https://s3.test/failed"]},
        ),
        patch("swanlab.sdk.internal.core_python.transport.sender.complete_save_files") as mock_complete,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch("swanlab.sdk.internal.core_python.transport.sender.upload_saves", side_effect=RuntimeError("s3 failed")),
        patch("swanlab.sdk.internal.pkg.safe.console.trace"),
    ):
        sender.upload_save([_make_save_record(source, name="failed.bin")])

    # 上传失败时仍应调用 complete，但标记为 FAILED
    mock_complete.assert_called_once_with("experiment-id", files=[{"path": "failed.bin", "state": "FAILED"}])


def test_upload_save_splits_pending_files_by_save_batch(tmp_path: Path):
    records = []
    for index in range(5):
        source = tmp_path / f"model-{index}.txt"
        source.write_text(f"weights-{index}", encoding="utf-8")
        records.append(_make_save_record(source, name=f"checkpoints/model-{index}.txt"))

    sender = _make_sender(tmp_path, save_batch=2)
    session = MagicMock()

    def _fake_prepare(_, *, files):
        return {"urls": [f"https://s3.test/{item['path']}" for item in files]}

    def _fake_upload(_, *, resources, tracker=None):
        return {r["source_path"] for r in resources}

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.prepare_save_files",
            side_effect=_fake_prepare,
        ) as mock_prepare,
        patch("swanlab.sdk.internal.core_python.transport.sender.complete_save_files") as mock_complete,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=session),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.upload_saves",
            side_effect=_fake_upload,
        ) as mock_upload_saves,
    ):
        sender.upload_save(records)

    assert mock_prepare.call_count == 3
    assert mock_upload_saves.call_count == 3
    assert mock_complete.call_count == 3
    assert [len(call.kwargs["files"]) for call in mock_prepare.call_args_list] == [2, 2, 1]
    assert [len(call.kwargs["files"]) for call in mock_complete.call_args_list] == [2, 2, 1]


def test_upload_save_tracks_small_file_progress_once_across_retries(tmp_path: Path):
    source = tmp_path / "model.txt"
    source.write_text("weights", encoding="utf-8")
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)

    def _fake_upload(_, *, resources, tracker=None):
        assert tracker is not None
        # 模拟 upload_saves 内部行为：finish_file
        for resource in resources:
            tracker.finish_file(resource["tracker_key"])
        return {resource["source_path"] for resource in resources}

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.prepare_save_files",
            return_value={"urls": ["https://s3.test/model"]},
        ),
        patch("swanlab.sdk.internal.core_python.transport.sender.complete_save_files"),
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch("swanlab.sdk.internal.core_python.transport.sender.upload_saves", side_effect=_fake_upload),
    ):
        sender.upload_save([_make_save_record(source)])
        sender.upload_save([_make_save_record(source)])

    stats = tracker.snapshot()
    assert stats.total_size == source.stat().st_size
    assert stats.uploaded_size == source.stat().st_size
    assert stats.total_number == 1
    assert stats.uploaded_number == 1
    assert list(stats.files) == []


def test_upload_save_tracks_multipart_parts_once_across_retries(tmp_path: Path):
    source = tmp_path / "large.bin"
    source.write_bytes(b"0123456789")
    tracker = UploadTracker()
    sender = _make_sender(tmp_path, save_split=5, save_part=4)
    sender.set_tracker(tracker)
    session = _FakeSession()

    multipart_response = {
        "files": [
            {
                "uploadId": "upload-id",
                "parts": [
                    {"partNumber": 2, "url": "https://s3.test/part-2"},
                    {"partNumber": 1, "url": "https://s3.test/part-1"},
                    {"partNumber": 3, "url": "https://s3.test/part-3"},
                ],
            }
        ]
    }

    tracked_keys = []
    original_track_file = sender._track_file

    def _capture_track_file(key: str, path: str, size: int) -> None:
        tracked_keys.append(key)
        original_track_file(key, path, size)

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.prepare_multipart_save",
            return_value=multipart_response,
        ),
        patch("swanlab.sdk.internal.core_python.transport.sender.complete_multipart_save") as mock_complete,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=session),
        patch.object(sender, "_track_file", side_effect=_capture_track_file),
    ):
        sender.upload_save([_make_save_record(source, name="large.bin")])
        sender.upload_save([_make_save_record(source, name="large.bin")])

    stats = tracker.snapshot()
    assert stats.total_size == source.stat().st_size
    assert stats.uploaded_size == source.stat().st_size
    assert stats.total_number == 1
    assert stats.uploaded_number == 1
    assert sorted(len(item["data"]) for item in session.puts) == [2, 2, 4, 4, 4, 4]
    assert sorted(item["payload"] for item in session.puts) == [b"0123", b"0123", b"4567", b"4567", b"89", b"89"]
    assert all(key.startswith(f"large.bin:{source.stat().st_size}:") for key in tracked_keys)
    assert all(isinstance(item["data"], ProgressFileWrapper) for item in session.puts)
    assert all(not isinstance(item["data"]._file_obj, BytesIO) for item in session.puts)
    assert mock_complete.call_count == 2


def test_upload_media_finishes_file_inside_progress_callback(tmp_path: Path):
    media_dir = tmp_path / "media" / "image"
    media_dir.mkdir(parents=True)
    first_path = media_dir / "first.png"
    second_path = media_dir / "second.png"
    first_path.write_bytes(b"one")
    second_path.write_bytes(b"two2")
    inner_tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(inner_tracker)
    snapshots = []

    def _fake_upload_resource(_, __, *, paths, buffers, tracker=None):
        assert paths == ["media/image/first.png", "media/image/second.png"]
        assert tracker is not None
        tracker.finish_file("media/image/first.png:3")
        snapshots.append(inner_tracker.snapshot())
        tracker.finish_file("media/image/second.png:4")
        snapshots.append(inner_tracker.snapshot())

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.upload_resource",
            side_effect=_fake_upload_resource,
        ),
        patch("swanlab.sdk.internal.core_python.transport.sender.upload_media"),
    ):
        sender.upload_media([_make_media_record("first.png"), _make_media_record("second.png")])

    assert snapshots[0].uploaded_size == 3
    assert snapshots[0].uploaded_number == 1
    assert [Path(file.path) for file in snapshots[0].files] == [second_path]
    assert snapshots[1].uploaded_size == 7
    assert snapshots[1].uploaded_number == 2
    assert list(snapshots[1].files) == []


def test_upload_advances_records_for_non_retryable_api_error(tmp_path: Path):
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)
    sender._upload_handlers["scalar"] = MagicMock(
        side_effect=ApiError(
            _FakeApiErrorResponse(400),
            method="POST",
            trace_id="trace-id",
            code="bad-request",
            message="bad request",
        )
    )

    sender.upload("scalar", [Record()])

    assert tracker.snapshot().uploaded_records == 1


def test_upload_does_not_advance_records_for_retryable_api_error(tmp_path: Path):
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)
    sender._upload_handlers["scalar"] = MagicMock(
        side_effect=ApiError(
            _FakeApiErrorResponse(500),
            method="POST",
            trace_id="trace-id",
            code="internal-error",
            message="server error",
        )
    )

    try:
        sender.upload("scalar", [Record()])
    except ApiError:
        pass

    assert tracker.snapshot().uploaded_records == 0
