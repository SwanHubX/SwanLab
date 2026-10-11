from io import BytesIO
from pathlib import Path
from typing import Optional, cast
from unittest.mock import ANY, MagicMock, patch

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnRecord, ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem, MediaRecord, MediaValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord, SaveType
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender, encode_column
from swanlab.sdk.internal.core_python.transport.sender import files as sender_files
from swanlab.sdk.internal.core_python.transport.sender.utils import resolve_save_source
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
from swanlab.sdk.internal.core_python.utils import ProgressFileWrapper
from swanlab.sdk.typings.core_python.transport import UploadResult


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
    save_size: int = 50 * 1024 * 1024 * 1024,
    skip_store: bool = False,
) -> HttpRecordSender:
    ctx = CoreContext(
        config=CoreConfig(
            run_id="test-run-id",
            run_dir=tmp_path,
            section_rule=0,
            record_batch=10000,
            record_interval=5.0,
            save_split=save_split,
            save_size=save_size,
            save_part=save_part,
            save_batch=save_batch,
            skip_store=skip_store,
        )
    )
    ctx.set_online_params(
        username="alice",
        project="demo",
        project_id="project-id",
        project_version=1,
        experiment_id="experiment-id",
    )
    return HttpRecordSender(ctx=ctx)


def _make_save_record(
    source: Path,
    name: str = "checkpoints/model.txt",
    save_type: SaveType = SaveType.SAVE_TYPE_CUSTOM,
) -> Record:
    return Record(save=SaveRecord(name=name, source_path=str(source), target_path=str(source), type=save_type))


def _make_media_record(
    filename: str, media_type=ColumnType.COLUMN_TYPE_IMAGE, payload: Optional[bytes] = None
) -> Record:
    timestamp = Timestamp()
    timestamp.GetCurrentTime()
    # payload=None 表示未提供 payload（默认模式，从本地路径读取）；
    # payload=b"" 表示合法的空文件（skip_store），两者必须可区分。
    item = MediaItem(filename=filename)
    if payload is not None:
        item.payload = payload
    return Record(
        media=MediaRecord(
            key="examples/image",
            step=1,
            type=media_type,
            timestamp=timestamp,
            value=MediaValue(items=[item]),
        )
    )


def test_upload_column_batches_all_columns_without_dropping_tail(tmp_path: Path):
    sender = _make_sender(tmp_path)
    records = [
        Record(column=ColumnRecord(column_key=f"metric-{index}", column_type=ColumnType.COLUMN_TYPE_SCALAR))
        for index in range(3)
    ]

    with patch("swanlab.sdk.internal.core_python.transport.sender.upload_columns") as mock_upload_columns:
        sender.upload_column(records, batch_size=2)

    assert mock_upload_columns.call_count == 2
    assert [call.args[:2] for call in mock_upload_columns.call_args_list] == [("alice", "demo"), ("alice", "demo")]
    batches = [call.kwargs["columns"]["series"] for call in mock_upload_columns.call_args_list]
    assert [[column["key"] for column in batch] for batch in batches] == [["metric-0", "metric-1"], ["metric-2"]]


def test_upload_save_delegates_small_file_s3_put_to_upload_resources(tmp_path: Path):
    source = tmp_path / "model.txt"
    source.write_text("weights", encoding="utf-8")
    sender = _make_sender(tmp_path)
    session = MagicMock()

    def _fake_upload(session, *, resources, tracker=None):
        return {}

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.prepare_save_files",
            return_value={"urls": ["https://s3.test/model"]},
        ) as mock_prepare,
        patch("swanlab.sdk.internal.core_python.transport.sender.files.complete_save_files") as mock_complete,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=session),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.upload_saves",
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


def test_upload_save_retains_record_when_put_fails(tmp_path: Path):
    """PUT 失败的文件不向后端补报 FAILED（后端会删对象），记录保留在结果中等待重试。"""
    source = tmp_path / "failed.bin"
    source.write_bytes(b"failed")
    sender = _make_sender(tmp_path)

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.prepare_save_files",
            return_value={"urls": ["https://s3.test/failed"]},
        ),
        patch("swanlab.sdk.internal.core_python.transport.sender.files.complete_save_files") as mock_complete,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.upload_saves",
            return_value={0: RuntimeError("s3 failed")},
        ),
    ):
        record = _make_save_record(source, name="failed.bin")
        result = sender.upload_save([record])

    # 不调用 complete，记录保留并携带对象存储失败
    mock_complete.assert_not_called()
    assert result.confirmed_records == []
    assert result.remaining_records == [record]
    failure = result.failures[0]
    assert failure.source == "object_storage"
    assert failure.stage == "transfer"


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
        return {}

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.prepare_save_files",
            side_effect=_fake_prepare,
        ) as mock_prepare,
        patch("swanlab.sdk.internal.core_python.transport.sender.files.complete_save_files") as mock_complete,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=session),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.upload_saves",
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
        return {}

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.prepare_save_files",
            return_value={"urls": ["https://s3.test/model"]},
        ),
        patch("swanlab.sdk.internal.core_python.transport.sender.files.complete_save_files"),
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch("swanlab.sdk.internal.core_python.transport.sender.files.upload_saves", side_effect=_fake_upload),
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
    original_track_file = sender_files._track_file

    def _capture_track_file(tracker, key: str, path: str, size: int) -> None:
        tracked_keys.append(key)
        original_track_file(tracker, key, path, size)

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.prepare_multipart_save",
            return_value=multipart_response,
        ),
        patch("swanlab.sdk.internal.core_python.transport.sender.files.complete_multipart_save") as mock_complete,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=session),
        patch("swanlab.sdk.internal.core_python.transport.sender.files._track_file", side_effect=_capture_track_file),
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


def test_upload_media_passes_guessed_content_type(tmp_path: Path):
    media_dir = tmp_path / "media" / "html"
    media_dir.mkdir(parents=True)
    html_path = media_dir / "report.html"
    html_path.write_text("<h1>report</h1>", encoding="utf-8")
    sender = _make_sender(tmp_path)

    def _fake_upload_resource(_, __, *, paths, buffers, content_types=None, tracker=None):
        assert paths == ["media/html/report.html"]
        assert buffers == [html_path]
        assert content_types == ["text/html"]

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.upload_resource",
            side_effect=_fake_upload_resource,
        ) as mock_upload_resource,
        patch("swanlab.sdk.internal.core_python.transport.sender.upload_media") as mock_upload_media,
    ):
        sender.upload_media([_make_media_record("report.html", ColumnType.COLUMN_TYPE_HTML)])

    mock_upload_resource.assert_called_once()
    mock_upload_media.assert_called_once()


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

    def _fake_upload_resource(_, __, *, paths, buffers, content_types=None, tracker=None):
        assert paths == ["media/image/first.png", "media/image/second.png"]
        assert content_types == ["image/png", "image/png"]
        assert tracker is not None
        tracker.finish_file("media/image/first.png:3")
        snapshots.append(inner_tracker.snapshot())
        tracker.finish_file("media/image/second.png:4")
        snapshots.append(inner_tracker.snapshot())

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.upload_resource",
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


def test_upload_media_uses_payload_without_touching_disk(tmp_path: Path):
    """payload 走内存上传、不读本地文件；缺 payload 时告警丢弃。"""
    sender = _make_sender(tmp_path, skip_store=True)
    captured = {}

    def _fake_upload_resource(_, __, *, paths, buffers, content_types=None, tracker=None):
        captured["paths"] = paths
        captured["content_types"] = content_types
        captured["content"] = buffers[0].read()

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.upload_resource",
            side_effect=_fake_upload_resource,
        ),
        patch("swanlab.sdk.internal.core_python.transport.sender.upload_media") as mock_upload_media,
        patch("swanlab.sdk.internal.core_python.transport.sender.console.warning") as mock_warning,
    ):
        sender.upload_media([_make_media_record("ok.png", payload=b"ok"), _make_media_record("missing.png")])

    assert captured["paths"] == ["media/image/ok.png"]
    assert captured["content_types"] == ["image/png"]
    assert captured["content"] == b"ok"
    mock_upload_media.assert_called_once()
    assert "Media payload missing" in mock_warning.call_args.args[0]
    # 不得访问/创建本地 media 目录
    assert not (tmp_path / "media").exists()


def test_upload_media_uploads_empty_payload_as_zero_byte_file(tmp_path: Path):
    """skip_store 下合法的空文件（payload 存在但为 b""）必须上传，不能被当作缺失丢弃。"""
    sender = _make_sender(tmp_path, skip_store=True)
    captured = {}

    def _fake_upload_resource(_, __, *, paths, buffers, content_types=None, tracker=None):
        captured["paths"] = paths
        captured["content_types"] = content_types
        captured["content"] = buffers[0].read()

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.upload_resource",
            side_effect=_fake_upload_resource,
        ),
        patch("swanlab.sdk.internal.core_python.transport.sender.upload_media") as mock_upload_media,
        patch("swanlab.sdk.internal.core_python.transport.sender.console.warning") as mock_warning,
    ):
        sender.upload_media([_make_media_record("empty.txt", ColumnType.COLUMN_TYPE_TEXT, payload=b"")])

    assert captured["paths"] == ["media/text/empty.txt"]
    assert captured["content_types"] == ["text/plain"]
    assert captured["content"] == b""
    mock_upload_media.assert_called_once()
    mock_warning.assert_not_called()


def test_upload_does_not_advance_records(tmp_path: Path):
    """sender.upload 不推进统计：确认计数由 Dispatch 依据 confirmed_records 推进。"""
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)
    records = [Record(), Record()]
    sender._upload_handlers["scalar"] = MagicMock(return_value=UploadResult(confirmed_records=records))

    result = sender.upload("scalar", records)

    assert result.confirmed_records == records
    assert tracker.snapshot().uploaded_records == 0


def test_sender_resource_session_is_single_shot(tmp_path: Path):
    """对象存储 session 单次执行（default_retry=0）：媒体 PUT、save PUT 与分片子线程不自行重试。"""
    sender = _make_sender(tmp_path)

    with patch(
        "swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()
    ) as mock_create:
        session = sender._ensure_session()

    mock_create.assert_called_once_with(default_retry=0)
    assert session is sender._ensure_session()


def test_upload_maps_business_reject_to_failure_result(tmp_path: Path):
    """4xx 业务拒绝：转为结果内的失败，记录保留，不抛异常、不消费。"""
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)
    record = Record()
    sender._upload_handlers["scalar"] = MagicMock(
        side_effect=ApiError(
            _FakeApiErrorResponse(400),
            method="POST",
            trace_id="trace-id",
            code="bad-request",
            message="bad request",
        )
    )

    result = sender.upload("scalar", [record])

    assert result.confirmed_records == []
    assert result.remaining_records == [record]
    failure = result.failures[0]
    assert failure.status == 400
    assert failure.code == "bad-request"
    assert failure.source == "swanlab_api"
    assert tracker.snapshot().uploaded_records == 0


def test_upload_maps_server_error_to_failure_result(tmp_path: Path):
    """5xx 服务端错误：同样转为结果内的失败并保留记录，不向上抛出。"""
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)
    record = Record()
    sender._upload_handlers["scalar"] = MagicMock(
        side_effect=ApiError(
            _FakeApiErrorResponse(500),
            method="POST",
            trace_id="trace-id",
            code="internal-error",
            message="server error",
        )
    )

    result = sender.upload("scalar", [record])

    assert result.remaining_records == [record]
    assert result.failures[0].status == 500
    assert tracker.snapshot().uploaded_records == 0


def _make_internal_save_record(source_path: str, save_type: SaveType, name: str = "") -> Record:
    """构造内部保存（metadata/requirements/conda/config）记录，source_path 可指向不存在的训练机绝对路径。"""
    return Record(save=SaveRecord(name=name, source_path=source_path, target_path=source_path, type=save_type))


def test_resolve_save_source_prefers_primary_when_readable(tmp_path: Path):
    """本机运行时 source_path 可读 → 直接返回 primary，不走回退、不打 debug。"""
    source = tmp_path / "model.txt"
    source.write_text("weights", encoding="utf-8")
    rec = _make_save_record(source, name="checkpoints/model.txt")
    sender = _make_sender(tmp_path)

    with patch("swanlab.sdk.internal.core_python.transport.sender.console.debug") as mock_debug:
        result = resolve_save_source(sender._ctx, rec.save)

    assert result == source
    mock_debug.assert_not_called()


def test_resolve_save_source_falls_back_to_files_dir_for_custom(tmp_path: Path):
    """CUSTOM 跨挂载根：source_path 不可读 → 回退到 files_dir/save.name。"""
    # source_path 指向不存在路径，模拟训练机绝对路径在 sync 机不可读
    rec = _make_save_record(Path("/nonexistent/training-host/model.txt"), name="checkpoints/model.txt")
    sender = _make_sender(tmp_path)
    # files 子目录下有镜像（含子目录层级）
    fallback = tmp_path / "files" / "checkpoints" / "model.txt"
    fallback.parent.mkdir(parents=True)
    fallback.write_text("recovered", encoding="utf-8")

    with patch("swanlab.sdk.internal.core_python.transport.sender.console.debug") as mock_debug:
        result = resolve_save_source(sender._ctx, rec.save)

    assert result == fallback
    mock_debug.assert_called_once()


def test_resolve_save_source_falls_back_for_internal_metadata(tmp_path: Path):
    """内部 metadata：source_path 不可读 → 按 basename 回退到 files_dir/swanlab-metadata.json。"""
    rec = _make_internal_save_record("/nonexistent/host/swanlab-metadata.json", SaveType.SAVE_TYPE_METADATA)
    sender = _make_sender(tmp_path)
    fallback = tmp_path / "files" / "swanlab-metadata.json"
    fallback.parent.mkdir(parents=True)
    fallback.write_text("{}", encoding="utf-8")

    result = resolve_save_source(sender._ctx, rec.save)

    assert result == fallback


def test_resolve_save_source_uses_basename_not_name_for_config(tmp_path: Path):
    """config 的 save.name 为 "config"，但文件名是 config.yaml → 必须用 basename 而非 name 回退。"""
    rec = _make_internal_save_record("/nonexistent/host/config.yaml", SaveType.SAVE_TYPE_CONFIG, name="config")
    sender = _make_sender(tmp_path)
    fallback = tmp_path / "files" / "config.yaml"
    fallback.parent.mkdir(parents=True)
    fallback.write_text("key: value", encoding="utf-8")
    # files/config（即 save.name）不应被命中
    assert not (tmp_path / "files" / "config").exists()

    result = resolve_save_source(sender._ctx, rec.save)

    assert result == fallback


def test_resolve_save_source_returns_none_when_both_unreadable(tmp_path: Path):
    rec = _make_save_record(Path("/nonexistent/training-host/model.txt"), name="checkpoints/model.txt")
    sender = _make_sender(tmp_path)
    # 既无 primary 也无 files 镜像

    assert resolve_save_source(sender._ctx, rec.save) is None


def test_resolve_save_source_skips_stale_fallback_when_skip_store(tmp_path: Path):
    """skip_store：本 run 无镜像，源不可读时不回退——files 下的同名文件属于旧 run，不能上传。"""
    rec = _make_save_record(Path("/nonexistent/training-host/model.txt"), name="checkpoints/model.txt")
    sender = _make_sender(tmp_path, skip_store=True)
    fallback = tmp_path / "files" / "checkpoints" / "model.txt"
    fallback.parent.mkdir(parents=True)
    fallback.write_text("old-run-weights", encoding="utf-8")

    assert resolve_save_source(sender._ctx, rec.save) is None


def test_upload_save_recovers_internal_metadata_from_run_dir(tmp_path: Path):
    """内部 metadata 跨挂载根 sync：source_path 不可读，但 files/swanlab-metadata.json 有内容 → 成功上传。"""
    sender = _make_sender(tmp_path)
    # files 目录下放置真实 metadata（probe 直接写入）
    metadata_file = tmp_path / "files" / "swanlab-metadata.json"
    metadata_file.parent.mkdir(parents=True)
    metadata_file.write_text('{"hostname": "gpu-01"}', encoding="utf-8")
    # source_path 指向训练机绝对路径（sync 机不可读）
    rec = _make_internal_save_record("/nonexistent/host/swanlab-metadata.json", SaveType.SAVE_TYPE_METADATA)

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.upload_metadata") as mock_upload_meta,
        patch("swanlab.sdk.internal.core_python.transport.sender.console.debug"),
    ):
        sender.upload_save([rec])

    mock_upload_meta.assert_called_once()
    assert mock_upload_meta.call_args.kwargs["content"] == {"hostname": "gpu-01"}


def test_upload_save_recovers_custom_file_from_run_dir(tmp_path: Path):
    """用户保存跨挂载根 sync：source_path 不可读，files/checkpoints/model.txt 有镜像 → 走 save 上传。"""
    sender = _make_sender(tmp_path)
    files_dir = tmp_path / "files"
    files_dir.mkdir()
    recovered = files_dir / "checkpoints" / "model.txt"
    recovered.parent.mkdir(parents=True)
    recovered.write_text("recovered-weights", encoding="utf-8")
    # source_path 指向训练机绝对路径（含子目录），sync 机不可读
    rec = _make_save_record(Path("/nonexistent/host/checkpoints/model.txt"), name="checkpoints/model.txt")

    def _fake_upload(session, *, resources, tracker=None):
        # 断言读取自回退路径，而非训练机绝对路径
        assert all(r["source_path"] == str(recovered) for r in resources)
        return {}

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.prepare_save_files",
            return_value={"urls": ["https://s3.test/model"]},
        ) as mock_prepare,
        patch("swanlab.sdk.internal.core_python.transport.sender.files.complete_save_files") as mock_complete,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.upload_saves", side_effect=_fake_upload
        ) as mock_upload_saves,
        patch("swanlab.sdk.internal.core_python.transport.sender.console.debug"),
    ):
        sender.upload_save([rec])

    mock_prepare.assert_called_once()
    # 远端 path 仍为记录的相对名
    assert mock_prepare.call_args.kwargs["files"][0]["path"] == "checkpoints/model.txt"
    mock_upload_saves.assert_called_once()
    mock_complete.assert_called_once_with(
        "experiment-id", files=[{"path": "checkpoints/model.txt", "state": "UPLOADED"}]
    )


def test_resolve_save_source_handles_windows_separators_on_posix(tmp_path: Path):
    """跨平台：训练在 Windows 记录（反斜杠），sync 在 POSIX 读取时回退正确。

    - 内部保存 source_path='C:\\host\\config.yaml'：POSIX 的 Path 无法切分反斜杠，
      直接取 .name 会得到整串；必须用 PureWindowsPath 取 basename。
    - 用户保存 name='checkpoints\\model.pt'：files 下真实镜像是分层子目录 trees，
      files_dir / name 会拼成单个含反斜杠的文件名，须用 PureWindowsPath 还原层级。
    """
    sender = _make_sender(tmp_path)

    # ── 内部保存（config）：basename 解析 ──
    config_fallback = tmp_path / "files" / "config.yaml"
    config_fallback.parent.mkdir(parents=True)
    config_fallback.write_text("key: value", encoding="utf-8")
    internal_rec = _make_internal_save_record(r"C:\host\config.yaml", SaveType.SAVE_TYPE_CONFIG, name="config")
    assert resolve_save_source(sender._ctx, internal_rec.save) == config_fallback

    # ── 用户保存（CUSTOM）：层级还原 ──
    custom_fallback = tmp_path / "files" / "checkpoints" / "model.pt"
    custom_fallback.parent.mkdir(parents=True)
    custom_fallback.write_text("weights", encoding="utf-8")
    custom_rec = Record(
        save=SaveRecord(
            name=r"checkpoints\model.pt",
            source_path=r"C:\host\checkpoints\model.pt",
            target_path=r"C:\host\checkpoints\model.pt",
            type=SaveType.SAVE_TYPE_CUSTOM,
        )
    )
    assert resolve_save_source(sender._ctx, custom_rec.save) == custom_fallback


# ============================================================
# skip_store：内部 save payload 双源
# ============================================================


_INTERNAL_PAYLOAD_CASES = [
    pytest.param(
        SaveType.SAVE_TYPE_METADATA,
        b'{"hostname": "gpu-01"}',
        "upload_metadata",
        {"hostname": "gpu-01"},
        id="metadata",
    ),
    pytest.param(
        SaveType.SAVE_TYPE_REQUIREMENTS,
        b"numpy==1.0",
        "upload_requirements",
        "numpy==1.0",
        id="requirements",
    ),
    pytest.param(
        SaveType.SAVE_TYPE_CONFIG,
        b"lr:\n  value: 0.01\n",
        "upload_config",
        {"lr": {"value": 0.01}},
        id="config",
    ),
]


@pytest.mark.parametrize("save_type, payload, api_func, expected_content", _INTERNAL_PAYLOAD_CASES)
def test_upload_save_uses_internal_payload_without_disk(tmp_path, save_type, payload, api_func, expected_content):
    """内部 save payload 非空：直接解析上传，不走 files 回退。"""
    sender = _make_sender(tmp_path)
    record = Record(save=SaveRecord(name="internal", type=save_type, payload=payload))

    with (
        patch(f"swanlab.sdk.internal.core_python.transport.sender.{api_func}") as mock_api,
        patch("swanlab.sdk.internal.core_python.transport.sender.resolve_save_source") as mock_resolve,
    ):
        sender.upload_save([record])

    mock_resolve.assert_not_called()
    mock_api.assert_called_once()
    assert mock_api.call_args.kwargs["content"] == expected_content


def test_upload_save_rejects_invalid_payload_usage(tmp_path: Path):
    """CUSTOM 带 payload、skip_store 下内部 save 缺 payload，均告警跳过、不走上传。"""
    source = tmp_path / "model.pt"
    source.write_bytes(b"weights")
    sender = _make_sender(tmp_path, skip_store=True)
    custom_with_payload = Record(
        save=SaveRecord(
            name="model.pt",
            source_path=str(source),
            type=SaveType.SAVE_TYPE_CUSTOM,
            payload=b"forbidden",
        )
    )
    internal_without_payload = _make_internal_save_record("", SaveType.SAVE_TYPE_METADATA)

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.files.prepare_save_files") as mock_prepare,
        patch("swanlab.sdk.internal.core_python.transport.sender.resolve_save_source") as mock_resolve,
        patch("swanlab.sdk.internal.core_python.transport.sender.upload_metadata") as mock_upload_meta,
        patch("swanlab.sdk.internal.core_python.transport.sender.console.warning") as mock_warning,
    ):
        sender.upload_save([custom_with_payload, internal_without_payload])

    mock_prepare.assert_not_called()
    mock_resolve.assert_not_called()
    mock_upload_meta.assert_not_called()
    messages = [call.args[0] for call in mock_warning.call_args_list]
    assert any("CUSTOM save must not carry payload" in message for message in messages)
    assert any("Internal save payload missing" in message for message in messages)


def test_upload_save_retains_unparseable_payload(tmp_path: Path):
    """payload 无法解析属确定性脏数据：保留记录并报告失败，不进入上传、不冒充完成。"""
    sender = _make_sender(tmp_path)
    record = Record(save=SaveRecord(name="metadata", type=SaveType.SAVE_TYPE_METADATA, payload=b"{not-json"))

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.upload_metadata") as mock_upload_meta,
        patch("swanlab.sdk.internal.core_python.transport.sender.console.warning") as mock_warning,
    ):
        result = sender.upload("save", [record])  # 不得抛异常

    mock_upload_meta.assert_not_called()
    assert result.remaining_records == [record]
    assert result.failures[0].source == "local_data"
    assert "Failed to parse internal save payload" in mock_warning.call_args.args[0]


def test_upload_save_payload_upload_error_is_retained(tmp_path: Path):
    """payload 解析成功但上传失败（5xx）：记录保留在结果中，不静默吞掉。"""
    sender = _make_sender(tmp_path)
    record = Record(
        save=SaveRecord(name="metadata", type=SaveType.SAVE_TYPE_METADATA, payload=b'{"hostname": "gpu-01"}')
    )
    error = ApiError(
        _FakeApiErrorResponse(502),
        method="PUT",
        trace_id="trace-id",
        code="bad-gateway",
        message="server error",
    )

    with patch(
        "swanlab.sdk.internal.core_python.transport.sender.upload_metadata",
        side_effect=error,
    ):
        result = sender.upload("save", [record])

    assert result.remaining_records == [record]
    assert result.failures[0].status == 502


def test_upload_save_disk_channel_upload_error_is_retained(tmp_path: Path):
    """默认模式磁盘回读：上传失败（5xx）同样保留记录，不被静默吞掉。"""
    sender = _make_sender(tmp_path)
    metadata_file = tmp_path / "files" / "swanlab-metadata.json"
    metadata_file.parent.mkdir(parents=True)
    metadata_file.write_text('{"hostname": "gpu-01"}', encoding="utf-8")
    record = _make_internal_save_record("/nonexistent/host/swanlab-metadata.json", SaveType.SAVE_TYPE_METADATA)
    error = ApiError(
        _FakeApiErrorResponse(502),
        method="PUT",
        trace_id="trace-id",
        code="bad-gateway",
        message="server error",
    )

    with patch(
        "swanlab.sdk.internal.core_python.transport.sender.upload_metadata",
        side_effect=error,
    ):
        result = sender.upload("save", [record])

    assert result.remaining_records == [record]
    assert result.failures[0].status == 502


def test_upload_save_warns_once_for_unknown_internal_type(tmp_path: Path):
    """未知内部 save 类型：单次告警跳过，不读盘、不上传。"""
    sender = _make_sender(tmp_path)
    record = Record(save=SaveRecord(name="unknown", type=cast(SaveType, 99), payload=b"opaque"))

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.resolve_save_source") as mock_resolve,
        patch("swanlab.sdk.internal.core_python.transport.sender.console.warning") as mock_warning,
    ):
        sender.upload_save([record])

    mock_resolve.assert_not_called()
    assert "Unknown internal save type" in mock_warning.call_args.args[0]


# ============================================================
# P2：部分成功与数据保留
# ============================================================


def test_upload_column_partial_success_keeps_failed_tail(tmp_path: Path):
    """列上传分批部分成功：前序批次确认，失败批次及后续保留。"""
    sender = _make_sender(tmp_path)
    records = [
        Record(column=ColumnRecord(column_key=f"metric-{index}", column_type=ColumnType.COLUMN_TYPE_SCALAR))
        for index in range(3)
    ]
    error = ApiError(
        _FakeApiErrorResponse(502),
        method="POST",
        trace_id="trace-id",
        code="bad-gateway",
        message="server error",
    )

    with patch(
        "swanlab.sdk.internal.core_python.transport.sender.upload_columns",
        side_effect=[None, error],
    ):
        result = sender.upload_column(records, batch_size=2)

    assert result.confirmed_records == records[:2]
    assert result.remaining_records == records[2:]
    assert result.failures[0].status == 502


def test_upload_media_transfer_failure_retains_all_records(tmp_path: Path):
    """媒体对象存储 PUT 失败：全部参与记录保留，失败来源为对象存储。"""
    media_dir = tmp_path / "media" / "image"
    media_dir.mkdir(parents=True)
    (media_dir / "a.png").write_bytes(b"image")
    sender = _make_sender(tmp_path)
    record = _make_media_record("a.png")

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.upload_resource",
            side_effect=ApiError(
                _FakeApiErrorResponse(503),
                method="PUT",
                trace_id="trace-id",
                code="slow-down",
                message="storage unavailable",
            ),
        ),
    ):
        result = sender.upload_media([record])

    assert result.remaining_records == [record]
    failure = result.failures[0]
    assert failure.source == "object_storage"
    assert failure.stage == "transfer"
    assert failure.status == 503


def test_upload_media_metrics_failure_keeps_file_rows_open(tmp_path: Path):
    """媒体 PUT 成功但指标提交失败：记录保留，文件行不得提前确认完成。"""
    media_dir = tmp_path / "media" / "image"
    media_dir.mkdir(parents=True)
    image_path = media_dir / "a.png"
    image_path.write_bytes(b"image")
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)
    record = _make_media_record("a.png")

    def _fake_upload_resource(_, __, *, paths, buffers, content_types=None, tracker=None):
        # 模拟 PUT 全部字节传输完成（finish 留给指标提交成功后的确认阶段）
        if tracker is not None:
            for path in paths:
                tracker.update_file_progress(f"{path}:5", f"{path}:5", 5)

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.upload_resource",
            side_effect=_fake_upload_resource,
        ),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.upload_media",
            side_effect=ApiError(
                _FakeApiErrorResponse(500),
                method="POST",
                trace_id="trace-id",
                code="internal-error",
                message="server error",
            ),
        ),
    ):
        result = sender.upload_media([record])

    assert result.remaining_records == [record]
    assert result.failures[0].stage == "submit_metrics"
    stats = tracker.snapshot()
    # 字节进度已全部上报，但文件完成计数为 0：finish_file 只在指标提交成功后调用，
    # 旧的"PUT 成功即确认"行为会让 uploaded_number 提前变为 1
    assert stats.uploaded_size == 5
    assert stats.uploaded_number == 0


def test_upload_media_confirms_records_after_metrics_success(tmp_path: Path):
    """媒体内容与指标均成功：文件行完成，Record 进入 confirmed。"""
    media_dir = tmp_path / "media" / "image"
    media_dir.mkdir(parents=True)
    image_path = media_dir / "a.png"
    image_path.write_bytes(b"image")
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)
    record = _make_media_record("a.png")

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch("swanlab.sdk.internal.core_python.transport.sender.files.upload_resource"),
        patch("swanlab.sdk.internal.core_python.transport.sender.upload_media"),
    ):
        result = sender.upload_media([record])

    assert result.ok
    assert result.confirmed_records == [record]
    stats = tracker.snapshot()
    assert stats.uploaded_size == 5
    assert stats.uploaded_number == 1
    assert list(stats.files) == []


def test_upload_save_small_files_partial_success_completes_only_uploaded(tmp_path: Path):
    """两个小文件其中一个 PUT 失败：complete 只含成功文件，失败记录保留。"""
    first = tmp_path / "ok.txt"
    second = tmp_path / "failed.txt"
    first.write_text("weights-1", encoding="utf-8")
    second.write_text("weights-2", encoding="utf-8")
    sender = _make_sender(tmp_path)

    def _fake_upload(session, *, resources, tracker=None):
        return {1: RuntimeError("s3 failed")}  # 第二个文件失败

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.files.prepare_save_files",
            side_effect=lambda _, *, files: {"urls": [f"https://s3.test/{item['path']}" for item in files]},
        ),
        patch("swanlab.sdk.internal.core_python.transport.sender.files.complete_save_files") as mock_complete,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch("swanlab.sdk.internal.core_python.transport.sender.files.upload_saves", side_effect=_fake_upload),
    ):
        records = [
            _make_save_record(first, name="checkpoints/ok.txt"),
            _make_save_record(second, name="checkpoints/failed.txt"),
        ]
        result = sender.upload_save(records)

    # complete 只上报成功的文件
    assert mock_complete.call_args.kwargs["files"] == [{"path": "checkpoints/ok.txt", "state": "UPLOADED"}]
    assert result.confirmed_records == [records[0]]
    assert result.remaining_records == [records[1]]
    assert result.failures[0].source == "object_storage"


def test_upload_save_oversize_file_is_retained_not_uploaded(tmp_path: Path):
    """超过 save_size 限制的文件：明确拒绝结果，记录保留，不伪报上传成功。"""
    source = tmp_path / "huge.bin"
    source.write_bytes(b"0" * 16)
    sender = _make_sender(tmp_path, save_split=100 * 1024 * 1024, save_size=8)  # 模拟超限
    record = _make_save_record(source, name="huge.bin")

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.files.prepare_save_files") as mock_prepare,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
    ):
        result = sender.upload_save([record])

    mock_prepare.assert_not_called()
    assert result.remaining_records == [record]
    assert result.failures[0].source == "local_data"
    assert "size limit" in result.failures[0].message


def test_upload_save_invalid_save_batch_retains_records(tmp_path: Path):
    """save_batch 配置非法：整批保留并报告，不发起请求。"""
    source = tmp_path / "model.txt"
    source.write_text("weights", encoding="utf-8")
    sender = _make_sender(tmp_path, save_batch=0)
    record = _make_save_record(source, name="checkpoints/model.txt")

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.files.prepare_save_files") as mock_prepare,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch("swanlab.sdk.internal.core_python.transport.sender.console.warning"),
    ):
        result = sender.upload_save([record])

    mock_prepare.assert_not_called()
    assert result.remaining_records == [record]
    assert "save batch size" in result.failures[0].message


def test_upload_save_missing_source_file_is_retained(tmp_path: Path):
    """源文件缺失：记录保留并报告 local_data 失败，不静默过滤。"""
    sender = _make_sender(tmp_path)
    record = _make_save_record(Path("/nonexistent/training-host/model.txt"), name="checkpoints/model.txt")

    with (
        patch("swanlab.sdk.internal.core_python.transport.sender.files.prepare_save_files") as mock_prepare,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch("swanlab.sdk.internal.core_python.transport.sender.console.warning"),
    ):
        result = sender.upload_save([record])

    mock_prepare.assert_not_called()
    assert result.remaining_records == [record]
    assert result.failures[0].source == "local_data"
    assert "not found" in result.failures[0].message


# ============================================================
# encode_column: presence-based encoding
# ============================================================


class TestEncodeColumn:
    """encode_column 的 presence-based 编码验证。"""

    def test_inferred_column_only_section_name(self):
        """inferred column 只发 sectionName，不发 xAxis/hidden。"""
        col = ColumnRecord(
            column_key="train/loss",
            column_type=ColumnType.COLUMN_TYPE_SCALAR,
            section_name="Training",
        )
        result = encode_column(col)
        assert result == {"key": "train/loss", "type": "FLOAT", "sectionName": "Training"}

    def test_inferred_column_no_section_name(self):
        """inferred column 无 sectionName 时只发 key/type。"""
        col = ColumnRecord(
            column_key="train/loss",
            column_type=ColumnType.COLUMN_TYPE_SCALAR,
        )
        result = encode_column(col)
        assert result == {"key": "train/loss", "type": "FLOAT"}

    def test_column_with_custom_xaxis(self):
        """ColumnRecord 带 x_axis → xAxis 发送。"""
        col = ColumnRecord(
            column_key="train/loss",
            column_type=ColumnType.COLUMN_TYPE_SCALAR,
            x_axis="train/epoch",
            section_name="Training",
        )
        result = encode_column(col)
        assert result == {
            "key": "train/loss",
            "type": "FLOAT",
            "xAxis": "train/epoch",
            "sectionName": "Training",
        }

    def test_column_without_xaxis(self):
        """ColumnRecord 无 x_axis → xAxis key 不发送。"""
        col = ColumnRecord(
            column_key="loss",
            column_type=ColumnType.COLUMN_TYPE_SCALAR,
            section_name="",
        )
        result = encode_column(col)
        # section_name="" → falsy → 不发送
        assert result == {"key": "loss", "type": "FLOAT"}

    def test_column_with_relative_time_xaxis(self):
        """xAxis=_relative_time → 字符串透传，由 Server 转换。"""
        col = ColumnRecord(
            column_key="loss",
            column_type=ColumnType.COLUMN_TYPE_SCALAR,
            x_axis="_relative_time",
        )
        result = encode_column(col)
        assert result == {
            "key": "loss",
            "type": "FLOAT",
            "xAxis": "_relative_time",
        }

    def test_column_hidden_true(self):
        """hidden=True → hidden key 发送。"""
        col = ColumnRecord(
            column_key="secret",
            column_type=ColumnType.COLUMN_TYPE_SCALAR,
            hidden=True,
        )
        result = encode_column(col)
        assert result == {
            "key": "secret",
            "type": "FLOAT",
            "hidden": True,
        }

    def test_column_hidden_false_not_sent(self):
        """hidden=False → hidden key 不发送（server 视 absent 为 false）。"""
        col = ColumnRecord(
            column_key="loss",
            column_type=ColumnType.COLUMN_TYPE_SCALAR,
            hidden=False,
        )
        result = encode_column(col)
        assert result == {"key": "loss", "type": "FLOAT"}

    def test_media_column(self):
        """media ColumnRecord 无 x_axis（resolver 对 MEDIA 不设 x_axis）。"""
        col = ColumnRecord(
            column_key="train/img",
            column_type=ColumnType.COLUMN_TYPE_IMAGE,
            section_name="Images",
        )
        result = encode_column(col)
        assert result == {
            "key": "train/img",
            "type": "IMAGE",
            "sectionName": "Images",
        }

    def test_system_class_overrides_type(self):
        """SYSTEM class → type 被覆写为 'SYSTEM'。"""
        col = ColumnRecord(
            column_key="system/cpu",
            column_type=ColumnType.COLUMN_TYPE_SCALAR,
            column_class=ColumnClass.COLUMN_CLASS_SYSTEM,
        )
        result = encode_column(col)
        assert result == {"key": "system/cpu", "type": "SYSTEM"}

    def test_upload_column_end_to_end(self, tmp_path: Path):
        """upload_column 端到端：带 xAxis/hidden 的 ColumnRecord 正确编码并发送到 HTTP。"""
        sender = _make_sender(tmp_path)
        records = [
            Record(
                column=ColumnRecord(
                    column_key="train/loss",
                    column_type=ColumnType.COLUMN_TYPE_SCALAR,
                    x_axis="train/epoch",
                    section_name="Training",
                    hidden=True,
                )
            )
        ]
        with patch("swanlab.sdk.internal.core_python.transport.sender.upload_columns") as mock_upload:
            sender.upload_column(records)
        assert mock_upload.call_count == 1
        series = mock_upload.call_args.kwargs["columns"]["series"]
        assert len(series) == 1
        assert series[0] == {
            "key": "train/loss",
            "type": "FLOAT",
            "xAxis": "train/epoch",
            "sectionName": "Training",
            "hidden": True,
        }
