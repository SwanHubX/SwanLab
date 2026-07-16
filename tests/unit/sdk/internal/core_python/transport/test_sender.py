from io import BytesIO
from pathlib import Path
from unittest.mock import ANY, MagicMock, patch

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord, ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem, MediaRecord, MediaValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord, SaveType
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


def _make_media_record(filename: str, media_type=ColumnType.COLUMN_TYPE_IMAGE) -> Record:
    timestamp = Timestamp()
    timestamp.GetCurrentTime()
    return Record(
        media=MediaRecord(
            key="examples/image",
            step=1,
            type=media_type,
            timestamp=timestamp,
            value=MediaValue(items=[MediaItem(filename=filename)]),
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
            "swanlab.sdk.internal.core_python.transport.sender.upload_resource",
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


def test_upload_advances_records_only_on_success(tmp_path: Path):
    """成功上传才递进 uploaded，进度条仅保留上传成功的进度。"""
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)
    sender._upload_handlers["scalar"] = MagicMock()

    sender.upload("scalar", [Record(), Record()])

    assert tracker.snapshot().uploaded_records == 2


def test_upload_does_not_advance_records_for_non_retryable_api_error(tmp_path: Path):
    """4xx 业务拒绝不计入 uploaded：失败原地不动，进度条不递进（也不抛异常，避免上层无意义重试）。"""
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

    assert tracker.snapshot().uploaded_records == 0


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
        result = sender._resolve_save_source(rec.save)

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
        result = sender._resolve_save_source(rec.save)

    assert result == fallback
    mock_debug.assert_called_once()


def test_resolve_save_source_falls_back_for_internal_metadata(tmp_path: Path):
    """内部 metadata：source_path 不可读 → 按 basename 回退到 files_dir/swanlab-metadata.json。"""
    rec = _make_internal_save_record("/nonexistent/host/swanlab-metadata.json", SaveType.SAVE_TYPE_METADATA)
    sender = _make_sender(tmp_path)
    fallback = tmp_path / "files" / "swanlab-metadata.json"
    fallback.parent.mkdir(parents=True)
    fallback.write_text("{}", encoding="utf-8")

    result = sender._resolve_save_source(rec.save)

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

    result = sender._resolve_save_source(rec.save)

    assert result == fallback


def test_resolve_save_source_returns_none_when_both_unreadable(tmp_path: Path):
    rec = _make_save_record(Path("/nonexistent/training-host/model.txt"), name="checkpoints/model.txt")
    sender = _make_sender(tmp_path)
    # 既无 primary 也无 files 镜像

    assert sender._resolve_save_source(rec.save) is None


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
        return {r["source_path"] for r in resources}

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.prepare_save_files",
            return_value={"urls": ["https://s3.test/model"]},
        ) as mock_prepare,
        patch("swanlab.sdk.internal.core_python.transport.sender.complete_save_files") as mock_complete,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.upload_saves", side_effect=_fake_upload
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
    assert sender._resolve_save_source(internal_rec.save) == config_fallback

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
    assert sender._resolve_save_source(custom_rec.save) == custom_fallback
