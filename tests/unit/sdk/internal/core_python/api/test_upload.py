from io import BytesIO
from pathlib import Path
from typing import Optional, cast

import pytest
from requests.sessions import Session

from swanlab.sdk.internal.core_python.api import upload as upload_api
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker


class _FakeResponse:
    def __init__(self, *, should_fail: bool = False) -> None:
        self._should_fail = should_fail

    def raise_for_status(self) -> None:
        if self._should_fail:
            raise RuntimeError("put failed")
        pass


class _FakeSession:
    def __init__(self, failing_urls: Optional[set[str]] = None) -> None:
        self._failing_urls = failing_urls or set()
        self.puts = []

    def put(self, url, *, data, headers):
        if hasattr(data, "read"):
            payload = data.read()
        else:
            payload = data.getvalue()
        self.puts.append(
            {
                "url": url,
                "data": payload,
                "headers": headers,
            }
        )
        return _FakeResponse(should_fail=url in self._failing_urls)


class _FakePresignResponse:
    def __init__(self, urls: list[str]) -> None:
        self.data = {"urls": urls}


class _RetryingReadSession:
    def put(self, url, *, data, headers):
        data.read()
        data.seek(0)
        data.read(1)
        data.read()
        return _FakeResponse()


def test_upload_resources_uploads_files_and_reports_progress(tmp_path: Path):
    first = tmp_path / "first.txt"
    second = tmp_path / "second.bin"
    first.write_text("hello", encoding="utf-8")
    second.write_bytes(b"\x00\x01")
    fake_session = _FakeSession()
    session = cast(Session, fake_session)
    tracker = UploadTracker()

    # 注册文件到 tracker（模拟 sender 侧 _track_file 的行为）
    tracker.add_file(key="first-key", path=str(first), total=5)
    tracker.add_file(key="second-key", path=str(second), total=2)

    uploaded = upload_api.upload_saves(
        session,
        resources=[
            {
                "url": "https://s3.test/first",
                "source_path": str(first),
                "content_type": "text/plain",
                "tracker_key": "first-key",
                "size": 5,
            },
            {
                "url": "https://s3.test/second",
                "source_path": str(second),
                "content_type": "application/octet-stream",
                "tracker_key": "second-key",
                "size": 2,
            },
        ],
        tracker=tracker,
    )

    puts_by_url = {item["url"]: item for item in fake_session.puts}
    assert puts_by_url == {
        "https://s3.test/first": {
            "url": "https://s3.test/first",
            "data": b"hello",
            "headers": {"Content-Type": "text/plain"},
        },
        "https://s3.test/second": {
            "url": "https://s3.test/second",
            "data": b"\x00\x01",
            "headers": {"Content-Type": "application/octet-stream"},
        },
    }
    assert uploaded == {str(first), str(second)}
    # tracker 应记录两个已完成文件
    stats = tracker.snapshot()
    assert stats.uploaded_number == 2
    assert stats.uploaded_size >= 7  # 5 + 2


def test_upload_saves_does_not_report_progress_when_put_fails(tmp_path: Path):
    source = tmp_path / "failed.bin"
    source.write_bytes(b"failed")
    fake_session = _FakeSession(failing_urls={"https://s3.test/failed"})
    tracker = UploadTracker()
    tracker.add_file(key="failed-key", path=str(source), total=6)

    uploaded = upload_api.upload_saves(
        cast(Session, fake_session),
        resources=[
            {
                "url": "https://s3.test/failed",
                "source_path": str(source),
                "content_type": "application/octet-stream",
                "tracker_key": "failed-key",
                "size": 6,
            }
        ],
        tracker=tracker,
    )

    assert uploaded == set()
    # tracker 不应记录完成
    stats = tracker.snapshot()
    assert stats.uploaded_number == 0


def test_upload_resource_reports_progress_after_successful_put():
    fake_session = _FakeSession()
    tracker = UploadTracker()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            upload_api.client,
            "post",
            lambda *_args, **_kwargs: _FakePresignResponse(["https://s3.test/image", "https://s3.test/text"]),
        )

        # 模拟 sender._track_file：上传前注册文件
        tracker.add_file(key="media/image/a.png:5", path="media/image/a.png", total=5)
        tracker.add_file(key="media/text/a.txt:4", path="media/text/a.txt", total=4)

        upload_api.upload_resource(
            cast(Session, fake_session),
            "experiment-id",
            paths=["media/image/a.png", "media/text/a.txt"],
            buffers=[BytesIO(b"image"), BytesIO(b"text")],
            tracker=tracker,
        )

    # tracker 应记录两个已完成文件（由 upload_resource 内部 finish_file）
    stats = tracker.snapshot()
    assert stats.uploaded_number == 2
    assert stats.uploaded_size >= 9  # 5 + 4


def test_upload_resource_uses_provided_content_type(tmp_path: Path):
    html_file = tmp_path / "report.html"
    html_file.write_text("<h1>report</h1>", encoding="utf-8")
    fake_session = _FakeSession()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            upload_api.client,
            "post",
            lambda *_args, **_kwargs: _FakePresignResponse(["https://s3.test/html"]),
        )

        upload_api.upload_resource(
            cast(Session, fake_session),
            "experiment-id",
            paths=["media/html/report.html"],
            buffers=[html_file],
            content_types=["text/html"],
        )

    assert fake_session.puts == [
        {
            "url": "https://s3.test/html",
            "data": b"<h1>report</h1>",
            "headers": {"Content-Type": "text/html"},
        }
    ]


def test_upload_resource_handles_standard_file_objects(tmp_path: Path):
    temp_file = tmp_path / "media_file.bin"
    temp_file.write_bytes(b"some media data")
    fake_session = _FakeSession()
    tracker = UploadTracker()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            upload_api.client,
            "post",
            lambda *_args, **_kwargs: _FakePresignResponse(["https://s3.test/media"]),
        )

        tracker.add_file(key="media/video/media_file.bin:15", path=str(temp_file), total=15)

        with open(temp_file, "rb") as f:
            upload_api.upload_resource(
                cast(Session, fake_session),
                "experiment-id",
                paths=["media/video/media_file.bin"],
                buffers=[f],
                tracker=tracker,
            )

    stats = tracker.snapshot()
    assert stats.uploaded_number == 1
    assert stats.uploaded_size == 15
    assert len(fake_session.puts) == 1
    assert fake_session.puts[0]["data"] == b"some media data"


def test_upload_resource_read_progress_stays_within_buffer_size_on_seek():
    """验证 ProgressFileWrapper 在 seek/retry 场景下进度不会超过 total_size。"""
    tracker = UploadTracker()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            upload_api.client,
            "post",
            lambda *_args, **_kwargs: _FakePresignResponse(["https://s3.test/media"]),
        )

        tracker.add_file(key="media/image/a.png:5", path="media/image/a.png", total=5)

        upload_api.upload_resource(
            cast(Session, _RetryingReadSession()),
            "experiment-id",
            paths=["media/image/a.png"],
            buffers=[BytesIO(b"image")],
            tracker=tracker,
        )

    stats = tracker.snapshot()
    # 文件大小 5，不应该超过
    assert stats.uploaded_size == 5
    assert stats.uploaded_number == 1
