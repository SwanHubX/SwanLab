from pathlib import Path
from typing import cast

from requests.sessions import Session

from swanlab.sdk.internal.core_python.api import upload as upload_api


class _FakeResponse:
    def raise_for_status(self) -> None:
        pass


class _FakeSession:
    def __init__(self) -> None:
        self.puts = []

    def put(self, url, *, data, headers):
        self.puts.append(
            {
                "url": url,
                "data": data.read(),
                "headers": headers,
            }
        )
        return _FakeResponse()


def test_upload_resources_uploads_files_and_reports_progress(tmp_path: Path):
    first = tmp_path / "first.txt"
    second = tmp_path / "second.bin"
    first.write_text("hello", encoding="utf-8")
    second.write_bytes(b"\x00\x01")
    fake_session = _FakeSession()
    session = cast(Session, fake_session)
    uploaded = upload_api.upload_saves(
        session,
        resources=[
            {"url": "https://s3.test/first", "source_path": str(first), "content_type": "text/plain"},
            {
                "url": "https://s3.test/second",
                "source_path": str(second),
                "content_type": "application/octet-stream",
            },
        ],
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
