import hashlib
from unittest.mock import MagicMock

import swanlab.core_python.save.manager as save_manager
from swanlab.core_python.save import SaveFileModel


def test_upload_single_passes_md5_and_mime_type(tmp_path, monkeypatch):
    source = tmp_path / "metrics.txt"
    source.write_text("hello world")
    file = SaveFileModel(
        source_path=str(source),
        name="logs/metrics.txt",
        target_path=str(tmp_path / "run" / "logs" / "metrics.txt"),
    )
    manager = save_manager.FileUploadManager(mode="disabled", file_dir=str(tmp_path))
    manager._client = MagicMock()

    captured = {}

    def fake_prepare_upload(client, exp_id, files):
        captured["prepare"] = (client, exp_id, files)
        return ["https://upload.example.com/file"]

    def fake_upload_file(*, url, buffer, max_retries=3, mime_type=None):
        captured["upload"] = (url, buffer.read(), max_retries, mime_type)

    def fake_complete_upload(client, exp_id, files):
        captured["complete"] = (client, exp_id, files)

    monkeypatch.setattr(save_manager, "prepare_upload", fake_prepare_upload)
    monkeypatch.setattr(save_manager, "upload_file", fake_upload_file)
    monkeypatch.setattr(save_manager, "complete_upload", fake_complete_upload)

    try:
        manager._upload_single(file, source.stat().st_size, "exp-123")
    finally:
        manager.close()

    assert captured["prepare"][1] == "exp-123"
    assert captured["prepare"][2] == [
        {
            "path": "logs/metrics.txt",
            "size": source.stat().st_size,
            "md5": hashlib.md5(b"hello world").hexdigest(),
            "mimeType": "text/plain",
        }
    ]
    assert captured["upload"] == (
        "https://upload.example.com/file",
        b"hello world",
        3,
        "text/plain",
    )
    assert captured["complete"][1:] == (
        "exp-123",
        [{"path": "logs/metrics.txt", "state": "UPLOADED"}],
    )


def test_upload_multipart_passes_md5_and_upload_id(tmp_path, monkeypatch):
    content = b"abcdefghij"
    source = tmp_path / "checkpoint.txt"
    source.write_bytes(content)
    file = SaveFileModel(
        source_path=str(source),
        name="checkpoints/checkpoint.txt",
        target_path=str(tmp_path / "run" / "checkpoints" / "checkpoint.txt"),
    )
    manager = save_manager.FileUploadManager(mode="disabled", file_dir=str(tmp_path))
    manager._client = MagicMock()

    captured = {}

    def fake_prepare_multipart(client, exp_id, file):
        captured["prepare"] = {"exp_id": exp_id, "file": file}
        return {
            "uploadId": "upload-1",
            "partSize": 4,
            "urls": [
                {"partNumber": 1, "url": "https://upload.example.com/part-1"},
                {"partNumber": 2, "url": "https://upload.example.com/part-2"},
                {"partNumber": 3, "url": "https://upload.example.com/part-3"},
            ],
        }

    def fake_upload_file(*, url, buffer, max_retries=3):
        if "parts" not in captured:
            captured["parts"] = []
        captured["parts"].append((url, buffer.read()))
        return {
            "https://upload.example.com/part-1": '"etag-1"',
            "https://upload.example.com/part-2": '"etag-2"',
            "https://upload.example.com/part-3": '"etag-3"',
        }[url]

    def fake_complete_multipart(client, exp_id, file):
        captured["complete"] = {"exp_id": exp_id, "file": file}

    monkeypatch.setattr(save_manager, "PART_SIZE", 4)
    monkeypatch.setattr(save_manager, "prepare_multipart", fake_prepare_multipart)
    monkeypatch.setattr(save_manager, "upload_file", fake_upload_file)
    monkeypatch.setattr(save_manager, "complete_multipart", fake_complete_multipart)

    try:
        manager._upload_multipart(file, len(content), "exp-123")
    finally:
        manager.close()

    assert captured["prepare"]["exp_id"] == "exp-123"
    assert captured["prepare"]["file"] == {
        "path": "checkpoints/checkpoint.txt",
        "size": len(content),
        "md5": hashlib.md5(content).hexdigest(),
        "mimeType": "text/plain",
        "count": 3,
    }
    assert captured["parts"] == [
        ("https://upload.example.com/part-1", b"abcd"),
        ("https://upload.example.com/part-2", b"efgh"),
        ("https://upload.example.com/part-3", b"ij"),
    ]
    assert captured["complete"]["exp_id"] == "exp-123"
    assert captured["complete"]["file"] == {
        "path": "checkpoints/checkpoint.txt",
        "uploadId": "upload-1",
        "parts": [
            {"partNumber": 1, "etag": "etag-1"},
            {"partNumber": 2, "etag": "etag-2"},
            {"partNumber": 3, "etag": "etag-3"},
        ],
    }


def test_upload_single_failure_reports_failed_state(tmp_path, monkeypatch):
    source = tmp_path / "broken.txt"
    source.write_text("broken")
    file = SaveFileModel(
        source_path=str(source),
        name="logs/broken.txt",
        target_path=str(tmp_path / "run" / "logs" / "broken.txt"),
    )
    manager = save_manager.FileUploadManager(mode="disabled", file_dir=str(tmp_path))
    manager._client = MagicMock()

    captured = {}

    def fake_prepare_upload(client, exp_id, files):
        return ["https://upload.example.com/broken"]

    def fake_upload_file(*, url, buffer, max_retries=3, mime_type=None):
        raise RuntimeError("Simulated upload failure")

    def fake_complete_upload(client, exp_id, files):
        captured["complete"] = files

    monkeypatch.setattr(save_manager, "prepare_upload", fake_prepare_upload)
    monkeypatch.setattr(save_manager, "upload_file", fake_upload_file)
    monkeypatch.setattr(save_manager, "complete_upload", fake_complete_upload)

    try:
        manager._upload_single(file, source.stat().st_size, "exp-456")
    finally:
        manager.close()

    assert captured["complete"] == [{"path": "logs/broken.txt", "state": "FAILED"}]
