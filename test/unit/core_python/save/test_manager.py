import hashlib
from unittest.mock import MagicMock

import swanlab.core_python.save.manager as save_manager
from swanlab.core_python.save import SaveFilePayload, SaveFileState, WatchSaveFileModel


def test_upload_single_passes_md5_and_mime_type(tmp_path, monkeypatch):
    source = tmp_path / "metrics.txt"
    source.write_text("hello world")
    file = WatchSaveFileModel(
        source_path=str(source),
        name="logs/metrics.txt",
        target_path=str(tmp_path / "run" / "logs" / "metrics.txt"),
    )
    manager = save_manager.FileUploadManager(mode="disabled", file_dir=str(tmp_path))
    manager._client = MagicMock()

    captured = {}

    def fake_prepare_upload(client, exp_id, files):
        captured["prepare"] = (client, exp_id, files)
        return [{"uploadUrl": "https://upload.example.com/file"}]

    def fake_upload_file(*, url, buffer, max_retries=3):
        captured["upload"] = (url, buffer.read(), max_retries)

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
        SaveFilePayload(
            name="logs/metrics.txt",
            size=source.stat().st_size,
            md5=hashlib.md5(b"hello world").hexdigest(),
            mime_type="text/plain",
        )
    ]
    assert captured["upload"] == ("https://upload.example.com/file", b"hello world", 3)
    assert captured["complete"][1:] == (
        "exp-123",
        [SaveFilePayload(name="logs/metrics.txt", state=SaveFileState.UPLOADED)],
    )


def test_upload_multipart_passes_md5_and_upload_id(tmp_path, monkeypatch):
    content = b"abcdefghij"
    source = tmp_path / "checkpoint.txt"
    source.write_bytes(content)
    file = WatchSaveFileModel(
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
            "uploadUrls": [
                "https://upload.example.com/part-1",
                "https://upload.example.com/part-2",
                "https://upload.example.com/part-3",
            ],
        }

    def fake_upload_file(*, url, buffer, max_retries=3):
        if "parts" not in captured:
            captured["parts"] = []
        captured["parts"].append((url, buffer.read()))

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
    assert captured["prepare"]["file"] == SaveFilePayload(
        name="checkpoints/checkpoint.txt",
        size=len(content),
        md5=hashlib.md5(content).hexdigest(),
        mime_type="text/plain",
        count=3,
    )
    assert captured["parts"] == [
        ("https://upload.example.com/part-1", b"abcd"),
        ("https://upload.example.com/part-2", b"efgh"),
        ("https://upload.example.com/part-3", b"ij"),
    ]
    assert captured["complete"]["exp_id"] == "exp-123"
    assert captured["complete"]["file"] == SaveFilePayload(
        name="checkpoints/checkpoint.txt",
        upload_id="upload-1",
        state=SaveFileState.UPLOADED,
    )
