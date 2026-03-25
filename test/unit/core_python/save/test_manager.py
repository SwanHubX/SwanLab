import hashlib
from unittest.mock import MagicMock

import swanlab.core_python.save.manager as save_manager
from swanlab.core_python.save import SaveFile


class _FakeFuture:
    def __init__(self, result=None, error=None):
        self._result = result
        self._error = error

    def result(self):
        if self._error is not None:
            raise self._error
        return self._result


class _FakeExecutor:
    def __init__(self, max_workers=None, thread_name_prefix=None):
        self.max_workers = max_workers
        self.thread_name_prefix = thread_name_prefix
        self.submits = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return False

    def submit(self, fn, *args, **kwargs):
        self.submits.append((fn, args, kwargs))
        try:
            return _FakeFuture(result=fn(*args, **kwargs))
        except Exception as e:
            return _FakeFuture(error=e)

    def shutdown(self, wait=True):
        return None


def test_upload_single_part_passes_md5_and_mime_type(tmp_path, monkeypatch):
    source = tmp_path / "metrics.txt"
    source.write_text("hello world")
    file = SaveFile(
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

    def fake_complete_upload(client, exp_id, names, state="UPLOADED"):
        captured["complete"] = (client, exp_id, names, state)

    monkeypatch.setattr(save_manager, "prepare_upload", fake_prepare_upload)
    monkeypatch.setattr(save_manager, "upload_file", fake_upload_file)
    monkeypatch.setattr(save_manager, "complete_upload", fake_complete_upload)

    try:
        manager._upload_single_part(file, source.stat().st_size, "exp-123")
    finally:
        manager.close()

    assert captured["prepare"][1] == "exp-123"
    assert captured["prepare"][2] == [
        {
            "name": "logs/metrics.txt",
            "size": source.stat().st_size,
            "md5": hashlib.md5(b"hello world").hexdigest(),
            "mimeType": "text/plain",
        }
    ]
    assert captured["upload"] == ("https://upload.example.com/file", b"hello world", 3)
    assert captured["complete"][1:] == ("exp-123", ["logs/metrics.txt"], "UPLOADED")


def test_upload_multipart_passes_md5_count_upload_id_and_mime_type(tmp_path, monkeypatch):
    content = b"abcdefghij"
    source = tmp_path / "checkpoint.txt"
    source.write_bytes(content)
    file = SaveFile(
        source_path=str(source),
        name="checkpoints/checkpoint.txt",
        target_path=str(tmp_path / "run" / "checkpoints" / "checkpoint.txt"),
    )
    manager = save_manager.FileUploadManager(mode="disabled", file_dir=str(tmp_path))
    manager._client = MagicMock()

    captured = {"parts": []}

    def fake_prepare_multipart(client, exp_id, name, size, part_count, md5, mime_type=None):
        captured["prepare"] = (client, exp_id, name, size, part_count, md5, mime_type)
        return {
            "uploadId": "upload-1",
            "partSize": 4,
            "uploadUrls": [
                "https://upload.example.com/part-1",
                "https://upload.example.com/part-2",
                "https://upload.example.com/part-3",
            ],
        }

    fake_executor = _FakeExecutor()

    def fake_executor_factory(max_workers=None, thread_name_prefix=None):
        fake_executor.max_workers = max_workers
        fake_executor.thread_name_prefix = thread_name_prefix
        return fake_executor

    def fake_upload_file(*, url, buffer, max_retries=3):
        captured["parts"].append((url, buffer.read(), max_retries))

    def fake_complete_multipart(client, exp_id, name, upload_id, state="UPLOADED"):
        captured["complete"] = (client, exp_id, name, upload_id, state)

    monkeypatch.setattr(save_manager, "PART_SIZE", 4)
    monkeypatch.setattr(save_manager, "prepare_multipart", fake_prepare_multipart)
    monkeypatch.setattr(save_manager, "ThreadPoolExecutor", fake_executor_factory)
    monkeypatch.setattr(save_manager, "upload_file", fake_upload_file)
    monkeypatch.setattr(save_manager, "complete_multipart", fake_complete_multipart)

    try:
        manager._upload_multipart(file, len(content), "exp-123")
    finally:
        manager.close()

    assert captured["prepare"][1:] == (
        "exp-123",
        "checkpoints/checkpoint.txt",
        len(content),
        3,
        hashlib.md5(content).hexdigest(),
        "text/plain",
    )
    assert fake_executor.max_workers == 3
    assert captured["parts"] == [
        ("https://upload.example.com/part-1", b"abcd", 3),
        ("https://upload.example.com/part-2", b"efgh", 3),
        ("https://upload.example.com/part-3", b"ij", 3),
    ]
    assert captured["complete"][1:] == (
        "exp-123",
        "checkpoints/checkpoint.txt",
        "upload-1",
        "UPLOADED",
    )
