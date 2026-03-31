import os
import hashlib
from unittest.mock import MagicMock

import swanlab.core_python.save.manager as save_manager
import swanlab.core_python.save.progress as save_progress
from swanlab.core_python.save import SaveFileModel


def test_iter_files_always_returns_list(tmp_path):
    """_iter_files should normalize both single items and iterables to a list."""
    file = SaveFileModel(
        source_path=str(tmp_path / "metrics.txt"),
        name="logs/metrics.txt",
        target_path=str(tmp_path / "run" / "logs" / "metrics.txt"),
    )

    single = save_manager._iter_files(file)
    multiple = save_manager._iter_files(item for item in [file])

    assert isinstance(single, list)
    assert isinstance(multiple, list)
    assert single == [file]
    assert multiple == [file]


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
            "parts": [
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


def test_do_upload_skips_file_exceeding_size_limit(tmp_path, monkeypatch):
    """文件超过 MAX_FILE_SIZE 时应跳过上传并 warning，不抛异常"""
    source = tmp_path / "huge.bin"
    source.write_text("data")
    file = SaveFileModel(
        source_path=str(source),
        name="data/huge.bin",
        target_path=str(tmp_path / "run" / "data" / "huge.bin"),
    )
    manager = save_manager.FileUploadManager(mode="local", file_dir=str(tmp_path))

    # 将阈值降到比文件实际大小更小以触发检查
    monkeypatch.setattr(save_manager, "MAX_FILE_SIZE", 1)

    warnings = []

    def fake_warning(msg):
        warnings.append(msg)

    monkeypatch.setattr(save_manager.swanlog, "warning", fake_warning)

    try:
        manager._do_upload(file)
    finally:
        manager.close()

    assert len(warnings) == 1
    assert "exceeds the size limit" in warnings[0]
    assert "data/huge.bin" in warnings[0]


def test_do_upload_skips_filesystem_checks_in_disabled_mode(tmp_path, monkeypatch):
    """Disabled mode should short-circuit before touching the filesystem."""
    file = SaveFileModel(
        source_path=str(tmp_path / "missing.bin"),
        name="data/missing.bin",
        target_path=str(tmp_path / "run" / "data" / "missing.bin"),
    )
    manager = save_manager.FileUploadManager(mode="disabled", file_dir=str(tmp_path))

    def raise_exists(path):
        raise AssertionError("os.path.exists should not be called in disabled mode")

    def raise_getsize(path):
        raise AssertionError("os.path.getsize should not be called in disabled mode")

    monkeypatch.setattr(save_manager.os.path, "exists", raise_exists)
    monkeypatch.setattr(save_manager.os.path, "getsize", raise_getsize)

    try:
        manager._do_upload(file)
    finally:
        manager.close()


def test_do_upload_skips_file_at_exact_size_limit(tmp_path, monkeypatch):
    """文件恰好等于 MAX_FILE_SIZE 时不触发阈值（仅 > 时才跳过）"""
    source = tmp_path / "exact.bin"
    source.write_bytes(b"\x00" * 10)
    file = SaveFileModel(
        source_path=str(source),
        name="data/exact.bin",
        target_path=str(tmp_path / "run" / "data" / "exact.bin"),
    )
    monkeypatch.setattr(save_manager, "get_client", lambda: MagicMock())
    manager = save_manager.FileUploadManager(mode="cloud", file_dir=str(tmp_path))

    monkeypatch.setattr(save_manager, "MAX_FILE_SIZE", 10)

    warnings = []

    def fake_warning(msg):
        warnings.append(msg)

    monkeypatch.setattr(save_manager.swanlog, "warning", fake_warning)

    fake_upload_cloud = MagicMock()
    monkeypatch.setattr(manager, "_upload_cloud", fake_upload_cloud)

    try:
        manager._do_upload(file)
    finally:
        manager.close()

    assert len(warnings) == 0
    fake_upload_cloud.assert_called_once_with(file)


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


def test_upload_multipart_failure_reports_failed_state(tmp_path, monkeypatch):
    """Multipart upload failures should be reported as FAILED instead of bubbling silently."""
    source = tmp_path / "broken.bin"
    source.write_bytes(b"abcdefghij")
    file = SaveFileModel(
        source_path=str(source),
        name="logs/broken.bin",
        target_path=str(tmp_path / "run" / "logs" / "broken.bin"),
    )
    manager = save_manager.FileUploadManager(mode="disabled", file_dir=str(tmp_path))
    manager._client = MagicMock()

    captured = {}

    def fake_prepare_multipart(client, exp_id, payload):
        return {
            "uploadId": "upload-1",
            "parts": [
                {"partNumber": 1, "url": "https://upload.example.com/part-1"},
                {"partNumber": 2, "url": "https://upload.example.com/part-2"},
            ],
        }

    def fake_upload_file(*, url, buffer, max_retries=3, mime_type=None):
        raise RuntimeError("Simulated multipart upload failure")

    def fake_complete_upload(client, exp_id, files):
        captured["complete"] = files

    monkeypatch.setattr(save_manager, "PART_SIZE", 5)
    monkeypatch.setattr(save_manager, "prepare_multipart", fake_prepare_multipart)
    monkeypatch.setattr(save_manager, "upload_file", fake_upload_file)
    monkeypatch.setattr(save_manager, "complete_upload", fake_complete_upload)

    try:
        manager._upload_multipart(file, source.stat().st_size, "exp-789")
    finally:
        manager.close()

    assert captured["complete"] == [{"path": "logs/broken.bin", "state": "FAILED"}]


def test_upload_multipart_missing_upload_id_reports_failed_state(tmp_path, monkeypatch):
    """Missing uploadId should report FAILED with an explicit warning message."""
    source = tmp_path / "broken.bin"
    source.write_bytes(b"abcdefghij")
    file = SaveFileModel(
        source_path=str(source),
        name="logs/broken.bin",
        target_path=str(tmp_path / "run" / "logs" / "broken.bin"),
    )
    manager = save_manager.FileUploadManager(mode="disabled", file_dir=str(tmp_path))
    manager._client = MagicMock()

    captured = {}
    warnings = []

    def fake_prepare_multipart(client, exp_id, payload):
        return {
            "parts": [
                {"partNumber": 1, "url": "https://upload.example.com/part-1"},
            ],
        }

    def fake_complete_upload(client, exp_id, files):
        captured["complete"] = files

    monkeypatch.setattr(save_manager, "prepare_multipart", fake_prepare_multipart)
    monkeypatch.setattr(save_manager, "complete_upload", fake_complete_upload)
    monkeypatch.setattr(save_manager.swanlog, "warning", warnings.append)

    try:
        manager._upload_multipart(file, source.stat().st_size, "exp-790")
    finally:
        manager.close()

    assert captured["complete"] == [{"path": "logs/broken.bin", "state": "FAILED"}]
    assert any("missing uploadId" in warning for warning in warnings)


def test_save_progress_starts_without_touching_rich_private_state(monkeypatch):
    """Save progress should rely on Status.start instead of private rich internals."""

    class FakeStatus:
        def __init__(self, *args, **kwargs):
            self.started = 0
            self.stopped = 0
            self.messages = []

        def start(self):
            self.started += 1

        def stop(self):
            self.stopped += 1

        def update(self, message):
            self.messages.append(message)

    monkeypatch.setattr(save_progress, "Status", FakeStatus)

    progress = save_progress._SaveProgress(1)
    progress.add(1)
    progress.done()

    assert progress._status.started >= 1
    assert progress._status.messages[-1] == "Uploading files (1/2)..."


def test_dir_watcher_cloud_falls_back_to_copy_when_symlink_unavailable(tmp_path, monkeypatch):
    source = tmp_path / "metrics.txt"
    source.write_text("hello world")
    target = tmp_path / "run" / "logs" / "metrics.txt"
    file = SaveFileModel(
        source_path=str(source),
        name="logs/metrics.txt",
        target_path=str(target),
    )
    watcher = save_manager.DirWatcher(
        mode="cloud",
        file_dir=str(tmp_path / "run"),
        on_change=MagicMock(),
    )

    def raise_symlink_unavailable(src, dst):
        raise OSError("symlink unavailable")

    monkeypatch.setattr(watcher, "_start", lambda: None)
    monkeypatch.setattr(save_manager.os, "symlink", raise_symlink_unavailable)

    watcher.watch(file)

    assert target.exists()
    assert not target.is_symlink()
    assert target.read_text() == "hello world"
    assert watcher._target_modes[file.name] == "copy"


def test_dir_watcher_cloud_refreshes_copied_target_on_change(tmp_path, monkeypatch):
    source = tmp_path / "metrics.txt"
    source.write_text("hello world")
    target = tmp_path / "run" / "logs" / "metrics.txt"
    file = SaveFileModel(
        source_path=str(source),
        name="logs/metrics.txt",
        target_path=str(target),
    )
    on_change = MagicMock()
    watcher = save_manager.DirWatcher(
        mode="cloud",
        file_dir=str(tmp_path / "run"),
        on_change=on_change,
    )

    def raise_symlink_unavailable(src, dst):
        raise OSError("symlink unavailable")

    monkeypatch.setattr(watcher, "_start", lambda: None)
    monkeypatch.setattr(save_manager.os, "symlink", raise_symlink_unavailable)

    watcher.watch(file)

    source.write_text("hello world updated")
    os.utime(source, None)
    watcher._poll_task()

    assert target.read_text() == "hello world updated"
    on_change.assert_called_once_with(file)


def test_dir_watcher_cloud_falls_back_to_copy_when_relpath_raises_value_error(tmp_path, monkeypatch):
    source = tmp_path / "metrics.txt"
    source.write_text("hello world")
    target = tmp_path / "run" / "logs" / "metrics.txt"
    file = SaveFileModel(
        source_path=str(source),
        name="logs/metrics.txt",
        target_path=str(target),
    )
    watcher = save_manager.DirWatcher(
        mode="cloud",
        file_dir=str(tmp_path / "run"),
        on_change=MagicMock(),
    )
    symlink = MagicMock()

    def raise_cross_drive_relpath(path, start):
        raise ValueError("path is on mount 'D:', start on mount 'C:'")

    monkeypatch.setattr(watcher, "_start", lambda: None)
    monkeypatch.setattr(save_manager.os.path, "relpath", raise_cross_drive_relpath)
    monkeypatch.setattr(save_manager.os, "symlink", symlink)

    watcher.watch(file)

    symlink.assert_not_called()
    assert target.exists()
    assert not target.is_symlink()
    assert target.read_text() == "hello world"
    assert watcher._target_modes[file.name] == "copy"


def test_dir_watcher_cloud_rechecks_cached_symlink_target(tmp_path, monkeypatch):
    source = tmp_path / "metrics.txt"
    source.write_text("hello world")
    target = tmp_path / "run" / "logs" / "metrics.txt"
    file = SaveFileModel(
        source_path=str(source),
        name="logs/metrics.txt",
        target_path=str(target),
    )
    watcher = save_manager.DirWatcher(
        mode="cloud",
        file_dir=str(tmp_path / "run"),
        on_change=MagicMock(),
    )
    copy_file = MagicMock()
    resolve_sync_mode = MagicMock(return_value="copy")

    watcher._target_modes[file.name] = "symlink"
    monkeypatch.setattr(watcher, "_is_symlink_target", lambda current: False)
    monkeypatch.setattr(watcher, "_resolve_sync_mode", resolve_sync_mode)
    monkeypatch.setattr(save_manager, "copy_file", copy_file)

    watcher._sync_live_target(file)

    resolve_sync_mode.assert_called_once_with(file)
    copy_file.assert_called_once_with(file.source_path, file.target_path)
    assert watcher._target_modes.get(file.name) != "symlink"


def test_dir_watcher_guards_target_modes_access_with_lock(tmp_path, monkeypatch):
    """DirWatcher should guard _target_modes reads and writes with the dedicated lock."""

    class GuardedLock:
        def __init__(self):
            self.held = False

        def __enter__(self):
            self.held = True
            return self

        def __exit__(self, exc_type, exc, tb):
            self.held = False

    class GuardedDict(dict):
        def __init__(self, lock):
            super().__init__()
            self._lock = lock

        def _assert_locked(self):
            assert self._lock.held, "target_modes lock not held"

        def get(self, key, default=None):
            self._assert_locked()
            return super().get(key, default)

        def pop(self, key, default=None):
            self._assert_locked()
            return super().pop(key, default)

        def __setitem__(self, key, value):
            self._assert_locked()
            return super().__setitem__(key, value)

    source = tmp_path / "metrics.txt"
    source.write_text("hello world")
    target = tmp_path / "run" / "logs" / "metrics.txt"
    file = SaveFileModel(
        source_path=str(source),
        name="logs/metrics.txt",
        target_path=str(target),
    )
    watcher = save_manager.DirWatcher(
        mode="cloud",
        file_dir=str(tmp_path / "run"),
        on_change=MagicMock(),
    )
    lock = GuardedLock()
    watcher._target_modes_lock = lock
    watcher._target_modes = GuardedDict(lock)

    def raise_symlink_unavailable(src, dst):
        raise OSError("symlink unavailable")

    monkeypatch.setattr(watcher, "_is_symlink_target", lambda current: False)
    monkeypatch.setattr(save_manager.os, "symlink", raise_symlink_unavailable)
    monkeypatch.setattr(save_manager, "copy_file", MagicMock())

    watcher._sync_live_target(file)

    with watcher._target_modes_lock:
        assert watcher._target_modes[file.name] == "copy"
