from pathlib import Path
from unittest.mock import MagicMock

from swanlab.sdk.internal.core_python.watcher import FileWatcher


def _make_watcher() -> FileWatcher:
    return FileWatcher(on_change=MagicMock(), debounce_delay=0.1)


# ── watch idempotency ──


def test_watch_skips_already_registered_file(tmp_path: Path):
    watcher = _make_watcher()
    f = tmp_path / "model.pt"
    f.write_bytes(b"v1")

    watcher.watch(str(tmp_path), ["model.pt"])
    assert len(watcher._registered) == 1

    # 修改文件内容后再 watch 同一文件
    f.write_bytes(b"v2")
    watcher.watch(str(tmp_path), ["model.pt"])

    # 仍然只有一条记录，且签名是第一次注册时的
    assert len(watcher._registered) == 1
    abs_path = str((tmp_path / "model.pt").resolve())
    assert watcher._registered[abs_path].signature is not None


def test_watch_registers_different_files(tmp_path: Path):
    watcher = _make_watcher()
    (tmp_path / "a.pt").write_bytes(b"a")
    (tmp_path / "b.pt").write_bytes(b"b")

    watcher.watch(str(tmp_path), ["a.pt"])
    watcher.watch(str(tmp_path), ["b.pt"])

    assert len(watcher._registered) == 2


def test_watch_batch_idempotent(tmp_path: Path):
    watcher = _make_watcher()
    (tmp_path / "a.pt").write_bytes(b"a")
    (tmp_path / "b.pt").write_bytes(b"b")

    watcher.watch(str(tmp_path), ["a.pt", "b.pt"])
    watcher.watch(str(tmp_path), ["a.pt", "b.pt"])

    assert len(watcher._registered) == 2


# ── observer started once ──


def test_observer_starts_only_once(tmp_path: Path):
    watcher = _make_watcher()
    (tmp_path / "a.pt").write_bytes(b"a")
    (tmp_path / "b.pt").write_bytes(b"b")

    watcher.watch(str(tmp_path), ["a.pt"])
    assert watcher._started is True

    # 第二次 watch 不再调用 observer.schedule
    observer_mock = MagicMock()
    watcher._observer = observer_mock
    watcher.watch(str(tmp_path), ["b.pt"])
    observer_mock.schedule.assert_not_called()
