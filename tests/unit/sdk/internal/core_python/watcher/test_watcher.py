import threading
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
from watchdog.events import FileMovedEvent

from swanlab.proto.swanlab.save.v1.save_pb2 import SavePolicy, SaveRecord
from swanlab.sdk.internal.core_python import watcher as watcher_module
from swanlab.sdk.internal.core_python.watcher import FileWatcher, NullFileWatcher
from swanlab.sdk.internal.core_python.watcher import helper as watcher_helper
from swanlab.sdk.internal.core_python.watcher.helper import _Handler


@pytest.fixture
def mock_watcher_threads(monkeypatch):
    monkeypatch.setattr(watcher_module, "Observer", MagicMock())
    monkeypatch.setattr(watcher_module.threading, "Timer", MagicMock())


def _make_watcher() -> FileWatcher:
    return FileWatcher(on_change=MagicMock(), debounce_delay=0.1)


# ── watch idempotency ──


@pytest.mark.usefixtures("mock_watcher_threads")
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
    assert watcher._registered[abs_path][0].signature is not None


@pytest.mark.usefixtures("mock_watcher_threads")
def test_watch_registers_different_files(tmp_path: Path):
    watcher = _make_watcher()
    (tmp_path / "a.pt").write_bytes(b"a")
    (tmp_path / "b.pt").write_bytes(b"b")

    watcher.watch(str(tmp_path), ["a.pt"])
    watcher.watch(str(tmp_path), ["b.pt"])

    assert len(watcher._registered) == 2


@pytest.mark.usefixtures("mock_watcher_threads")
def test_watch_batch_idempotent(tmp_path: Path):
    watcher = _make_watcher()
    (tmp_path / "a.pt").write_bytes(b"a")
    (tmp_path / "b.pt").write_bytes(b"b")

    watcher.watch(str(tmp_path), ["a.pt", "b.pt"])
    watcher.watch(str(tmp_path), ["a.pt", "b.pt"])

    assert len(watcher._registered) == 2


# ── observer started once ──


@pytest.mark.usefixtures("mock_watcher_threads")
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


# ── 事件匹配与签名 ──


@pytest.mark.usefixtures("mock_watcher_threads")
def test_on_moved_matches_registered_dest_only(tmp_path: Path):
    """on_moved 按 dest_path 匹配注册表：tmp 源路径与未注册目标都不触发。"""
    watcher = _make_watcher()
    source = tmp_path / "model.pt"
    source.write_bytes(b"v1")
    watcher.watch(str(tmp_path), ["model.pt"])
    key = str((tmp_path / "model.pt").resolve())

    handler = _Handler(watcher)
    handler.on_moved(FileMovedEvent(src_path=str(tmp_path / "model.pt.tmp"), dest_path=key))
    assert key in watcher._timers

    # 同目录其他文件的原子替换：dest 未注册，不触发
    watcher._timers.pop(key).cancel()
    handler.on_moved(FileMovedEvent(src_path=str(tmp_path / "a.tmp"), dest_path=str(tmp_path / "other.pt")))
    assert watcher._timers == {}


@pytest.mark.usefixtures("mock_watcher_threads")
def test_stat_error_keeps_signature(tmp_path: Path, monkeypatch):
    on_change = MagicMock()
    watcher = FileWatcher(on_change=on_change, debounce_delay=0.1)
    source = tmp_path / "model.pt"
    source.write_bytes(b"v1")
    watcher.watch(str(tmp_path), ["model.pt"])
    key = str((tmp_path / "model.pt").resolve())
    original_signature = watcher._registered[key][0].signature
    stat = MagicMock(side_effect=PermissionError("access denied"))
    monkeypatch.setattr(watcher_helper, "os", SimpleNamespace(stat=stat))

    watcher._process_change(key)

    assert watcher._registered[key][0].signature == original_signature
    on_change.assert_not_called()


def test_signature_detects_inode_change(monkeypatch):
    stat = MagicMock(
        side_effect=[
            SimpleNamespace(st_mtime_ns=100, st_size=2, st_ino=7),
            SimpleNamespace(st_mtime_ns=100, st_size=2, st_ino=8),
        ]
    )
    monkeypatch.setattr(watcher_helper, "os", SimpleNamespace(stat=stat))

    assert watcher_helper.compute_signature("model.pt") != watcher_helper.compute_signature("model.pt")


# ── NullFileWatcher（skip_store）──


@pytest.mark.usefixtures("mock_watcher_threads")
def test_null_watcher_is_noop(tmp_path: Path):
    """空监听器：注册 LIVE 记录与 stop 均为无操作，不创建监听线程。"""
    watcher = NullFileWatcher()
    save = SaveRecord(name="model.pt", source_path=str(tmp_path / "model.pt"), policy=SavePolicy.SAVE_POLICY_LIVE)

    threads_before = threading.active_count()
    watcher.register_live_watches([save], tmp_path)
    watcher.stop()

    assert threading.active_count() == threads_before
