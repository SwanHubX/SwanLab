from pathlib import Path
from unittest.mock import MagicMock

from swanlab.proto.swanlab.save.v1.save_pb2 import SavePolicy, SaveRecord
from swanlab.sdk.internal.core_python.watcher import FileWatcher


def _make_watcher() -> FileWatcher:
    return FileWatcher(on_change=MagicMock(), debounce_delay=0.1)


def _make_save(name: str, source: Path, policy: SavePolicy = SavePolicy.SAVE_POLICY_LIVE) -> SaveRecord:
    return SaveRecord(name=name, source_path=str(source), policy=policy)


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
    assert watcher._registered[abs_path][0].signature is not None


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


# ── direct-source（skip_store，无本地镜像）──


def test_watch_sources_and_notifies_each_name(tmp_path: Path):
    """同一源文件保存为多个 name：一对多注册，变化时每个 name 各回调一次。"""
    calls = []
    watcher = FileWatcher(on_change=calls.append, debounce_delay=0.1)
    source = tmp_path / "model.pt"
    source.write_bytes(b"v1")

    watcher.watch_sources([_make_save("model.pt", source), _make_save("weights/model.pt", source)])

    key = str(source.resolve())
    assert list(watcher._registered) == [key]
    assert [e.name for e in watcher._registered[key]] == ["model.pt", "weights/model.pt"]
    assert watcher._registered[key][0].source_path == key
    assert watcher._registered[key][0].target_path == ""
    # 重复注册幂等
    watcher.watch_sources([_make_save("model.pt", source)])
    assert len(watcher._registered[key]) == 2

    for entry in watcher._registered[key]:
        entry.signature = "stale"
    source.write_bytes(b"v2")
    watcher._process_change(key)

    assert sorted(r.name for r in calls) == ["model.pt", "weights/model.pt"]
    assert all(r.source_path == key and r.target_path == "" for r in calls)


def test_watch_sources_ignores_other_files_in_same_dir(tmp_path: Path):
    on_change = MagicMock()
    watcher = FileWatcher(on_change=on_change, debounce_delay=0.1)
    source = tmp_path / "model.pt"
    other = tmp_path / "other.pt"
    source.write_bytes(b"v1")
    other.write_bytes(b"v1")

    watcher.watch_sources([_make_save("model.pt", source)])

    other_abs = str(other.resolve())
    watcher._schedule_debounce(other_abs)
    assert watcher._timers == {}
    watcher._process_change(other_abs)
    on_change.assert_not_called()


def test_watch_sources_deleted_source_removes_registration(tmp_path: Path):
    on_change = MagicMock()
    watcher = FileWatcher(on_change=on_change, debounce_delay=0.1)
    source = tmp_path / "model.pt"
    source.write_bytes(b"v1")

    watcher.watch_sources([_make_save("model.pt", source)])
    key = str(source.resolve())
    source.unlink()
    watcher._process_change(key)

    assert key not in watcher._registered
    on_change.assert_not_called()
