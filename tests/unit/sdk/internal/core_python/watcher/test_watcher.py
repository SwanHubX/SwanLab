import os
import time
from pathlib import Path
from typing import Callable, List
from unittest.mock import MagicMock

from watchdog.events import FileMovedEvent

from swanlab.proto.swanlab.save.v1.save_pb2 import SavePolicy, SaveRecord
from swanlab.sdk.internal.core_python.watcher import FileWatcher
from swanlab.sdk.internal.core_python.watcher.helper import _Handler


def _make_watcher() -> FileWatcher:
    return FileWatcher(on_change=MagicMock(), debounce_delay=0.1)


def _make_save(name: str, source: Path, policy: SavePolicy = SavePolicy.SAVE_POLICY_LIVE) -> SaveRecord:
    return SaveRecord(name=name, source_path=str(source), policy=policy)


def _wait_until(predicate: Callable[[], bool], timeout: float = 15.0) -> bool:
    """轮询等待条件成立，用于真实 observer 测试（事件到达时间因平台而异）。"""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return True
        time.sleep(0.05)
    return False


# Windows 的 os.stat 读 NTFS 目录项缓存，句柄关闭后约 1s 内文件签名仍然是旧值；
# 若 debounce 定时在该窗口内触发，签名比较会读到旧签名而误判"未变化"。
REAL_OBSERVER_DEBOUNCE = 1.5


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


def test_watch_sources_missing_source_keeps_registration(tmp_path: Path):
    """文件暂时缺失时保留注册：不回调，且删除后重建仍能触发回调。"""
    on_change = MagicMock()
    watcher = FileWatcher(on_change=on_change, debounce_delay=0.1)
    source = tmp_path / "model.pt"
    source.write_bytes(b"v1")

    watcher.watch_sources([_make_save("model.pt", source)])
    key = str(source.resolve())

    # debounce 窗口内文件被删除：不回调，也不移除注册
    source.unlink()
    watcher._process_change(key)
    assert key in watcher._registered
    on_change.assert_not_called()

    # 删除后重建：签名变化触发回调
    source.write_bytes(b"v2")
    watcher._process_change(key)
    on_change.assert_called_once()
    assert on_change.call_args[0][0].source_path == key


def test_on_moved_matches_registered_dest_only(tmp_path: Path):
    """on_moved 按 dest_path 匹配注册表：tmp 源路径与未注册目标都不触发。"""
    watcher = _make_watcher()
    source = tmp_path / "model.pt"
    source.write_bytes(b"v1")
    watcher.watch_sources([_make_save("model.pt", source)])
    key = str(source.resolve())

    handler = _Handler(watcher)
    handler.on_moved(FileMovedEvent(src_path=str(tmp_path / "model.pt.tmp"), dest_path=key))
    assert key in watcher._timers

    # 同目录其他文件的原子替换：dest 未注册，不触发
    watcher._timers.clear()
    handler.on_moved(FileMovedEvent(src_path=str(tmp_path / "a.tmp"), dest_path=str(tmp_path / "other.pt")))
    assert watcher._timers == {}


# ── 真实 observer（跨平台事件语义）──


def test_watch_sources_real_observer_atomic_replace(tmp_path: Path):
    """真实 observer 下 tmp + os.replace 的原子替换必须触发回调。

    Linux/Windows 将原子替换上报为 moved 事件而非 modified，依赖 on_moved；
    macOS/FSEvents 可能附带 modified 事件，因此断言只针对最终回调。
    """
    calls: List[SaveRecord] = []
    watcher = FileWatcher(on_change=calls.append, debounce_delay=REAL_OBSERVER_DEBOUNCE)
    source = tmp_path / "model.pt"
    source.write_bytes(b"v1")
    watcher.watch_sources([_make_save("model.pt", source)])

    try:
        # in-place 写作为 warm-up，确认 observer 已开始接收事件
        source.write_bytes(b"v2-warm-up")
        assert _wait_until(lambda: len(calls) >= 1), "observer warm-up write not detected"

        tmp = tmp_path / "model.pt.tmp"
        tmp.write_bytes(b"v3-atomic-replace")
        os.replace(tmp, source)
        assert _wait_until(lambda: len(calls) >= 2), f"atomic replace not detected, calls: {len(calls)}"
        assert calls[-1].source_path == str(source.resolve())
    finally:
        watcher.stop()


def test_watch_sources_real_observer_delete_and_recreate(tmp_path: Path):
    """删除后重建的文件仍在监听范围内：debounce 窗口内删除（定时器在文件缺失时触发）不得移除注册。"""
    calls: List[SaveRecord] = []
    watcher = FileWatcher(on_change=calls.append, debounce_delay=REAL_OBSERVER_DEBOUNCE)
    source = tmp_path / "model.pt"
    source.write_bytes(b"v1")
    watcher.watch_sources([_make_save("model.pt", source)])

    try:
        source.write_bytes(b"v2-warm-up")
        assert _wait_until(lambda: len(calls) >= 1), "observer warm-up write not detected"

        # 写入后立刻删除：pending 的 debounce 定时器将在文件缺失时触发，
        # 注册不得因此被移除，重建后仍要能触发回调
        source.write_bytes(b"v3-doomed")
        source.unlink()
        time.sleep(REAL_OBSERVER_DEBOUNCE + 0.3)
        source.write_bytes(b"v4-recreated")
        assert _wait_until(lambda: len(calls) >= 2), f"recreated file not detected, calls: {len(calls)}"
    finally:
        watcher.stop()
