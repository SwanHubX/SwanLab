"""
@author: Nexisato
@file: test_lock.py
@time: 2026/11/16
@description: per-user runtime 文件启动锁（core.lock）测试。
"""

import os
import sys
import threading
import time
from pathlib import Path

import pytest

from swanlab.sdk.internal.core_go.config import LOCK_FILE_NAME, RUNTIME_DIR_MODE
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.process import acquire_start_lock, try_acquire_start_lock
from swanlab.sdk.internal.core_go.process.lock import _create_lock_file, _parse_lock


def _make_runtime_dir(tmp_path: Path) -> Path:
    runtime_dir = tmp_path / "swanlab-core-test"
    runtime_dir.mkdir(mode=RUNTIME_DIR_MODE)
    if sys.platform != "win32":
        os.chmod(runtime_dir, RUNTIME_DIR_MODE)
    return runtime_dir


def test_try_acquire_and_release_lock(tmp_path):
    runtime_dir = _make_runtime_dir(tmp_path)
    now_ns = time.time_ns()

    lock = try_acquire_start_lock(runtime_dir, now_ns_fn=lambda: now_ns)
    assert lock is not None
    assert lock.path == runtime_dir / LOCK_FILE_NAME
    assert lock.pid == os.getpid()
    assert lock.path.is_file()

    lock.release()
    lock.release()
    assert not lock.path.exists()


def test_release_does_not_delete_replaced_lock(tmp_path):
    runtime_dir = _make_runtime_dir(tmp_path)
    base_ns = time.time_ns()

    old_lock = try_acquire_start_lock(runtime_dir, now_ns_fn=lambda: base_ns)
    assert old_lock is not None

    old_lock.path.unlink()
    assert _create_lock_file(old_lock.path, os.getpid() + 1, base_ns + 1_000_000_000)

    old_lock.release()
    assert old_lock.path.exists()


def test_try_acquire_respects_live_holder_and_reclaims_stale(tmp_path):
    runtime_dir = _make_runtime_dir(tmp_path)
    now_ns = time.time_ns()
    lock_path = runtime_dir / LOCK_FILE_NAME

    assert _create_lock_file(lock_path, 424242, now_ns)

    # 持有者存活且未超时：不能抢占。
    assert (
        try_acquire_start_lock(
            runtime_dir,
            stale_seconds=30.0,
            now_ns_fn=lambda: now_ns,
            is_alive_fn=lambda _: True,
        )
        is None
    )

    # 持有者已死：立即回收并抢锁成功。
    reclaimed = try_acquire_start_lock(
        runtime_dir,
        stale_seconds=30.0,
        now_ns_fn=lambda: now_ns,
        is_alive_fn=lambda _: False,
    )
    assert reclaimed is not None
    reclaimed.release()

    # 持有者报告存活但超过 stale_seconds：超时回收。
    assert _create_lock_file(lock_path, 424242, now_ns)
    later_ns = now_ns + 31_000_000_000
    expired_lock = try_acquire_start_lock(
        runtime_dir,
        stale_seconds=30.0,
        now_ns_fn=lambda: later_ns,
        is_alive_fn=lambda _: True,
    )
    assert expired_lock is not None
    expired_lock.release()


def test_try_acquire_concurrent_single_winner(tmp_path):
    runtime_dir = _make_runtime_dir(tmp_path)
    barrier = threading.Barrier(8)
    winners = []
    errors = []
    guard = threading.Lock()

    def worker():
        try:
            barrier.wait()
            lock = try_acquire_start_lock(
                runtime_dir,
                is_alive_fn=lambda _: True,
            )
            if lock is not None:
                with guard:
                    winners.append(lock)
        except BaseException as exc:  # pragma: no cover
            with guard:
                errors.append(exc)

    threads = [threading.Thread(target=worker) for _ in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert errors == []
    assert len(winners) == 1
    winners[0].release()


def test_acquire_start_lock_waits_and_times_out(tmp_path):
    runtime_dir = _make_runtime_dir(tmp_path)
    lock_path = runtime_dir / LOCK_FILE_NAME
    now = [time.time_ns()]
    assert _create_lock_file(lock_path, 999999, now[0])

    # 第一次等待后前一个持有者释放锁，第二次尝试应立即获取。
    def release_on_sleep(seconds: float) -> None:
        now[0] += int(seconds * 1_000_000_000)
        lock_path.unlink()

    with acquire_start_lock(
        runtime_dir,
        timeout=1.0,
        stale_seconds=60.0,
        poll_interval=0.01,
        now_ns_fn=lambda: now[0],
        is_alive_fn=lambda _: True,
        sleep_fn=release_on_sleep,
    ) as lock:
        assert lock.path == lock_path
    assert not lock_path.exists()

    # 持续被占用时按注入时钟超时报错。
    assert _create_lock_file(lock_path, 999999, now[0])

    def advance_clock(seconds: float) -> None:
        now[0] += int(seconds * 1_000_000_000)

    with pytest.raises(CoreGoError) as caught:
        acquire_start_lock(
            runtime_dir,
            timeout=0.05,
            stale_seconds=60.0,
            poll_interval=0.02,
            now_ns_fn=lambda: now[0],
            is_alive_fn=lambda _: True,
            sleep_fn=advance_clock,
        )
    assert caught.value.stage is Stage.RUNTIME
    assert caught.value.reason is ReasonCode.RUNTIME_UNAVAILABLE


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX symlink check")
def test_try_acquire_rejects_symlink_lock(tmp_path):
    runtime_dir = _make_runtime_dir(tmp_path)
    outside = tmp_path / "outside"
    outside.write_text("x")
    (runtime_dir / LOCK_FILE_NAME).symlink_to(outside)

    with pytest.raises(CoreGoError) as caught:
        try_acquire_start_lock(runtime_dir)

    assert caught.value.stage is Stage.RUNTIME
    assert caught.value.reason is ReasonCode.RUNTIME_UNAVAILABLE


@pytest.mark.parametrize(
    "content",
    [
        "pid=12\ncreated_at=100\n",
        "created_at=100\nEOF\n",
        "pid=12\nEOF\n",
        "pid=0\ncreated_at=100\nEOF\n",
        "pid=12\ncreated_at=0\nEOF\n",
        "pid=12\npid=12\ncreated_at=100\nEOF\n",
        "pid=12\ncreated_at=100\nextra=1\nEOF\n",
    ],
)
def test_parse_lock_rejects_malformed(content):
    with pytest.raises(ValueError):
        _parse_lock(content.encode("utf-8"))


def test_create_lock_file_uses_binary_flag_and_writes_lf(tmp_path, monkeypatch):
    from swanlab.sdk.internal.core_go.process import lock as lock_module

    runtime_dir = _make_runtime_dir(tmp_path)
    lock_path = runtime_dir / LOCK_FILE_NAME
    seen_flags = []
    real_open = os.open

    has_real_binary = hasattr(os, "O_BINARY")
    binary_flag = getattr(os, "O_BINARY", 0x8000)
    if not has_real_binary:
        monkeypatch.setattr(lock_module.os, "O_BINARY", binary_flag, raising=False)

    def spy_open(path, flags, mode=0o777):
        seen_flags.append(flags)
        real_flags = flags if has_real_binary else (flags & ~binary_flag)
        return real_open(path, real_flags, mode)

    monkeypatch.setattr(lock_module.os, "open", spy_open)

    assert _create_lock_file(lock_path, 123, 456)
    assert seen_flags and all(flag & binary_flag for flag in seen_flags)
    raw = lock_path.read_bytes()
    assert b"\r" not in raw
    assert _parse_lock(raw) == (123, 456)
