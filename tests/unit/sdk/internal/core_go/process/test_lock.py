"""
@author: Nexisato
@file: test_lock.py
@time: 2026/11/16
@description: per-user runtime 内核独占启动锁（core.lock）测试。
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

    # 持锁期间二次抢锁必须失败。
    assert try_acquire_start_lock(runtime_dir) is None

    lock.release()
    lock.release()

    # 释放后新持有者可重新加锁，且旧锁重复 release 不破坏新锁。
    next_lock = try_acquire_start_lock(runtime_dir)
    assert next_lock is not None
    try:
        lock.release()
        assert try_acquire_start_lock(runtime_dir) is None
    finally:
        next_lock.release()


def test_lock_interleaved_no_double_acquire(tmp_path):
    runtime_dir = _make_runtime_dir(tmp_path)
    first = try_acquire_start_lock(runtime_dir)
    assert first is not None

    # 交错场景：A 未真正释放前，B 多次抢锁都必须失败；仅 A 释放后 B 才成功。
    assert try_acquire_start_lock(runtime_dir) is None

    first.release()
    second = try_acquire_start_lock(runtime_dir)
    assert second is not None
    try:
        first.release()
        assert try_acquire_start_lock(runtime_dir) is None
    finally:
        second.release()


def test_try_acquire_concurrent_single_winner(tmp_path):
    runtime_dir = _make_runtime_dir(tmp_path)
    barrier = threading.Barrier(8)
    winners = []
    errors = []
    guard = threading.Lock()

    def worker():
        try:
            barrier.wait()
            lock = try_acquire_start_lock(runtime_dir)
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

    try:
        assert errors == []
        assert len(winners) == 1
    finally:
        for w in winners:
            w.release()


def test_acquire_start_lock_waits_and_times_out(tmp_path):
    runtime_dir = _make_runtime_dir(tmp_path)
    holder = try_acquire_start_lock(runtime_dir)
    assert holder is not None
    now = [time.time_ns()]

    # 第一次等待时前一个持有者释放锁，第二次尝试应立即获取。
    def release_on_sleep(seconds: float) -> None:
        now[0] += int(seconds * 1_000_000_000)
        holder.release()

    with acquire_start_lock(
        runtime_dir,
        timeout=1.0,
        poll_interval=0.01,
        now_ns_fn=lambda: now[0],
        sleep_fn=release_on_sleep,
    ) as lock:
        assert lock.path == runtime_dir / LOCK_FILE_NAME
        assert try_acquire_start_lock(runtime_dir) is None

    # 持续被占用时按注入时钟超时报错。
    holder2 = try_acquire_start_lock(runtime_dir)
    assert holder2 is not None
    try:

        def advance_clock(seconds: float) -> None:
            now[0] += int(seconds * 1_000_000_000)

        with pytest.raises(CoreGoError) as caught:
            acquire_start_lock(
                runtime_dir,
                timeout=0.05,
                poll_interval=0.02,
                now_ns_fn=lambda: now[0],
                sleep_fn=advance_clock,
            )
        assert caught.value.stage is Stage.RUNTIME
        assert caught.value.reason is ReasonCode.RUNTIME_UNAVAILABLE
    finally:
        holder2.release()


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX symlink check")
def test_try_acquire_rejects_symlink_lock(tmp_path):
    runtime_dir = _make_runtime_dir(tmp_path)
    outside = tmp_path / "outside"
    outside.write_bytes(b"x")
    (runtime_dir / LOCK_FILE_NAME).symlink_to(outside)

    with pytest.raises(CoreGoError) as caught:
        try_acquire_start_lock(runtime_dir)

    assert caught.value.stage is Stage.RUNTIME
    assert caught.value.reason is ReasonCode.RUNTIME_UNAVAILABLE


def test_lock_metadata_written_as_lf(tmp_path):
    runtime_dir = _make_runtime_dir(tmp_path)
    lock = try_acquire_start_lock(runtime_dir)
    assert lock is not None
    # Windows LockFileEx 是强制字节范围锁，需先释放句柄再由外部句柄读取校验元数据。
    lock.release()
    raw = lock.path.read_bytes()
    assert b"\r" not in raw
    assert raw.endswith(b"\nEOF\n")
    assert f"pid={os.getpid()}\n".encode("ascii") in raw
