import threading
from unittest.mock import patch

from swanlab.sdk.internal.core_python.pkg.executor import SafeThreadPoolExecutor, _is_shutting_down


def test_submit_returns_future():
    """正常时 submit() 返回 Future，可获取结果。"""
    with SafeThreadPoolExecutor(max_workers=2) as e:
        future = e.submit(lambda x: x * 2, 21)
    assert future.result() == 42


def test_submit_sync_on_shutting_down():
    """shutting down 时 submit() 同步执行并返回已完成 Future。"""
    with patch("swanlab.sdk.internal.core_python.pkg.executor._is_shutting_down", return_value=True):
        e = SafeThreadPoolExecutor(max_workers=2)
        future = e.submit(lambda: 42)
        assert future.result() == 42
        e.shutdown()


def test_submit_captures_exception_on_shutting_down():
    """shutting down 时 submit() 同步执行，异常写入 Future。"""
    with patch("swanlab.sdk.internal.core_python.pkg.executor._is_shutting_down", return_value=True):
        e = SafeThreadPoolExecutor(max_workers=2)
        future = e.submit(lambda: 1 / 0)
        assert future.exception() is not None
        e.shutdown()


def test_run_async_on_normal():
    """正常时 run() 异步提交，shutdown 默认 wait=True 保证完成。"""
    result = []
    with SafeThreadPoolExecutor(max_workers=2) as e:
        e.run(lambda: result.append(1))
    # context manager __exit__ 调 shutdown(wait=True)，任务一定完成
    assert result == [1]


def test_run_sync_on_shutting_down():
    """shutting down 时 run() 同步执行。"""
    result = []
    with patch("swanlab.sdk.internal.core_python.pkg.executor._is_shutting_down", return_value=True):
        e = SafeThreadPoolExecutor(max_workers=2)
        e.run(lambda: result.append(1))
        e.shutdown()
    assert result == [1]


def test_shutdown_default_wait():
    """shutdown() 默认 wait=True，等待已提交任务完成。"""
    result = []
    e = SafeThreadPoolExecutor(max_workers=2)
    e.run(lambda: result.append(1))
    e.shutdown()  # wait=True by default
    assert result == [1]


def test_is_shutting_down():
    """_is_shutting_down() 读取 threading._SHUTTING_DOWN。"""
    assert _is_shutting_down() is False
    with patch.object(threading, "_SHUTTING_DOWN", True):
        assert _is_shutting_down() is True
