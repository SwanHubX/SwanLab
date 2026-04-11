"""
@author: cunyue
@file: test_async_task.py
@time: 2026/4/11
@description: AsyncTaskManager 单元测试
"""

import asyncio
import sys
import time

import pytest

from swanlab.sdk.internal.run.async_task import AsyncTaskManager


# spawn/fork 模式要求 func 可 pickle，lambda 和局部闭包不可 pickle，
# 因此进程模式测试使用模块级顶层函数。
def _double(x):
    return x * 2


def _ret42():
    return 42


@pytest.fixture
def mgr():
    m = AsyncTaskManager()
    yield m
    # fixture 结束时如果还没 shutdown，补一个（防止泄漏）
    # noinspection PyBroadException
    try:
        m.shutdown()
    except Exception:
        pass


# ----------------------------------
# 懒初始化
# ----------------------------------


class TestLazyInit:
    def test_no_pool_on_init(self, mgr):
        assert mgr._thread_pool is None
        assert mgr._spawn_pool is None
        assert mgr._fork_pool is None

    def test_thread_pool_created_on_threading(self, mgr):
        mgr.submit(lambda: 1, mode="threading")
        assert mgr._thread_pool is not None
        assert mgr._spawn_pool is None

    def test_spawn_pool_created_on_spawn(self, mgr):
        mgr.submit(_ret42, mode="spawn")
        assert mgr._spawn_pool is not None
        assert mgr._thread_pool is None

    @pytest.mark.skipif(sys.platform == "win32", reason="fork is not available on Windows")
    def test_fork_pool_created_on_fork(self, mgr):
        mgr.submit(_ret42, mode="fork")
        assert mgr._fork_pool is not None
        assert mgr._thread_pool is None
        assert mgr._spawn_pool is None


# ----------------------------------
# submit + 回调
# ----------------------------------


class TestSubmit:
    def test_threading_success(self, mgr):
        results = []
        f = mgr.submit(lambda: {"loss": 0.5}, step=1, mode="threading", on_success=lambda r, s: results.append((r, s)))
        assert f.result() == {"loss": 0.5}
        assert results == [({"loss": 0.5}, 1)]

    def test_threading_error(self, mgr):
        errors = []

        def boom():
            raise ValueError("kaboom")

        f = mgr.submit(boom, mode="threading", on_error=lambda tb: errors.append(tb))
        with pytest.raises(ValueError, match="kaboom"):
            f.result()
        assert len(errors) == 1
        assert "kaboom" in errors[0]

    def test_spawn_success(self, mgr):
        # func 必须是模块级顶层函数（可 pickle）
        results = []
        f = mgr.submit(_double, args=(3,), mode="spawn", on_success=lambda r, s: results.append(r))
        assert f.result() == 6
        assert results == [6]

    @pytest.mark.skipif(sys.platform == "win32", reason="fork is not available on Windows")
    def test_fork_success(self, mgr):
        # func 必须是模块级顶层函数（可 pickle）
        results = []
        f = mgr.submit(_ret42, mode="fork", on_success=lambda r, s: results.append(r))
        assert f.result() == 42
        assert results == [42]

    def test_asyncio_coroutine(self):
        # asyncio 模式需要运行中的事件循环，用 asyncio.run 提供
        mgr = AsyncTaskManager()
        results = []

        async def coro():
            await asyncio.sleep(0.01)
            return {"v": 9}

        async def main():
            f = mgr.submit(coro, step=2, mode="asyncio", on_success=lambda r, s: results.append((r, s)))
            await asyncio.wrap_future(f)
            assert results == [({"v": 9}, 2)]

        asyncio.run(main())
        mgr.shutdown()

    def test_asyncio_no_loop_raises(self, mgr):
        # 无事件循环时提交 asyncio 模式应抛 RuntimeError
        with pytest.raises(RuntimeError):
            mgr.submit(lambda: 1, mode="asyncio")

    def test_invalid_mode(self, mgr):
        with pytest.raises(ValueError, match="Invalid mode"):
            mgr.submit(lambda: 1, mode="invalid")  # type: ignore


# ----------------------------------
# shutdown
# ----------------------------------


class TestShutdown:
    def test_waits_for_all_and_closes(self, mgr):
        order = []

        def slow(n):
            time.sleep(0.05 * n)
            order.append(n)
            return n

        mgr.submit(slow, args=(3,), mode="threading")
        mgr.submit(slow, args=(1,), mode="threading")
        mgr.shutdown()
        assert sorted(order) == [1, 3]
        # shutdown 后池已关闭，不可再 submit
        with pytest.raises(RuntimeError):
            mgr.submit(lambda: 1, mode="threading")

    def test_error_does_not_block_others(self, mgr):
        # 某个任务异常不应中断对其他任务的等待
        results = []

        def ok():
            return "done"

        def fail():
            raise RuntimeError("oops")

        mgr.submit(ok, mode="threading", on_success=lambda r, s: results.append(r))
        mgr.submit(fail, mode="threading", on_error=lambda tb: results.append("err"))
        mgr.shutdown()
        assert "done" in results
        assert "err" in results

    def test_handles_self_cleanup(self, mgr):
        # 回调完成后 handle 自动从 _handles 中移除，不会无限增长
        for i in range(10):
            mgr.submit(lambda n=i: n, mode="threading")
        # 等待所有回调完成
        for h in list(mgr._handles):
            h.callback_done.wait(timeout=5)
        assert len(mgr._handles) == 0


# ----------------------------------
# args / kwargs / step
# ----------------------------------


class TestArgsKwargs:
    def test_args_and_kwargs(self, mgr):
        results = []

        def add(a, b=0):
            return a + b

        f = mgr.submit(add, args=(1,), kwargs={"b": 2}, mode="threading", on_success=lambda r, s: results.append(r))
        assert f.result() == 3
        assert results == [3]

    def test_step_none(self, mgr):
        results = []
        f = mgr.submit(lambda: 1, step=None, mode="threading", on_success=lambda r, s: results.append(s))
        f.result()
        assert results == [None]
