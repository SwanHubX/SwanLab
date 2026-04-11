"""
@author: cunyue
@file: test_async_log_e2e.py
@time: 2026/4/11
@description: swanlab.async_log() + swanlab.finish() 端到端测试

全部使用 mode='disabled'，避免文件系统和网络依赖。
"""

import asyncio
import sys
import time
from unittest.mock import MagicMock

import pytest

import swanlab
from swanlab.sdk.cmd.init import init


# spawn 模式要求 func 可 pickle，局部函数不可 pickle，使用模块级顶层函数
def _compute_score():
    return {"score": 42}


def _slow(n):
    time.sleep(0.05 * n)
    return {"v": n}


class TestAsyncLogE2E:
    """端到端：通过 swanlab 顶层 API 调用 async_log，验证与 Run 生命周期的集成。"""

    def test_threading_mode_auto_logs(self):
        """async_log(threading) 返回值自动传入 _log_impl()，finish() 等待完成"""
        run = init(mode="disabled")
        run._log_impl = MagicMock(wraps=run._log_impl)

        def compute():
            return {"loss": 0.5}

        f = swanlab.async_log(compute, step=1, mode="threading")
        swanlab.finish()

        assert f.done()
        assert f.result() == {"loss": 0.5}
        run._log_impl.assert_called_once_with({"loss": 0.5}, step=1)

    def test_asyncio_mode_auto_logs(self):
        """async_log(asyncio) 协程返回值自动传入 _log_impl()"""
        run = init(mode="disabled")
        run._log_impl = MagicMock(wraps=run._log_impl)

        async def compute():
            await asyncio.sleep(0.01)
            return {"acc": 0.95}

        async def main():
            f = swanlab.async_log(compute, step=2, mode="asyncio")
            await asyncio.wrap_future(f)
            return f

        f = asyncio.run(main())
        swanlab.finish()

        assert f.done()
        assert f.result() == {"acc": 0.95}
        run._log_impl.assert_called_once_with({"acc": 0.95}, step=2)

    @pytest.mark.skipif(sys.platform == "win32", reason="fork is not available on Windows")
    def test_spawn_mode_auto_logs(self):
        """async_log(spawn) 进程池返回值自动传入 _log_impl()"""
        run = init(mode="disabled")
        run._log_impl = MagicMock(wraps=run._log_impl)

        f = swanlab.async_log(_compute_score, step=3, mode="spawn")
        swanlab.finish()

        assert f.done()
        run._log_impl.assert_called_once_with({"score": 42}, step=3)

    def test_error_does_not_crash_finish(self):
        """async_log 任务异常不会阻塞 finish()"""
        run = init(mode="disabled")
        run._log_impl = MagicMock(wraps=run._log_impl)

        def boom():
            raise ValueError("kaboom")

        swanlab.async_log(boom, step=1, mode="threading")
        swanlab.finish()

        # 异常任务不应触发 _log_impl
        run._log_impl.assert_not_called()

    def test_multiple_tasks_all_complete(self):
        """多个 async_log 任务在 finish() 时全部完成"""
        run = init(mode="disabled")
        run._log_impl = MagicMock(wraps=run._log_impl)

        swanlab.async_log(_slow, 2, step=1, mode="threading")
        swanlab.async_log(_slow, 1, step=2, mode="threading")

        swanlab.finish()

        # 两个任务都应完成并调用 _log_impl
        assert run._log_impl.call_count == 2

    def test_asyncio_no_loop_raises(self):
        """无事件循环时 async_log(asyncio) 抛出 RuntimeError"""
        init(mode="disabled")

        with pytest.raises(RuntimeError):
            swanlab.async_log(lambda: {"x": 1}, mode="asyncio")

        swanlab.finish()

    def test_args_kwargs_forwarded(self):
        """*args 和 **kwargs 正确透传到 func"""
        run = init(mode="disabled")
        run._log_impl = MagicMock(wraps=run._log_impl)

        def add(a, b=0):
            return {"sum": a + b}

        f = swanlab.async_log(add, 1, b=9, step=1, mode="threading")
        swanlab.finish()

        assert f.result() == {"sum": 10}
        run._log_impl.assert_called_once_with({"sum": 10}, step=1)
