"""
@author: cunyue
@file: test_pkg_fork.py
@time: 2026/4/19
@description: 测试 swanlab.sdk.internal.pkg.fork 模块
"""

import multiprocessing
import os

import pytest

from swanlab.sdk.internal.pkg import fork


class TestCurrentPid:
    def test_returns_current_pid(self):
        """current_pid() 应返回当前进程 PID"""
        assert fork.current_pid() == os.getpid()


class TestIsForked:
    def test_returns_false_when_pid_matches(self):
        """pre_pid 与当前 PID 一致时返回 False"""
        assert not fork.is_forked(fork.current_pid())

    def test_returns_true_when_pid_differs(self):
        """pre_pid 与当前 PID 不一致时返回 True"""
        assert fork.is_forked(0)

    @pytest.mark.skipif(not hasattr(os, "register_at_fork"), reason="fork not available on this platform")
    def test_returns_true_in_forked_child(self):
        """真实 fork 后，子进程中 is_forked(parent_pid) 应返回 True"""
        parent_pid = fork.current_pid()

        ctx = multiprocessing.get_context("fork")
        result = ctx.Queue()

        def child():
            result.put(fork.is_forked(parent_pid))

        p = ctx.Process(target=child)
        p.start()
        p.join(timeout=5)
        assert result.get() is True

    @pytest.mark.skipif(not hasattr(os, "register_at_fork"), reason="fork not available on this platform")
    def test_returns_false_with_own_pid_in_child(self):
        """子进程中 is_forked(child_pid) 应返回 False"""
        ctx = multiprocessing.get_context("fork")
        result = ctx.Queue()

        def child():
            result.put(fork.is_forked(fork.current_pid()))

        p = ctx.Process(target=child)
        p.start()
        p.join(timeout=5)
        assert result.get() is False


class TestRegister:
    @pytest.mark.skipif(not hasattr(os, "register_at_fork"), reason="fork not available on this platform")
    def test_callback_is_called_after_fork(self):
        """注册的回调应在 fork 后的子进程中被调用"""
        results = []

        fork.register(lambda: results.append("called"))

        ctx = multiprocessing.get_context("fork")
        result = ctx.Queue()

        def child():
            result.put(list(results))

        p = ctx.Process(target=child)
        p.start()
        p.join(timeout=5)
        assert result.get() == ["called"]

    @pytest.mark.skipif(not hasattr(os, "register_at_fork"), reason="fork not available on this platform")
    def test_multiple_callbacks_execute_in_order(self):
        """多个回调按注册顺序执行"""
        order = []

        fork.register(lambda: order.append(1))
        fork.register(lambda: order.append(2))
        fork.register(lambda: order.append(3))

        ctx = multiprocessing.get_context("fork")
        result = ctx.Queue()

        def child():
            result.put(list(order))

        p = ctx.Process(target=child)
        p.start()
        p.join(timeout=5)
        assert result.get() == [1, 2, 3]


class TestUnregister:
    def test_unregister_removes_callback(self):
        """unregister 后回调不再存在于列表中"""
        calls = []

        def cb():
            calls.append(1)

        fork.register(cb)
        fork.unregister(cb)
        # 直接检查 _callbacks 列表
        from swanlab.sdk.internal.pkg.fork import _callbacks

        assert cb not in _callbacks

    def test_unregister_nonexistent_silently_ignores(self):
        """unregister 不存在的回调不抛异常"""
        fork.unregister(lambda: None)  # 不应抛异常

    @pytest.mark.skipif(not hasattr(os, "register_at_fork"), reason="fork not available on this platform")
    def test_unregistered_callback_not_called_after_fork(self):
        """unregister 后的回调在 fork 子进程中不再被调用"""
        calls = []

        def cb():
            calls.append("should_not_be_called")

        fork.register(cb)
        fork.unregister(cb)

        ctx = multiprocessing.get_context("fork")
        result = ctx.Queue()

        def child():
            result.put(list(calls))

        p = ctx.Process(target=child)
        p.start()
        p.join(timeout=5)
        assert result.get() == []
