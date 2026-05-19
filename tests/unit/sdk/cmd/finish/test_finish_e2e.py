"""
@author: cunyue
@file: test_finish_e2e.py
@time: 2026/3/14
@description: 测试 swanlab.finish() 的端到端行为（配合真实 init() 启动 Run）
全部使用 mode='disabled'，避免文件系统和网络依赖。
"""

import pytest

import swanlab
from swanlab.sdk.cmd.init import init
from swanlab.sdk.cmd.run import finish
from swanlab.sdk.internal.run import clear_run, has_run


class TestFinishE2E:
    def test_finish_success_after_init(self):
        """init() → finish() 完整生命周期：Run 创建后正常结束，has_run() 变为 False"""
        run = init(mode="disabled")
        assert run is not None
        assert has_run()

        swanlab.finish()

        assert not has_run()

    def test_finish_crashed_state(self):
        """finish(state='crashed', error=...) 应正常执行，不抛出异常"""
        init(mode="disabled")
        assert has_run()

        finish(state="crashed", error="something went wrong in training")

        assert not has_run()

    def test_finish_without_init_warns(self, monkeypatch):
        """未调用 init() 直接调用 finish() 应告警并返回"""
        clear_run()
        assert not has_run()

        warnings = []
        monkeypatch.setattr("swanlab.sdk.cmd.guard.console.warning", lambda msg, *a, **kw: warnings.append(msg))

        assert finish() is None

        assert not has_run()
        assert warnings == ["SwanLab Run has already finished or has not started."]

    def test_double_finish_warns(self, monkeypatch):
        """连续两次顶层 finish()：第二次只告警，不抛出异常，同时保持全局 Run 已清理"""
        warnings = []
        monkeypatch.setattr(
            "swanlab.sdk.internal.run.console.warning", lambda msg, *args, **kwargs: warnings.append(msg)
        )

        init(mode="disabled")

        finish()
        assert not has_run()
        assert swanlab.run is None

        finish()
        assert not has_run()
        assert swanlab.run is None
        assert warnings == ["SwanLab Run has already finished or has not started."]

    def test_run_finish_warns_when_already_finished(self, monkeypatch):
        """直接持有 Run 对象时，重复 finish() 也只告警，不抛出异常"""
        warnings = []
        monkeypatch.setattr(
            "swanlab.sdk.internal.run.console.warning", lambda msg, *args, **kwargs: warnings.append(msg)
        )

        run = init(mode="disabled")

        run.finish()
        assert not has_run()

        run.finish()
        assert warnings == ["SwanLab Run has already finished or has not started."]

    def test_other_cmd_raises_after_finish(self):
        """finish 后调用 log（非 finish 命令）应抛出 RuntimeError，证明 allow_finished 仅对 finish 生效"""
        init(mode="disabled")
        finish()

        with pytest.raises(RuntimeError, match="`swanlab.log` requires an active Run"):
            swanlab.log({"x": 1})
