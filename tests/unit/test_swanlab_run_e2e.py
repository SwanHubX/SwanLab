"""
@file: test_run_e2e.py
@description: 测试 swanlab.run 属性的端到端行为
全部使用 mode='disabled'，避免文件系统和网络依赖。
"""

import swanlab
from swanlab.sdk.cmd.init import init
from swanlab.sdk.cmd.run import finish
from swanlab.sdk.internal.run import SwanLabRun


class TestRunAttributeE2E:
    def test_run_is_none_before_init(self):
        """未调用 init() 时，swanlab.run 应为 None"""
        assert swanlab.run is None

    def test_run_is_swanlab_run_after_init(self):
        """init() 之后，swanlab.run 应返回活跃的 SwanLabRun 实例"""
        run = init(mode="disabled")
        assert swanlab.run is not None
        assert isinstance(swanlab.run, SwanLabRun)
        assert swanlab.run is run
        finish()

    def test_run_equals_get_run(self):
        """init() 之后，swanlab.run 应与 swanlab.get_run() 返回同一对象"""
        init(mode="disabled")
        assert swanlab.run is swanlab.get_run()
        finish()

    def test_run_is_none_after_finish(self):
        """finish() 之后，swanlab.run 应重新变为 None"""
        run1 = init(mode="disabled")
        assert swanlab.run is not None
        finish()
        assert swanlab.run is None
        run2 = init(mode="disabled")
        assert swanlab.run is not None
        assert swanlab.run is run2
        assert run1 is not run2
