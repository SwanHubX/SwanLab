"""
@author: cunyue
@file: test_finish_e2e.py
@time: 2026/3/14
@description: 测试 swanlab.finish() 的端到端行为（配合真实 init() 启动 Run）
全部使用 mode='disabled'，避免文件系统和网络依赖。
"""

import swanlab
from swanlab.sdk.cmd.finish import finish
from swanlab.sdk.cmd.init import init
from swanlab.sdk.internal.run import has_run


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

    def test_finish_without_init_is_noop(self):
        """未调用 init() 直接调用 finish() 不应抛出异常"""
        assert not has_run()

        finish()  # 应打印 error 日志，但不 raise

        assert not has_run()

    def test_double_finish_is_noop(self):
        """连续两次 finish()：第二次应为 no-op，不抛出异常"""
        init(mode="disabled")

        finish()
        assert not has_run()

        finish()  # 第二次：无活跃 Run，打印 error 日志，不 raise
        assert not has_run()
