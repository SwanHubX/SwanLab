"""
@author: cunyue
@file: test_utils.py
@time: 2024/12/5 13:27
@description: 硬件信息采集工具测试
"""

from swanlab.data.run.metadata.hardware.type import HardwareCollector, CollectGuard
from swanlab.data.run.metadata.hardware.type import HardwareInfo


def test_hardware():
    class TestCollector(HardwareCollector):
        def collect(self) -> HardwareInfo:
            return {"key": "test", "value": 1, "name": "test", "config": None}

    t = TestCollector()
    assert t.collect() == {"key": "test", "value": 1, "name": "test", "config": None}

    class TestErrorCollector(HardwareCollector):
        def collect(self) -> HardwareInfo:
            raise Exception("test")

    t = TestErrorCollector()
    assert t() is None


def test_collect_guard():
    class TestGuard(CollectGuard):
        def __init__(self):
            super().__init__()
            self.before_count = 0
            self.after_count = 0

        def before_collect_impl(self):
            self.before_count += 1
            return "before"

        def after_collect_impl(self):
            self.after_count += 1
            return "after"

    t = TestGuard()
    t.before_collect()
    assert t.before_count == 1
    assert t.after_count == 0
    assert t.collect_num == 1
    t.after_collect()
    assert t.before_count == 1
    assert t.after_count == 1
    assert t.collect_num == 1
    # 60>count>0, 不执行
    t.collect_num = 30
    t.before_collect()
    t.after_collect()
    assert t.before_count == 1
    assert t.after_count == 1
    assert t.collect_num == 31
    # 60>count>0, 不执行
    t.collect_num = 59
    t.before_collect()
    t.after_collect()
    assert t.before_count == 1
    assert t.after_count == 1
    assert t.collect_num == 60
    t.before_collect()
    t.after_collect()
    assert t.before_count == 2
    assert t.after_count == 2
    assert t.collect_num == 61
