"""
@author: cunyue
@file: test_utils.py
@time: 2024/12/5 13:27
@description: 硬件信息采集工具测试
"""

from swanlab.data.run.metadata.hardware.type import HardwareCollector, CollectGuard, HardwareInfoList
from swanlab.data.run.metadata.hardware.type import HardwareInfo


def test_try_run():
    """
    测试try_run包装器
    """
    data = {"key": "test", "value": 1, "name": "test", "config": None}

    class TestTryRun(HardwareCollector):
        @HardwareCollector.try_run()
        def collect(self) -> HardwareInfo:
            return data

    t = TestTryRun()
    assert t.collect() == data

    class TestErrorTryRun(HardwareCollector):
        @HardwareCollector.try_run()
        def collect(self) -> HardwareInfo:
            raise Exception("test")

    t = TestErrorTryRun()
    assert t() is None


def test_hardware():
    data = {"key": "test", "value": 1, "name": "test", "config": None}

    class TestCollector(HardwareCollector):
        def collect(self) -> HardwareInfo:
            return data

    t = TestCollector()
    assert t.collect() == data

    class TestErrorCollector(HardwareCollector):
        def collect(self) -> HardwareInfo:
            raise Exception("test")

    t = TestErrorCollector()
    assert t() is None

    # 采集任务中有部分采集失败，但是不影响其他采集任务
    class TestErrorCollector(HardwareCollector):
        @staticmethod
        def collect_first() -> HardwareInfo:
            return data

        @HardwareCollector.try_run()
        def collect_second(self) -> HardwareInfo:
            raise Exception("test")

        def collect(self) -> HardwareInfoList:
            return [
                self.collect_first(),
                self.collect_second(),
            ]

    t = TestErrorCollector()
    assert t.collect() == [data, None]
    assert t() == [data]


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
    assert t.after_count == 0
    assert t.collect_num == 1
    # 60>count>0, 不执行
    t.collect_num = 30
    t.before_collect()
    t.after_collect()
    assert t.before_count == 1
    assert t.after_count == 0
    assert t.collect_num == 31
    # 60>count>0, 不执行
    t.collect_num = 59
    t.before_collect()
    t.after_collect()
    assert t.before_count == 1
    assert t.after_count == 1
    assert t.collect_num == 60
    # count=60，不执行，但是after执行
    t.before_collect()
    t.after_collect()
    assert t.before_count == 2
    assert t.after_count == 2
    assert t.collect_num == 61
