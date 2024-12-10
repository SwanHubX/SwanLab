"""
@author: cunyue
@file: test_utils.py
@time: 2024/12/5 13:27
@description: 硬件信息采集工具测试
"""

from swanlab.data.run.metadata.hardware.type import HardwareCollector
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
