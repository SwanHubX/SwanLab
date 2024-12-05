"""
@author: cunyue
@file: test_utils.py
@time: 2024/12/5 13:27
@description: 硬件信息采集工具测试
"""

from swanlab.data.run.metadata.hardware.type import HardwareInfo
from swanlab.data.run.metadata.hardware.utils import hardware


def test_hardware():

    @hardware
    def func() -> HardwareInfo:
        return {"value": 1, "name": "test", "key": "12345"}

    assert func() == {"key": "12345", "value": 1, "name": "test"}


def test_hardware_err():

    @hardware
    def func() -> HardwareInfo:
        raise Exception("test error")

    assert func() is None
