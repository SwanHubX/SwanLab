"""
@author: cunyue
@file: utils.py
@time: 2024/12/5 13:09
@description: 硬件信息采集工具
"""

from swanlab.log import swanlog
from .type import HardwareMonitorFunc


def hardware(func: HardwareMonitorFunc) -> HardwareMonitorFunc:
    """
    硬件信息采集函数装饰器
    如果函数执行失败返回None
    """

    def wrapper():
        try:
            return func()
        except Exception as e:
            swanlog.error("Hardware info collection failed: %s, %s", func.__name__, str(e))
            return None

    return wrapper
