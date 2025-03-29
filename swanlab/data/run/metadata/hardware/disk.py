"""
@author: cunyue
@file: disk.py
@time: 2024/12/3 20:12
@description: 磁盘信息采集
"""

import time
from typing import List

import psutil

from swanlab.swanlab_settings import get_settings
from .type import HardwareFuncResult, HardwareCollector, HardwareInfo, HardwareConfig
from .utils import random_index, generate_key


def get_disk_info() -> HardwareFuncResult:
    """获取磁盘信息"""
    return None, DiskCollector()


class DiskCollector(HardwareCollector):
    def __init__(self):
        super().__init__()
        self.last_disk_io = psutil.disk_io_counters()
        self.last_time = time.time()

        self.disk_path = get_settings().disk_io_dir

        disk_base_config = HardwareConfig(
            y_range=(0, None),
            chart_name="Disk I/O Utilization (MB)",
            chart_index=random_index(),
        ).clone()

        # 磁盘读取速度
        self.read_key = generate_key("disk.read")
        self.read_config = disk_base_config.clone(metric_name="read")

        # 磁盘写入速度
        self.write_key = generate_key("disk.write")
        self.write_config = disk_base_config.clone(metric_name="write")

        # 磁盘使用率
        self.usage_key = generate_key("disk.usage")
        self.usage_config = HardwareConfig(
            y_range=(0, 100),
            chart_name="Disk Utilization (%)",
        ).clone()

        self.unit = 1024**2

    def get_disk_read_speed(self, read_bytes: int, time_diff: float) -> HardwareInfo:
        """
        获取磁盘读取速度 (MB/s)
        """
        diff = read_bytes - self.last_disk_io.read_bytes
        return {
            "key": self.read_key,
            "name": "MB read from disk",
            "value": self.division_guard(diff, time_diff * self.unit),
            "config": self.read_config,
        }

    def get_disk_write_speed(self, write_bytes: int, time_diff: float) -> HardwareInfo:
        """
        获取磁盘写入速度 (MB/s)
        """
        diff = write_bytes - self.last_disk_io.write_bytes
        return {
            "key": self.write_key,
            "name": "MB written to disk",
            "value": self.division_guard(diff, time_diff * self.unit),
            "config": self.write_config,
        }

    def get_disk_usage(self) -> HardwareInfo:
        """
        获取磁盘使用率
        """
        result = psutil.disk_usage(self.disk_path).percent
        return {
            "key": self.usage_key,
            "name": "Disk Utilization (%)",
            "value": result,
            "config": self.usage_config,
        }

    def collect(self) -> List[HardwareInfo]:
        current_disk_io = psutil.disk_io_counters()
        current_time = time.time()
        time_diff = current_time - self.last_time

        disk_result = [
            self.get_disk_read_speed(current_disk_io.read_bytes, time_diff),
            self.get_disk_write_speed(current_disk_io.write_bytes, time_diff),
            self.get_disk_usage(),
        ]
        # 更新上次的值
        self.last_disk_io = current_disk_io
        self.last_time = current_time
        return disk_result
