"""
@author: cunyue
@file: network.py
@time: 2024/12/10 14:30
@description: 网络信息采集
"""

import time
from typing import List

import psutil

from .type import HardwareFuncResult, HardwareCollector, HardwareInfo
from .utils import generate_key, random_index, HardwareConfig


def get_network_info() -> HardwareFuncResult:
    """获取网络信息"""
    return None, NetworkCollector()


class NetworkCollector(HardwareCollector):
    """网络信息采集器"""

    def __init__(self):
        super().__init__()
        self.last_net_io = psutil.net_io_counters()
        self.last_time = time.time()
        network_base_config = HardwareConfig(
            y_range=(0, None),
            chart_name="Network Traffic (KB)",
            chart_index=random_index(),
        ).clone()

        # 网络发送速度
        self.sent_key = generate_key("network.sent")
        self.sent_config = network_base_config.clone(metric_name="sent")

        # 网络接收速度
        self.recv_key = generate_key("network.recv")
        self.recv_config = network_base_config.clone(metric_name="received")

        self.unit = 1024

    def collect(self) -> List[HardwareInfo]:
        current_net_io = psutil.net_io_counters()
        current_time = time.time()  # noqa
        time_diff = current_time - self.last_time

        results = [
            self.get_network_sent_speed(current_net_io.bytes_sent, time_diff),
            self.get_network_recv_speed(current_net_io.bytes_recv, time_diff),
        ]

        # 更新上次的值
        self.last_net_io = current_net_io
        self.last_time = current_time
        return results

    def get_network_sent_speed(self, bytes_sent: int, time_diff: float) -> HardwareInfo:
        """
        获取网络发送速度 (KB/s)
        """
        diff = max(0, bytes_sent - self.last_net_io.bytes_sent)
        return {
            "key": self.sent_key,
            "name": "Network Traffic Sent (KB)",
            "value": self.division_guard(diff, time_diff * self.unit),
            "config": self.sent_config,
        }

    def get_network_recv_speed(self, bytes_recv: int, time_diff: float) -> HardwareInfo:
        """
        获取网络接收速度 (KB/s)
        """
        diff = max(0, bytes_recv - self.last_net_io.bytes_recv)
        return {
            "key": self.recv_key,
            "name": "Network Traffic Received (KB)",
            "value": self.division_guard(diff, time_diff * self.unit),
            "config": self.recv_config,
        }
