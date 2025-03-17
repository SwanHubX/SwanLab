"""
@author: cunyue
@file: network.py
@time: 2024/12/10 14:30
@description: 网络信息采集
"""
from typing import List

import psutil

from swanlab.data.run.metadata.hardware.type import HardwareFuncResult, HardwareCollector, HardwareInfo
from .utils import generate_key, random_index, HardwareConfig


def get_network_info() -> HardwareFuncResult:
    """获取网络信息"""
    return None, NetworkCollector()


class NetworkCollector(HardwareCollector):
    """网络信息采集器"""
    
    NETWORK_BASE_CONFIG = HardwareConfig(
        y_range=(0, None),
        chart_name="Network Traffic (KB)",
        chart_index=random_index(),
    ).clone()
    
    # 网络发送速度
    NETWORK_SENT_KEY = generate_key("network.sent")
    NETWORK_SENT_CONFIG = NETWORK_BASE_CONFIG.clone(metric_name="sent")
    
    # 网络接收速度
    NETWORK_RECV_KEY = generate_key("network.recv")
    NETWORK_RECV_CONFIG = NETWORK_BASE_CONFIG.clone(metric_name="received")
    
    def __init__(self):
        super().__init__()
        self.last_net_io = psutil.net_io_counters()
        self.last_time = psutil.time.time()
    
    def collect(self) -> List[HardwareInfo]:
        current_net_io = psutil.net_io_counters()
        current_time = psutil.time.time()
        
        # 计算时间差
        time_diff = current_time - self.last_time
        
        # 防止除零错误
        if time_diff <= 0:
            sent_speed = 0
            recv_speed = 0
        else:
            # 计算发送和接收速度 (bytes/s)
            sent_bytes_diff = current_net_io.bytes_sent - self.last_net_io.bytes_sent
            recv_bytes_diff = current_net_io.bytes_recv - self.last_net_io.bytes_recv
            
            sent_speed = sent_bytes_diff / time_diff
            recv_speed = recv_bytes_diff / time_diff
        
        # 更新上次的值
        self.last_net_io = current_net_io
        self.last_time = current_time
        
        return [
            self.get_network_sent_speed(sent_speed),
            self.get_network_recv_speed(recv_speed)
        ]
    
    @staticmethod
    def get_network_sent_speed(speed: float) -> HardwareInfo:
        """
        获取网络发送速度 (KB/s)
        """
        return {
            "key": NetworkCollector.NETWORK_SENT_KEY,
            "name": "Network Sent (KB)",
            "value": speed / 1024,  # 转换为KB/s
            "config": NetworkCollector.NETWORK_SENT_CONFIG,
        }
    
    @staticmethod
    def get_network_recv_speed(speed: float) -> HardwareInfo:
        """
        获取网络接收速度 (KB/s)
        """
        return {
            "key": NetworkCollector.NETWORK_RECV_KEY,
            "name": "Network Received (KB)",
            "value": speed / 1024,  # 转换为KB/s
            "config": NetworkCollector.NETWORK_RECV_CONFIG,
        }
