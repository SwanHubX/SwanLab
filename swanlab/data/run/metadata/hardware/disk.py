"""
@author: cunyue
@file: disk.py
@time: 2024/12/3 20:12
@description: 磁盘信息采集
"""
from typing import List
import psutil

from swanlab.data.run.metadata.hardware.type import HardwareFuncResult, HardwareCollector, HardwareInfo
from .utils import DiskCollector as D


def get_disk_info() -> HardwareFuncResult:
    """获取磁盘信息"""
    return None, DiskCollector()


class DiskCollector(HardwareCollector, D):
    def __init__(self):
        super().__init__()
        self.last_disk_io = psutil.disk_io_counters()
        self.last_time = psutil.time.time()
    
    def collect(self) -> List[HardwareInfo]:
        current_disk_io = psutil.disk_io_counters()
        current_time = psutil.time.time()
        
        # 计算时间差
        time_diff = current_time - self.last_time
        
        # 防止除零错误
        if time_diff <= 0:
            read_speed = 0
            write_speed = 0
        else:
            # 计算读写速度 (bytes/s)
            read_bytes_diff = current_disk_io.read_bytes - self.last_disk_io.read_bytes
            write_bytes_diff = current_disk_io.write_bytes - self.last_disk_io.write_bytes
            
            read_speed = read_bytes_diff / time_diff
            write_speed = write_bytes_diff / time_diff
        
        # 更新上次的值
        self.last_disk_io = current_disk_io
        self.last_time = current_time
        
        return [
            self.get_disk_read_speed(read_speed),
            self.get_disk_write_speed(write_speed),
            self.get_disk_usage(psutil.disk_usage('/').percent)
        ]