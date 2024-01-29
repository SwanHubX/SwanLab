from pynvml import *
import psutil


class Monitor:
    """监控硬件数据类"""

    deviceCount: int

    def __init__(self) -> None:
        nvmlInit()  # 初始化
        self.deviceCount = nvmlDeviceGetCount()

    def get_mem_usage(self):
        gpu_0_mem_usage = 0
        return gpu_0_mem_usage
