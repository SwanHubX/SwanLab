from pynvml import *
import psutil
import time
import asyncio


class Monitor:
    """监控硬件数据类"""

    deviceCount: int

    def __init__(self) -> None:
        """初始化nvml实例，并计算GPU数量"""

        nvmlInit()  # 初始化
        self.deviceCount = nvmlDeviceGetCount()

    async def run():
        """异步等待，避免阻塞事件循环"""

        while True:
            await asyncio.sleep(1)

    def get_gpu_mem_usage(self):
        """监控GPU内存的使用情况，返回字典，键：第几个GPU，值：内存使用情况 单位为GB，可观测到MB"""

        self.gpu_mem_usage_gbs = dict()
        for i in range(self.deviceCount):
            handle = nvmlDeviceGetHandleByIndex(i)
            info = nvmlDeviceGetMemoryInfo(handle)
            self.gpu_mem_usage_gbs[f"gpu_{i}_mem_usage_gb"] = f"{info.used/1024**3:.3f}"

    def get_gpu_mem_free(self):
        """监控GPU内存的剩余情况，返回字典，键：第几个GPU，值：内存剩余情况 单位为GB，可观测到MB"""

        self.gpu_mem_free_gbs = dict()
        for i in range(self.deviceCount):
            handle = nvmlDeviceGetHandleByIndex(i)
            info = nvmlDeviceGetMemoryInfo(handle)
            self.gpu_mem_free_gbs[f"gpu_{i}_mem_free_gb"] = f"{info.free/1024**3:.3f}"

    def get_gpu_mem_utilization(self):
        """监控GPU内存的使用率，返回字典，键：第几个GPU，值：内存的使用百分比 单位为%,精确到小数点后三位"""
        self.gpu_mem_utilization = dict()
        for i in range(self.deviceCount):
            handle = nvmlDeviceGetHandleByIndex(i)
            info = nvmlDeviceGetMemoryInfo(handle)
            self.gpu_mem_utilization[f"gpu_{i}_mem_utilization"] = f"{info.used*100/info.total:.3f}"

    def get_gpu_power(self):
        """监控GPU的功耗情况，返回字典，键：第几个GPU，值：功耗，单位为瓦"""

        self.gpu_power = dict()
        for i in range(self.deviceCount):
            handle = nvmlDeviceGetHandleByIndex(i)
            self.gpu_power[f"gpu_{i}_power"] = f"{nvmlDeviceGetPowerState(handle)/ 1000.0}"

    def get_gpu_temperature(self):
        """监控GPU的温度，返回字典，键：第几个GPU，值：温度 单位为摄氏度"""

        self.gpu_temperature = dict()
        for i in range(self.deviceCount):
            handle = nvmlDeviceGetHandleByIndex(i)
            self.gpu_temperature[f"gpu_{i}_temperature"] = f"{nvmlDeviceGetTemperature(handle,NVML_TEMPERATURE_GPU)}"

    def get_gpu_utilization(self):
        """监控GPU的负载率，返回字典，键：第几个GPU，值：负载率 单位为%"""

        self.gpu_utilization = dict()
        for i in range(self.deviceCount):
            handle = nvmlDeviceGetHandleByIndex(i)
            self.gpu_utilization[f"gpu_{i}_utilization"] = f"{nvmlDeviceGetUtilizationRates(handle).gpu}"

    def get_network_io(self):
        """监控网络IO数据，单位为KB，可观测到B"""

        def get_net_io_speed(interval=1):
            # 获取初始时间和网络IO统计数据

            initial_net_io = psutil.net_io_counters()

            initial_time = time.time()

            # 等待一段时间

            time.sleep(interval)

            # 再次获取当前时间和网络IO统计数据

            current_net_io = psutil.net_io_counters()
            current_time = time.time()
            # 计算时间间隔

            time_elapsed = current_time - initial_time

            # 计算接收和发送的字节数

            bytes_recv = current_net_io.bytes_recv - initial_net_io.bytes_recv

            bytes_sent = current_net_io.bytes_sent - initial_net_io.bytes_sent

            # 计算速度（字节/秒）

            recv_speed = bytes_recv / time_elapsed

            sent_speed = bytes_sent / time_elapsed

            return recv_speed, sent_speed

        # 获取网络接收和发送速度

        self.network_receive_kbs = f"{get_net_io_speed()[0]/1024**3:.3f}"
        self.network_send_kbs = f"{get_net_io_speed()[1]/1024**3:.3f}"

    def get_disk_io(self):
        """监控磁盘读写数据，单位为MB，可观测到KB"""

        disk_io_msg = psutil.disk_io_counters()
        self.disk_read_mbs = f"{disk_io_msg.read_bytes/1024**2:.3f}"
        self.disk_write_mbs = f"{disk_io_msg.write_bytes/1024**2:.3f}"

    def get_disk_usage(self):
        """监控磁盘占用率，单位为%，精确到小数点后三位"""
        disk_usage = psutil.disk_usage("/")
        self.disk_used_percent = f"{disk_usage.percent:.3f}"

    def get_memory_usage(self):
        """获取内存剩余量，使用量和使用率，单位为GB，可观测到MB，使用率单位为%，精确到小数点后3位"""

        memory = psutil.virtual_memory()
        self.memory_free_gbs = f"{memory.available / (1024 ** 3):.3f}"
        self.memory_used_gbs = f"{memory.used / (1024 ** 3):.3f}"
        self.memory_utilization = f"{memory.percent:.3f}"

    def get_system_cpu(self):
        """监控系统CPU的使用率，单位为%，精确到小数点后一位"""

        self.system_cpu_usage = psutil.cpu_percent(interval=1)

    def get_process_cpu(self):
        """监控进程的CPU使用率，单位为%，精确到小数点后一位"""

        current_process = psutil.Process()
        self.process_cpu_usage = current_process.cpu_percent(interval=1)
