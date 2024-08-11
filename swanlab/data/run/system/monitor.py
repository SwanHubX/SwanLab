#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-29 19:21:38
@File: swanlab\system\monitor.py
@IDE: vscode
@Description:
    监控硬件数据类，包含16种硬件参数
"""
from pynvml import *
import psutil
import time
import asyncio


class SwanSystemMonitor:
    """监控硬件数据类
    使用方法：
        m = SwanSystemMonitor()
        m.get_all()
    返回值为一个dict，包含16种硬件参数：
        gpu: GPU信息，包含GPU数量，每个GPU的显存使用情况，显存剩余情况，显存使用率，功耗，温度，负载率
        network: 网络IO数据，包含接收速度和发送速度
        disk_io: 磁盘IO数据，包含读速度和写速度
        disk_usage: 磁盘占用率
        memory: 内存数据，包含内存剩余量，使用量和使用率
        system_cpu_usage: 系统CPU使用率
        process_cpu_usage: 当前进程CPU使用率
        system_temperature: 当前系统所有传感器的温度
        system_cpu_freq: 当前系统CPU的频率
        timestamp: 当前时间
    """

    def __init__(self) -> None:
        """初始化nvml实例，并计算GPU数量"""
        try:
            nvmlInit()  # 初始化
            self.is_gpu_available = True
            self.gpu_device_count = nvmlDeviceGetCount()
        except:
            self.is_gpu_available = False
            self.gpu_device_count = 0

    async def run():
        """异步等待，避免阻塞事件循环"""
        while True:
            await asyncio.sleep(1)

    def get_gpu_info(self):
        """将当前时间的GPU信息全部收集，返回1个字典，字典的键为第几个GPU，值为该GPU的信息"""
        if self.is_gpu_available:
            gpu_info = dict()
            gpu_info["gpu_count"] = self.gpu_device_count
            for i in range(self.gpu_device_count):
                gpu_info[i] = {
                    "gpu_mem_usage_gbs": self.get_gpu_mem_usage(i),
                    "gpu_mem_free_gbs": self.get_gpu_mem_free(i),
                    "gpu_mem_utilization": self.get_gpu_mem_utilization(i),
                    "gpu_power": self.get_gpu_power(i),
                    "gpu_temperature": self.get_gpu_temperature(i),
                    "gpu_utilization": self.get_gpu_utilization(i),
                }

            return gpu_info
        else:
            return None

    def get_gpu_mem_usage(self, gpu_index):
        """得到当前时间GPU显存的使用情况，单位为GB，可观测到MB"""
        handle = nvmlDeviceGetHandleByIndex(gpu_index)
        gpu_mem_usage = f"{nvmlDeviceGetMemoryInfo(handle).used/1024**3:.3f}"

        return gpu_mem_usage

    def get_gpu_mem_free(self, gpu_index):
        """得到当前时间GPU显存的剩余情况，单位为GB，可观测到MB"""
        handle = nvmlDeviceGetHandleByIndex(gpu_index)
        gpu_mem_free = f"{nvmlDeviceGetMemoryInfo(handle).free/1024**3:.3f}"

        return gpu_mem_free

    def get_gpu_mem_utilization(self, gpu_index):
        """得到当前时间GPU内存的使用率, 单位为%,精确到小数点后三位"""
        handle = nvmlDeviceGetHandleByIndex(gpu_index)
        gpu_mem_utilization = f"{nvmlDeviceGetMemoryInfo(handle).used*100/nvmlDeviceGetMemoryInfo(handle).total:.3f}"

        return gpu_mem_utilization

    def get_gpu_power(self, gpu_index):
        """得到当前时间GPU的功耗情况, 单位为瓦"""
        handle = nvmlDeviceGetHandleByIndex(gpu_index)
        gpu_power = f"{nvmlDeviceGetPowerState(handle)/ 1000.0}"

        return gpu_power

    def get_gpu_temperature(self, gpu_index):
        """得到当前时间GPU的温度，返回字典，键：第几个GPU，值：温度 单位为摄氏度"""
        handle = nvmlDeviceGetHandleByIndex(gpu_index)
        gpu_temperature = f"{nvmlDeviceGetTemperature(handle,NVML_TEMPERATURE_GPU)}"

        return gpu_temperature

    def get_gpu_utilization(self, gpu_index):
        """得到当前时间GPU的负载率，返回字典，键：第几个GPU，值：负载率 单位为%"""
        handle = nvmlDeviceGetHandleByIndex(gpu_index)
        gpu_utilization = f"{nvmlDeviceGetUtilizationRates(handle).gpu}"

        return gpu_utilization

    def get_network_io(self):
        """得到当前时间网络IO数据，单位为KB，可观测到B"""

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
        network_receive_kbs = f"{get_net_io_speed()[0]/1024**3:.3f}"
        network_send_kbs = f"{get_net_io_speed()[1]/1024**3:.3f}"

        return {"network_receive_kbs": network_receive_kbs, "network_send_kbs": network_send_kbs}

    def get_disk_io(self):
        """得到当前时间磁盘读写数据，单位为MB，可观测到KB"""

        disk_io_msg = psutil.disk_io_counters()
        disk_read_mbs = f"{disk_io_msg.read_bytes/1024**2:.3f}"
        disk_write_mbs = f"{disk_io_msg.write_bytes/1024**2:.3f}"

        return {"disk_read_mbs": disk_read_mbs, "disk_write_mbs": disk_write_mbs}

    def get_disk_usage(self):
        """得到当前时间磁盘占用率，单位为%，精确到小数点后三位"""
        disk_usage = psutil.disk_usage("/")
        disk_used_percent = f"{disk_usage.percent:.3f}"

        return disk_used_percent

    def get_memory_usage(self):
        """得到当前时间内存剩余量，使用量和使用率，单位为GB，可观测到MB，使用率单位为%，精确到小数点后3位"""

        memory = psutil.virtual_memory()
        memory_free_gbs = f"{memory.available / (1024 ** 3):.3f}"
        memory_used_gbs = f"{memory.used / (1024 ** 3):.3f}"
        memory_utilization = f"{memory.percent:.3f}"

        return {
            "memory_free_gbs": memory_free_gbs,
            "memory_used_gbs": memory_used_gbs,
            "memory_utilization": memory_utilization,
        }

    def get_system_cpu_usage(self):
        """得到当前时间系统CPU的使用率，单位为%，精确到小数点后一位"""
        system_cpu_usage = psutil.cpu_percent(interval=1)

        return system_cpu_usage

    def get_process_cpu_usage(self):
        """得到当前时间当前进程的CPU使用率，单位为%，精确到小数点后一位"""
        current_process = psutil.Process()
        process_cpu_usage = current_process.cpu_percent(interval=1)

        return process_cpu_usage


    def get_system_cpu_freq(self):
        """
        得到当前时间系统CPU的频率
        ---
        信息:
            单位:   Mhz
            精确度: 不确定
            可用性: Linux,macOS,Windows,FreeBSD,OpenBSD
        
        WARNING:在linux上,current获取实时值,其他平台上则是固定值
        """
        cpu_freq = psutil.cpu_freq(percpu=False)
        return cpu_freq

    def get_system_temperature(self):
        """
        得到当前时间系统各传感器的温度
        ---
        信息:
            单位:   摄氏度℃
            精确度: 小数点后两位
            可用性: Linux,FreeBSD
        
        Return:
            temperatures: {"device1":[],"device2":[],...}
        
        TODO:传感器可能的返回参数https://www.kernel.org/doc/html/latest/subsystem-apis.html
        已知:k10temp(AMD CPU 10th~16th Opteron~zen3)
        """
        if not hasattr(psutil, "sensors_temperatures"):return "platform not supported" 
        temperatures = psutil.sensors_temperatures()
        if not temperatures:return "can't read any temperature" 
        return temperatures
                    
    def get_all(self):
        """获取全部硬件数据"""
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

        return {
            "gpu": self.get_gpu_info(),
            "network": self.get_network_io(),
            "disk_io": self.get_disk_io(),
            "disk_usage": self.get_disk_usage(),
            "memory": self.get_memory_usage(),
            "system_cpu_usage": self.get_system_cpu_usage(),
            "process_cpu_usage": self.get_process_cpu_usage(),
            "system_temperature":self.get_system_temperature(),
            'system_cpu_freq':self.get_system_cpu_freq(),
            "timestamp": timestamp,
        }


if __name__ == "__main__":
    # 测试脚本
    m = SwanSystemMonitor()
    # 获取当前时间的全部硬件数据
    # 例如：{'gpu': None, 'network': {'network_receive_kbs': '0.000', 'network_send_kbs': '0.000'}, 'disk_io': {'disk_read_mbs': '224801.238', 'disk_write_mbs': '199629.906'}, 'disk_usage': '2.700', 'memory': {'memory_free_gbs': '28.493', 'memory_used_gbs': '30.112', 'memory_utilization': '55.500'}, 'system_cpu_usage': 26.4, 'process_cpu_usage': 0.0, 'timestamp': '2024-07-29 01:52:05'}
    print(m.get_all())
