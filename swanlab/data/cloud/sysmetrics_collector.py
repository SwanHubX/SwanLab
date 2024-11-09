#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/10/31 21:30
@File: sysmetrics_collector.py
@IDE: pycharm
@Description:
    系统指标收集和上传器
@Author: kites262
"""
import os
from typing import Optional

import psutil
import GPUtil

from swanlab.data.modules import DataWrapper, Line
from swanlab.data.run.main import get_run
from .utils import ThreadTaskABC, ThreadUtil
from ..run.exp import SwanLabExp


class SysmetricsCollectorTask(ThreadTaskABC):
    """
    系统指标收集器，负责收集系统指标并上传到实验实例
    """

    PERIOD_TIME = 60

    def __init__(self):
        self.exp_instance: Optional[SwanLabExp] = None
        self.gpus = None
        self.metrics_dict = {}
        self.process = psutil.Process(os.getpid())
        pass

    # 任务函数，会被周期性调用
    def task(self, u: ThreadUtil, *args):
        if u.timer.can_run(self.PERIOD_TIME, cancel=False):
            self.exp_instance = args[0]
            self.record_and_upload_sysmetrics()

    def record_and_upload_sysmetrics(self):
        # 收集指标
        self.record_os_metrics()
        self.record_swanlab_metrics()
        self.record_gpu_metrics()

        # 上传指标
        self.upload_metrics_dict()

    # 添加一条记录到指标字典，等待上传
    def add_metric(self, key, value):
        prefixed_key = f'__swanlab/{key}'
        self.metrics_dict[prefixed_key] = value

    # 上传指标字典，将记录的所有指标添加到实验实例中，统一上传
    def upload_metrics_dict(self):
        for key, value in self.metrics_dict.items():
            self.exp_instance.add(key=key, data=DataWrapper(key, [Line(value)]))
        self.metrics_dict.clear()

    # 记录操作系统相关指标
    def record_os_metrics(self):
        cpu_percent = psutil.cpu_percent()
        cpus_percent = psutil.cpu_percent(percpu=True)
        net_io = psutil.net_io_counters()
        self.add_metric('OS.CPU (%)', cpu_percent)
        self.add_metric('OS.UsedMemory (%)', psutil.virtual_memory().percent)
        self.add_metric('OS.UsedMemory (MB)', psutil.virtual_memory().used >> 20)
        self.add_metric('OS.NetIO.BytesSent (MB)', net_io.bytes_sent >> 20)
        self.add_metric('OS.NetIO.BytesReceived (MB)', net_io.bytes_recv >> 20)
        for i, cpu_percent in enumerate(cpus_percent):
            self.add_metric(f'OS.CPU-{i} (%)', cpu_percent)

    # 记录 GPU 相关指标
    def record_gpu_metrics(self):
        self.gpus = GPUtil.getGPUs()
        for gpu in self.gpus:
            self.add_metric(f'GPU-{gpu.id}.Load (%)', gpu.load * 100)
            self.add_metric(f'GPU-{gpu.id}.UsedMemory (MB)', gpu.memoryUsed)
            self.add_metric(f'GPU-{gpu.id}.Temperature (°C)', gpu.temperature)

    # 记录本次 SwanLab Run 对应进程的指标
    # 在训练过程中，SwanLab Run 进程（即Python进程）的指标可以反映训练进程的资源占用情况
    def record_swanlab_metrics(self):
        self.process = psutil.Process(os.getpid())
        process_cpu = self.process.cpu_percent(interval=0.1)
        swanlab_io_counters = self.process.io_counters()
        self.add_metric('swanlab.CPU (%)', process_cpu)
        self.add_metric('swanlab.Memory (%)', self.process.memory_percent())
        self.add_metric('swanlab.Memory (MB)', self.process.memory_info().rss >> 20)
        self.add_metric('swanlab.IO.Read (MB)', swanlab_io_counters.read_bytes >> 20)
        self.add_metric('swanlab.IO.Write (MB)', swanlab_io_counters.write_bytes >> 20)

    # 主线程结束时，触发此回调
    def callback(self, u: ThreadUtil, *args):
        pass
