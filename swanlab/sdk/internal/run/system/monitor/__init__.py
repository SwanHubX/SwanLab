"""
@author: cunyue
@file: __init__.py
@time: 2026/4/9 16:41
@description: SwanLab 硬件监控对象
"""

import math
from concurrent.futures.thread import ThreadPoolExecutor
from typing import Dict, List, Optional

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.sdk.internal.bus.emitter import RunEmitter
from swanlab.sdk.internal.bus.events import MetricLogEvent, ScalarDefineEvent
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.pkg.timer import Timer
from swanlab.sdk.internal.run.system.hardware_vendor.apple import Apple
from swanlab.sdk.internal.run.system.hardware_vendor.cpu import CPU
from swanlab.sdk.internal.run.system.hardware_vendor.memory import Memory
from swanlab.sdk.typings.run.system import SystemScalars, SystemShim
from swanlab.sdk.typings.run.system.hardware_vendor import CollectorProtocol, CollectResult
from swanlab.utils.experiment import generate_id

__all__ = ["Monitor", "fmt_system_key", "is_system_key"]


_SYSTEM_KEY_PREFIX = "__swanlab__."


def fmt_system_key(key: str):
    """
    格式化系统指标key，如果已经是系统指标key，报错
    """
    if key.startswith(_SYSTEM_KEY_PREFIX):
        raise ValueError(f"System metric key '{key}' is already a system metric key")
    return f"{_SYSTEM_KEY_PREFIX}{key}"


def is_system_key(key: str):
    """
    判断是否是系统指标key
    :param key: 指标key
    :return: 是否是系统指标key
    """
    return key.startswith(_SYSTEM_KEY_PREFIX)


class Monitor:
    """
    硬件监控定时任务
    """

    def __init__(self, shim: SystemShim):
        self._timer: Optional[Timer] = None
        self._executor: Optional[ThreadPoolExecutor] = None
        self._shim = shim

    @property
    def is_running(self):
        return self._timer is not None and self._timer.is_running

    def start(self, ctx: RunContext, emitter: RunEmitter):
        # 1. 收集采集器
        # 注意 scalars 设计上越靠前的越先发送、在前端显示越靠前
        collectors: List[CollectorProtocol] = []
        scalars: SystemScalars = []
        # 1.1 基础硬件信息
        if (cpu := CPU.new(self._shim)) is not None:
            collectors.append(cpu[0])
            scalars = cpu[1] + scalars
        if (memory := Memory.new(self._shim)) is not None:
            collectors.append(memory[0])
            scalars = memory[1] + scalars

        # 1.2 苹果芯片信息
        if (apple := Apple.new(self._shim)) is not None:
            collectors.append(apple[0])
            scalars = apple[1] + scalars

        # 1.3 加速器信息

        # 3. 定义指标
        # 有部分指标会聚合到一个图表中，所以需要一个缓存来记录每个指标对应的图表索引
        cache_chart_index: Dict[str, str] = {}
        for scalar in scalars:
            key = fmt_system_key(scalar.key)
            if (chart_index := cache_chart_index.get(scalar.chart_name)) is None:
                chart_index = generate_id(8)
                cache_chart_index[scalar.chart_name] = chart_index
            emitter.emit(
                ScalarDefineEvent(
                    key=key,
                    name=scalar.name,
                    color=scalar.color,
                    system=True,
                    x_axis=scalar.x_axis,
                    chart=chart_index,
                    chart_name=scalar.chart_name,
                )
            )
        # 4. 定义任务
        all_handlers = [(type(c).__name__, c.collect) for c in collectors]
        if len(all_handlers) == 0:
            console.debug("No hardware monitor collectors found, skipping creating monitor task")
            return
        # 使用线程池执行任务，并发数为采集器数量的1/4，最多4个线程
        self._executor = ThreadPoolExecutor(max_workers=min(max(math.ceil(len(collectors) / 4), 1), 4))

        def task():
            assert self._executor is not None, "Monitor Executor is not initialized"
            futures = [(n, self._executor.submit(fn)) for n, fn in all_handlers]
            results: List[CollectResult] = []
            for n, f in futures:
                try:
                    result = f.result()
                    results.extend(result)
                except Exception as e:
                    console.error(f"Error collecting metric via {n}: {e}")
            ts = Timestamp()
            ts.GetCurrentTime()
            step = ctx.metrics.next_system_step()
            data = {fmt_system_key(k): v for k, v in results}
            emitter.emit(MetricLogEvent(step=step, data=data, timestamp=ts))

        self._timer = Timer(task, interval=ctx.config.settings.monitor.interval, immediate=True, name="SwanLab·Monitor")
        self._timer.start()

    def stop(self) -> None:
        assert self._timer is not None, "HardwareMonitor is not running"
        self._timer.cancel()
        self._timer.join()
        if self._executor is not None:
            self._executor.shutdown(wait=False)
            self._executor = None
