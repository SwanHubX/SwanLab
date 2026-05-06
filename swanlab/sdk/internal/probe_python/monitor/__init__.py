"""
@author: cunyue
@file: __init__.py
@time: 2026/4/9 16:41
@description: SwanLab 硬件监控对象
"""

from concurrent.futures.thread import ThreadPoolExecutor
from typing import Dict, List, Optional

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console, helper, safe, timer
from swanlab.sdk.internal.probe_python.hardware_vendor.apple import Apple
from swanlab.sdk.internal.probe_python.hardware_vendor.cpu import CPU
from swanlab.sdk.internal.probe_python.hardware_vendor.memory import Memory
from swanlab.sdk.protocol import CoreProtocol
from swanlab.sdk.typings.probe_python import SystemScalars, SystemShim
from swanlab.sdk.typings.probe_python.hardware_vendor import CollectorProtocol, CollectResult
from swanlab.utils.experiment import generate_id

__all__ = ["Monitor"]


class Monitor:
    """
    硬件监控定时任务
    """

    def __init__(self, shim: SystemShim, core: CoreProtocol):
        self._core = core
        self._shim = shim
        self._timer: Optional[timer.Timer] = None
        self._executor: Optional[ThreadPoolExecutor] = None

    def start(self, ctx: RunContext) -> Optional["Monitor"]:
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

        # 2. 定义指标
        # 有部分指标会聚合到一个图表中，所以需要一个缓存来记录每个指标对应的图表索引
        cache_chart_index: Dict[str, str] = {}
        for scalar in scalars:
            key = helper.fmt_system_key(scalar.key)
            if (chart_index := cache_chart_index.get(scalar.chart_name)) is None:
                chart_index = generate_id(8)
                cache_chart_index[scalar.chart_name] = chart_index
            _ = (key, chart_index)
            # TODO 向Core发送指标定义
            # emitter.emit(
            #     ScalarDefineEvent(
            #         key=key,
            #         name=scalar.name,
            #         color=scalar.color,
            #         system=True,
            #         x_axis=scalar.x_axis,
            #         chart=chart_index,
            #         chart_name=scalar.chart_name,
            #     )
            # )
        # 3. 定义并启动任务
        all_handlers = [(type(c).__name__, c.collect) for c in collectors]
        if len(all_handlers) == 0:
            console.debug("No hardware monitor collectors found, skipping creating monitor task")
            return None
        # 使用线程池执行任务
        self._executor = ThreadPoolExecutor(max_workers=2)

        def task():
            assert self._executor is not None, "Monitor Executor is not initialized"
            futures = [(n, self._executor.submit(fn)) for n, fn in all_handlers]
            results: List[CollectResult] = []
            for n, f in futures:
                with safe.block(message=f"Error collecting metric via {n}"):
                    result = f.result()
                    results.extend(result)
            ts = Timestamp()
            ts.GetCurrentTime()
            _ = {helper.fmt_system_key(k): v for k, v in results}
            # TODO 向Core发送指标数据
            # emitter.emit(MetricLogEvent(step=step, data=data, timestamp=ts))

        # 不设置立即执行，以避免产生一些无用的数据
        self._timer = timer.Timer(
            task,
            interval=ctx.config.settings.monitor.interval,
            immediate=False,
            name="SwanLab·Monitor",
        )
        self._timer.start()
        return self

    def stop(self) -> None:
        assert self._timer is not None, "HardwareMonitor is not running"
        self._timer.cancel()
        self._timer.join()
        if self._executor is not None:
            self._executor.shutdown(wait=True)
            self._executor = None
