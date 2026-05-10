"""
@author: cunyue
@file: __init__.py
@time: 2026/4/9 16:41
@description: SwanLab 硬件监控对象
"""

from concurrent.futures.thread import ThreadPoolExecutor
from typing import Dict, List, Optional

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScalarRecord
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console, safe, timer
from swanlab.sdk.internal.probe_python.hardware_vendor.apple import Apple
from swanlab.sdk.internal.probe_python.hardware_vendor.cpu import CPU
from swanlab.sdk.internal.probe_python.hardware_vendor.memory import Memory
from swanlab.sdk.internal.probe_python.monitor import builder
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
        now_step = ctx.global_system_step
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
        column_records: List[ColumnRecord] = []
        for scalar in scalars:
            if (chart_index := cache_chart_index.get(scalar.chart_name)) is None:
                chart_index = generate_id(8)
                cache_chart_index[scalar.chart_name] = chart_index
            with safe.block(message=f"Failed to build column record for scalar {scalar.key}"):
                column_records.append(builder.build_probe_column(scalar, chart_index=chart_index))
        if len(column_records) == 0:
            console.debug("No hardware monitor columns found, skipping creating monitor task")
            return None
        self._core.upsert_columns(column_records)
        # 3. 定义并启动任务
        all_handlers = [(type(c).__name__, c.collect) for c in collectors]
        if len(all_handlers) == 0:
            console.debug("No hardware monitor collectors found, skipping creating monitor task")
            return None
        # 使用线程池执行任务
        self._executor = ThreadPoolExecutor(max_workers=2)

        def task():
            nonlocal now_step
            assert self._executor is not None, "Monitor Executor is not initialized"
            futures = [(n, self._executor.submit(fn)) for n, fn in all_handlers]
            results: List[CollectResult] = []
            for n, f in futures:
                with safe.block(message=f"Error collecting metric via {n}"):
                    result = f.result()
                    results.extend(result)
            ts = Timestamp()
            ts.GetCurrentTime()
            scalar_records: List[ScalarRecord] = []
            now_step += 1
            for k, v in results:
                scalar_records.append(builder.build_probe_scalar(k, value=v, timestamp=ts, step=now_step))
            self._core.upsert_scalars(scalar_records)

        # 不设置立即执行，以避免产生一些无用的数据
        self._timer = timer.Timer(
            task,
            interval=ctx.config.settings.probe.monitor_interval,
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
