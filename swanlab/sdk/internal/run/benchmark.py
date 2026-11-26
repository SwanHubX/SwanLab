"""
@author: mx
@file: benchmark.py
@time: 2026/11/25
@description: 评测项对象，负责声明评测项并上报分数
"""

from __future__ import annotations

import copy
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, Dict

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.sdk.internal.bus.events import ScoreLogEvent
from swanlab.sdk.internal.pkg import console, fork

if TYPE_CHECKING:
    from .eval import EvalRun

__all__ = ["Benchmark"]

# 评测项名称的长度上限，与后端字段约束一致
MAX_BENCHMARK_NAME_LENGTH = 512


class Benchmark:
    """一个评测项，承载该评测项的局部配置与分数上报。

    由 ``EvalRun.benchmark()`` 创建，评测项名称作为分数记录的 dataset 维度。
    """

    def __init__(self, run: EvalRun, name: str, config: Dict[str, Any]):
        self._run = run
        self._name = name
        self._config = config

    @property
    def name(self) -> str:
        """评测项名称。"""
        return self._name

    @property
    def config(self) -> Dict[str, Any]:
        """评测项的局部配置。"""
        return copy.deepcopy(self._config)

    def log(self, metric: Mapping[str, Any]) -> None:
        """上报该评测项的一组指标。

        :param metric: 指标名到数值的映射，每个键值对展开为一条分数记录。
        """
        if fork.is_forked(self._run._init_pid):
            raise RuntimeError(
                "SwanLab Eval does not support fork yet. Use `multiprocessing.set_start_method('spawn')` "
                "or call `swanlab.eval()` in the child process."
            )
        if not isinstance(metric, Mapping) or not metric:
            console.error(
                f"Score data must be a non-empty dict, but got {type(metric).__name__}. SwanLab will ignore this log."
            )
            return
        timestamp = Timestamp()
        timestamp.GetCurrentTime()
        with self._run._api_lock:
            if not self._run.alive:
                raise RuntimeError("`benchmark.log()` requires an active Eval, call `swanlab.eval()` first.")
            self._run._components.emitter.emit(
                ScoreLogEvent(dataset=self._name, data=dict(metric), timestamp=timestamp)
            )
