"""
@author: cunyue
@file: __init__.py
@time: 2026/5/8 20:43
@description: 运行定时心跳，告诉后端当前实验存在运行进程
"""

from functools import partial

from swanlab.sdk.internal.core_python.api.experiment import send_experiment_heartbeat
from swanlab.sdk.internal.pkg.timer import Timer


class Heartbeat:
    def __init__(self, experiment_id: str, interval: int = 10 * 60):
        task = partial(send_experiment_heartbeat, experiment_id=experiment_id)
        # immediate 为 True 时，会在 start 时立即执行一次任务，立即告诉后端有运行进程
        self._timer = Timer(task, interval=interval, immediate=True, name="SwanLab·Heartbeat")

    def start(self):
        self._timer.start()

    def stop(self):
        if self._timer.is_running:
            self._timer.cancel()
            self._timer.join()
