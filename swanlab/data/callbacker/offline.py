"""
@author: cunyue
@file: backup.py
@time: 2025/6/2 15:07
@description: 日志备份回调
"""

import random

from rich.text import Text

from swanlab.data.callbacker.callback import SwanLabRunCallback
from swanlab.log import swanlog
from swanlab.toolkit import ColumnInfo, MetricInfo, RuntimeInfo
from . import utils as U
from ..namer import generate_colors
from ..run import get_run


class OfflineCallback(SwanLabRunCallback):
    def __init__(self):
        super().__init__()

    def __str__(self) -> str:
        return "SwanLabOfflineCallback"

    def on_init(self, proj_name: str, workspace: str, public: bool = None, logdir: str = None, *args, **kwargs):
        self.run_store.project = proj_name
        self.run_store.workspace = workspace
        self.run_store.visibility = public
        self.run_store.tags = [] if self.run_store.tags is None else self.run_store.tags
        # 设置颜色，随机生成一个
        self.run_store.run_colors = generate_colors(random.randint(0, 20))
        self.transfer.open_for_trace(python_backend='none')

    def on_run(self, *args, **kwargs):
        self._start_terminal_proxy()
        self._register_sys_callback()
        U.print_train_begin(self.run_store.run_dir)
        swanlog.info("Backing up run", Text(self.run_store.run_name, "yellow"), "locally")
        U.print_sync(self.run_store.run_dir)

    def on_runtime_info_update(self, r: RuntimeInfo, *args, **kwargs):
        self.transfer.trace_runtime_info(r)

    def on_column_create(self, column_info: ColumnInfo, *args, **kwargs):
        self.transfer.trace_column(column_info)

    def on_metric_create(self, metric_info: MetricInfo, *args, **kwargs):
        self.transfer.trace_metric(metric_info)

    def on_stop(self, error: str = None, *args, **kwargs):
        U.print_sync(self.run_store.run_dir)
        success = get_run().success
        error_epoch = swanlog.epoch + 1
        self._unregister_sys_callback()
        self.transfer.close_trace(success, error=error, epoch=error_epoch)
