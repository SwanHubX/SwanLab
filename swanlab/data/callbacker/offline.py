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
from .. import namer as N
from ..run import get_run


class OfflineCallback(SwanLabRunCallback):
    def __init__(self):
        super().__init__()

    def __str__(self) -> str:
        return "SwanLabOfflineCallback"

    def on_init(self, proj_name: str, workspace: str, public: bool = None, logdir: str = None, *args, **kwargs):
        run_store = self.run_store
        run_store.project = proj_name
        run_store.workspace = workspace
        run_store.visibility = public
        run_store.tags = [] if run_store.tags is None else run_store.tags
        # 设置颜色，随机生成一个
        exp_count = random.randint(0, 20)
        run_store.run_colors = N.generate_colors(exp_count)
        # 设置名称，随机生成
        run_store.run_name = N.generate_name(exp_count) if run_store.run_name is None else run_store.run_name
        run_store.run_colors = N.generate_colors(exp_count)
        run_store.run_id = N.generate_run_id()
        run_store.new = True

    def on_run(self, *args, **kwargs):
        self.porter.open_for_trace(backend='none', sync=False)
        self._start_terminal_proxy()
        self._register_sys_callback()
        U.print_train_begin(self.run_store.run_dir)
        swanlog.info("Backing up run", Text(self.run_store.run_name, "yellow"), "locally")
        U.print_sync(self.run_store.run_dir)

    def on_runtime_info_update(self, r: RuntimeInfo, *args, **kwargs):
        self.porter.trace_runtime_info(r)

    def on_column_create(self, column_info: ColumnInfo, *args, **kwargs):
        self.porter.trace_column(column_info)

    def on_metric_create(self, metric_info: MetricInfo, *args, **kwargs):
        self.porter.trace_metric(metric_info)

    def on_stop(self, error: str = None, *args, **kwargs):
        U.print_sync(self.run_store.run_dir)
        success = get_run().success
        error_epoch = swanlog.epoch + 1
        self._unregister_sys_callback()
        self.porter.close_trace(success, error=error, epoch=error_epoch)
