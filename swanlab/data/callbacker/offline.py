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
from swanlab.log.type import LogData
from swanlab.toolkit import ColumnInfo, MetricInfo, RuntimeInfo
from . import utils as U
from ..namer import generate_colors
from ..store import get_run_store


class OfflineCallback(SwanLabRunCallback):
    def __init__(self):
        super().__init__()

    def __str__(self) -> str:
        return "SwanLabOfflineCallback"

    def _terminal_handler(self, log_data: LogData):
        """
        终端输出写入操作
        """
        pass

    def on_init(self, proj_name: str, workspace: str, public: bool = None, logdir: str = None, *args, **kwargs):
        # 设置项目缓存
        run_store = get_run_store()
        run_store.project = proj_name
        run_store.workspace = workspace
        run_store.visibility = public
        run_store.tags = [] if run_store.tags is None else run_store.tags
        # 设置颜色，随机生成一个
        run_store.run_colors = generate_colors(random.randint(0, 20))

    def on_run(self, *args, **kwargs):
        super().on_run(*args, **kwargs)
        U.print_train_begin(self.run_store.run_dir)
        swanlog.info("Backing up run", Text(self.run_store.run_name, "yellow"), "locally")
        U.print_sync(self.run_store.run_dir)

    def on_runtime_info_update(self, r: RuntimeInfo, *args, **kwargs):
        # 更新运行时信息
        self.device.write_runtime_info(r, self.run_store.file_dir)

    def on_column_create(self, column_info: ColumnInfo, *args, **kwargs):
        pass

    def on_metric_create(self, metric_info: MetricInfo, *args, **kwargs):
        pass

    def on_stop(self, error: str = None, *args, **kwargs):
        U.print_sync(self.run_store.run_dir)
        super().on_stop(error, *args, **kwargs)
