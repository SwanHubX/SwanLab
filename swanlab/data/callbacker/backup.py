"""
@author: cunyue
@file: backup.py
@time: 2025/6/2 15:07
@description: 日志备份回调
"""

from swankit.callback import RuntimeInfo, ColumnInfo, MetricInfo

from swanlab.data.run.callback import SwanLabRunCallback
from swanlab.log import swanlog
from swanlab.log.type import LogData
from swanlab.swanlab_settings import get_settings


class BackupCallback(SwanLabRunCallback):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __str__(self) -> str:
        return "SwanLabBackupCallback"

    # ---------------------------------- 辅助函数 ----------------------------------

    def _terminal_handler(self, log_data: LogData):
        """
        终端输出写入操作
        """
        pass

    # ---------------------------------- 事件回调 ----------------------------------

    def on_init(self, proj_name: str, workspace: str, public: bool = None, logdir: str = None, *args, **kwargs):
        # 设置项目缓存
        self.backup.cache_proj_name = proj_name
        self.backup.cache_workspace = workspace
        self.backup.cache_public = public

    def on_run(self, *args, **kwargs):
        # 1. 注册终端输出流代理
        settings = get_settings()
        swanlog.start_proxy(
            proxy_type=settings.log_proxy_type,
            max_log_length=settings.max_log_length,
            handler=self._terminal_handler,
        )

    def on_runtime_info_update(self, r: RuntimeInfo, *args, **kwargs):
        pass

    def on_column_create(self, column_info: ColumnInfo, *args, **kwargs):
        pass

    def on_metric_create(self, metric_info: MetricInfo, *args, **kwargs):
        pass

    def on_stop(self, error: str = None, *args, **kwargs):
        pass
