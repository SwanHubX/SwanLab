"""
@author: cunyue
@file: backup.py
@time: 2025/6/2 15:07
@description: 日志备份回调
"""

from rich.text import Text

from swanlab.data.run.callback import SwanLabRunCallback
from swanlab.log import swanlog
from swanlab.log.backup import backup
from swanlab.log.type import LogData
from swanlab.toolkit import ColumnInfo, MetricInfo, RuntimeInfo
from ..run import get_run


class OfflineCallback(SwanLabRunCallback):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __str__(self) -> str:
        return "SwanLabOfflineCallback"

    # ---------------------------------- 辅助函数 ----------------------------------
    def _sync_tip_print(self):
        """
        提示用户可以通过命令上传日志到远程服务器
        """
        swanlog.info(
            " ☁️ Run `",
            Text("swanlab sync {}".format(self.fmt_windows_path(self.settings.run_dir))),
            "` to sync logs to remote server",
            sep="",
        )

    @backup("terminal")
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
        self.handle_run()
        self._train_begin_print(self.settings.run_dir)
        swanlog.info("Backing up run", Text(self.settings.exp_name, "yellow"), "locally")
        self._sync_tip_print()

    @backup('runtime')
    def on_runtime_info_update(self, r: RuntimeInfo, *args, **kwargs):
        pass

    @backup("column")
    def on_column_create(self, column_info: ColumnInfo, *args, **kwargs):
        pass

    @backup("metric")
    def on_metric_create(self, metric_info: MetricInfo, *args, **kwargs):
        pass

    def on_stop(self, error: str = None, *args, **kwargs):
        self._sync_tip_print()
        self.backup.stop(error=error, epoch=get_run().swanlog_epoch + 1)
        self._unregister_sys_callback()
