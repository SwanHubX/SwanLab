"""
@author: cunyue
@file: backup.py
@time: 2025/6/2 15:07
@description: 日志备份回调
"""

import random

from rich.text import Text

from swanlab.data.run.callback import SwanLabRunCallback
from swanlab.log import swanlog
from swanlab.log.backup import BackupHandler
from swanlab.log.type import LogData
from swanlab.toolkit import ColumnInfo, MetricInfo, RuntimeInfo
from ..namer import generate_colors
from ..store import get_run_store


class OfflineCallback(SwanLabRunCallback):
    def __init__(self):
        self.device = BackupHandler()
        self.run_store = get_run_store()

    def __str__(self) -> str:
        return "SwanLabOfflineCallback"

    # ---------------------------------- 辅助函数 ----------------------------------
    def _sync_tip_print(self):
        """
        提示用户可以通过命令上传日志到远程服务器
        """
        swanlog.info(
            " ☁️ Run `",
            Text("swanlab sync {}".format(self.fmt_windows_path(self.run_store.run_dir))),
            "` to sync logs to remote server",
            sep="",
        )

    def _terminal_handler(self, log_data: LogData):
        """
        终端输出写入操作
        """
        pass

    # ---------------------------------- 事件回调 ----------------------------------

    def on_init(self, proj_name: str, workspace: str, public: bool = None, logdir: str = None, *args, **kwargs):
        # 设置项目缓存
        run_store = get_run_store()
        run_store.project = proj_name
        run_store.workspace = workspace
        run_store.visibility = public
        run_store.tags = [] if run_store.tags is None else run_store.tags
        # 设置颜色
        run_store.run_colors = generate_colors(random.randint(0, 20))

    def on_run(self, *args, **kwargs):
        self.device.start(
            file_dir=self.run_store.file_dir,
            backup_file=self.run_store.backup_file,
            run_name=self.run_store.run_name,
            workspace=self.run_store.workspace,
            visibility=self.run_store.visibility,
            description=self.run_store.description,
            tags=self.run_store.tags,
        )
        self.handle_run()
        self._train_begin_print(self.run_store.run_dir)
        swanlog.info("Backing up run", Text(self.run_store.run_name, "yellow"), "locally")
        self._sync_tip_print()

    def on_runtime_info_update(self, r: RuntimeInfo, *args, **kwargs):
        # 更新运行时信息
        self.device.write_runtime_info(r, self.run_store.file_dir)

    def on_column_create(self, column_info: ColumnInfo, *args, **kwargs):
        pass

    def on_metric_create(self, metric_info: MetricInfo, *args, **kwargs):
        pass

    def on_stop(self, error: str = None, *args, **kwargs):
        self._sync_tip_print()
        self.device.stop(error=error, epoch=swanlog.epoch + 1)
        self._unregister_sys_callback()
