"""
@author: cunyue
@file: local.py
@time: 2025/2/13 14:59
@description: 本地模式下的回调
local模式（目前）将自动调用swanboard，如果不存在则报错
"""

import random

from rich.text import Text

from swanlab.log.backup import BackupHandler
from swanlab.log.type import LogData
from swanlab.toolkit import ColumnInfo
from ..namer import generate_colors
from ..store import get_run_store
from ...proto.v0 import Column

try:
    # noinspection PyPackageRequirements
    import swanboard
except ImportError:
    raise ImportError("Please install swanboard to use 'local' mode: pip install 'swanlab[dashboard]'")

from importlib.metadata import version

package_version = version("swanboard")
if package_version != "0.1.8b1":
    raise ImportError(
        "Your swanboard version does not match, please use this command to install the matching version: pip install 'swanlab[dashboard]'"
    )


import json
import os
from datetime import datetime
from typing import Tuple, Optional, TextIO

from swanlab.toolkit import RuntimeInfo, MetricInfo
from swanlab.data.run.callback import SwanLabRunCallback
from swanlab.log import swanlog


class LocalRunCallback(SwanLabRunCallback):

    def __init__(self) -> None:
        self.device = BackupHandler()
        self.board = swanboard.SwanBoardCallback()
        self.run_store = get_run_store()
        # 当前日志写入文件的句柄
        self.file: Optional[TextIO] = None

    def __str__(self):
        return "SwanLabLocalRunCallback"

    def _watch_tip_print(self):
        """
        watch命令提示打印
        """
        run_store = get_run_store()
        swanlog.info(
            "🌟 Run `",
            Text("swanlab watch {}".format(self.fmt_windows_path(run_store.swanlog_dir)), "bold"),
            "` to view SwanLab Experiment Dashboard locally",
            sep="",
        )

    def _terminal_handler(self, log_data: LogData):
        log_name = f"{datetime.now().strftime('%Y-%m-%d')}.log"
        run_store = get_run_store()
        if self.file is None:
            # 如果句柄不存在，则创建
            self.file = open(os.path.join(run_store.console_dir, log_name), "a", encoding="utf-8")
        elif os.path.basename(self.file.name) != log_name:
            # 如果句柄存在，但是文件名不一样，则关闭句柄，重新打开
            self.file.close()
            self.file = open(os.path.join(run_store.console_dir, log_name), "a", encoding="utf-8")
        # 写入日志
        for content in log_data["contents"]:
            self.file.write(content['message'] + '\n')
            self.file.flush()

    def on_init(self, proj_name: str, workspace: str, public: bool = None, logdir: str = None, *args, **kwargs):
        self.board.on_init(proj_name)

        run_store = get_run_store()
        run_store.project = proj_name
        run_store.workspace = workspace
        run_store.visibility = public
        run_store.tags = [] if run_store.tags is None else run_store.tags
        # 设置颜色
        run_store.run_colors = generate_colors(random.randint(0, 20))

    def before_run(self, *args, **kwargs):
        run_store = get_run_store()
        # path 是否存在
        if not os.path.exists(run_store.run_dir):
            os.mkdir(run_store.run_dir)

    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        colors: Tuple[str, str],
        *args,
        **kwargs,
    ):
        run_store = get_run_store()
        #  FIXME num 在 dashboard 中被要求传递但是没用上 🤡
        self.board.before_init_experiment(
            os.path.basename(run_store.run_dir), exp_name, description, colors=colors, num=1
        )

    def on_run(self):
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
        run_store = get_run_store()
        # 打印信息
        self._train_begin_print(run_store.run_dir)
        self._watch_tip_print()

    def on_runtime_info_update(self, r: RuntimeInfo, *args, **kwargs):
        # 更新运行时信息
        self.device.write_runtime_info(r, self.run_store.file_dir)

    def on_log(self, *args, **kwargs):
        self.board.on_log(*args, **kwargs)

    def on_column_create(self, column_info: ColumnInfo, *args, **kwargs):
        column = Column.from_column_info(column_info)
        self.device.backup(column)
        # 屏蔽 board 不支持的图表类型和列类型
        if column_info.chart_type.value.chart_type not in ["line", "image", "audio", "text"]:
            return
        if column_info.cls != "CUSTOM":
            return
        self.board.on_column_create(column_info)

    def on_metric_create(self, metric_info: MetricInfo, *args, **kwargs):
        # 对于指标保存，可以随意保存，因为这里与 dashboard 没有直接交互
        # 出现任何错误直接返回
        if metric_info.error:
            return
        # ---------------------------------- 保存指标数据 ----------------------------------
        os.makedirs(os.path.dirname(metric_info.metric_file_path), exist_ok=True)
        os.makedirs(os.path.dirname(metric_info.summary_file_path), exist_ok=True)
        with open(metric_info.summary_file_path, "w+", encoding="utf-8") as f:
            f.write(json.dumps(metric_info.metric_summary, ensure_ascii=False))
        with open(metric_info.metric_file_path, "a", encoding="utf-8") as f:
            f.write(json.dumps(metric_info.metric, ensure_ascii=False) + "\n")
        # ---------------------------------- 保存媒体字节流数据 ----------------------------------
        self.device.write_media_buffer(metric_info)

    def on_stop(self, error: str = None, *args, **kwargs):
        """
        训练结束，取消系统回调
        此函数被`run.finish`调用
        """
        # 写入错误信息
        if error is not None:
            with open(os.path.join(get_run_store().console_dir, "error.log"), "a") as fError:
                print(datetime.now(), file=fError)
                print(error, file=fError)
        self.board.on_stop(error)
        # 打印信息
        self._watch_tip_print()
        self.device.stop(error=error, epoch=swanlog.epoch + 1)
        self._unregister_sys_callback()
