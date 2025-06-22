"""
@DATE: 2024/5/5 20:22
@File: callback_cloud.py
@IDE: pycharm
@Description:
    云端回调
"""

from concurrent.futures.thread import ThreadPoolExecutor

from rich.status import Status
from rich.text import Text

from swanlab.core_python import auth
from swanlab.data.callbacker.callback import SwanLabRunCallback
from swanlab.env import in_jupyter
from swanlab.log import swanlog
from swanlab.toolkit import (
    RuntimeInfo,
    MetricInfo,
    ColumnInfo,
)
from . import utils as U
from .. import namer as N
from ..run import get_run
from ..store import get_run_store
from ...core_python import *
from ...log.type import LogData


class CloudPyCallback(SwanLabRunCallback):

    def __init__(self):
        super().__init__()
        self.executor = ThreadPoolExecutor(max_workers=1)

    def __str__(self):
        return "SwanLabCloudPyCallback"

    def on_init(self, *args, **kwargs):
        try:
            http = get_client()
        except ValueError:
            swanlog.debug("Login info is None, get login info.")
            http = create_client(auth.create_login_info())
        # 检测是否有最新的版本
        U.check_latest_version()
        run_store = get_run_store()
        with Status("Creating experiment...", spinner="dots"):
            http.mount_project(run_store.project, run_store.workspace, run_store.visibility)
            exp_count = http.history_exp_count
            run_store.run_name = N.generate_name(exp_count) if run_store.run_name is None else run_store.run_name
            run_store.description = "" if run_store.description is None else run_store.description
            run_store.run_colors = N.generate_colors(exp_count)
            run_store.tags = [] if run_store.tags is None else run_store.tags
            http.mount_exp(
                exp_name=run_store.run_name,
                colors=run_store.run_colors,
                description=run_store.description,
                tags=run_store.tags,
            )
        self.porter.open_for_trace(sync=True)

    def _terminal_handler(self, log_data: LogData):
        self.porter.trace_log(log_data)

    def on_run(self, *args, **kwargs):
        # 注册终端代理和系统回调
        self._start_terminal_proxy()
        self._register_sys_callback()
        # 打印实验开始信息，在 cloud 模式下如果没有开启 backup 的话不打印“数据保存在 xxx”的信息
        U.print_train_begin(run_dir=self.run_store.run_dir)
        http = get_client()
        swanlog.info("👋 Hi ", Text(http.username, "bold default"), ",welcome to swanlab!", sep="")
        swanlog.info("Syncing run", Text(self.run_store.run_name, "yellow"), "to the cloud")
        experiment_url = U.print_cloud_web()
        # 在Jupyter Notebook环境下，显示按钮
        if in_jupyter():
            U.show_button_html(experiment_url)

    def on_runtime_info_update(self, r: RuntimeInfo, *args, **kwargs):
        self.porter.trace_runtime_info(r)

    def on_column_create(self, column_info: ColumnInfo, *args, **kwargs):
        self.porter.trace_column(column_info)

    def on_metric_create(self, metric_info: MetricInfo, *args, **kwargs):
        # 有错误就不上传
        if metric_info.error:
            return
        self.porter.trace_metric(metric_info)

    def on_stop(self, error: str = None, *args, **kwargs):
        success = get_run().success
        # 打印信息
        U.print_cloud_web()
        error_epoch = swanlog.epoch + 1
        self._unregister_sys_callback()
        self.porter.close_trace(success, error=error, epoch=error_epoch)
        get_client().update_state(success)
        reset_client()
