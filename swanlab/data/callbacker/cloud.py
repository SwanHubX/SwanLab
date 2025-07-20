"""
@DATE: 2024/5/5 20:22
@File: callback_cloud.py
@IDE: pycharm
@Description:
    云端回调
"""

import shutil
from concurrent.futures.thread import ThreadPoolExecutor
from typing import Optional

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
from ..porter import Mounter
from ..run import get_run
from ...core_python import *
from ...log.type import LogData


class CloudPyCallback(SwanLabRunCallback):
    login_info: Optional[auth.LoginInfo] = None

    def __init__(self):
        super().__init__()
        self.executor = ThreadPoolExecutor(max_workers=1)

    def __str__(self):
        return "SwanLabCloudPyCallback"

    @staticmethod
    def _create_client():
        try:
            http = get_client()
        except ValueError:
            swanlog.debug("Login info is None, get login info.")
            login_info = CloudPyCallback.login_info
            if login_info is not None:
                # 如果有登录信息，则使用该信息创建客户端
                http = create_client(login_info)
                CloudPyCallback.login_info = None
            else:
                # 如果没有登录信息，则需要用户登录
                # 但是不保存登录信息到本地
                http = create_client(auth.create_login_info(save=False))
        return http

    @staticmethod
    def _converter_summarise_metric():
        pass

    def on_init(self, *args, **kwargs):
        _ = self._create_client()
        # 检测是否有最新的版本
        U.check_latest_version()
        with Status("Creating experiment...", spinner="dots"):
            with Mounter() as mounter:
                mounter.execute()

    def _terminal_handler(self, log_data: LogData):
        self.porter.trace_log(log_data)

    def on_run(self, *args, **kwargs):
        self.porter.open_for_trace(sync=False)
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
        http = get_client()
        if http.pending:
            swanlog.warning("This run was destroyed but it is pending!")
        # 打印信息
        U.print_cloud_web()
        error_epoch = swanlog.epoch + 1
        self._unregister_sys_callback()
        self.porter.close_trace(success, error=error, epoch=error_epoch)
        # 更新实验状态，在此之后实验会话关闭
        http.update_state(success)
        reset_client()
        if not self.user_settings.backup:
            shutil.rmtree(self.run_store.run_dir, ignore_errors=True)
