"""
@DATE: 2024/5/5 20:22
@File: callback_cloud.py
@IDE: pycharm
@Description:
    云端回调
"""

import sys
from concurrent.futures.thread import ThreadPoolExecutor

from rich.status import Status
from rich.text import Text

from swanlab.core_python import auth
from swanlab.env import in_jupyter
from swanlab.log import swanlog
from swanlab.log.backup import BackupHandler
from swanlab.toolkit import (
    RuntimeInfo,
    MetricInfo,
    ColumnInfo,
    create_time,
)
from swanlab.transfers import ProtoV0Transfer
from .utils import show_button_html, check_latest_version, view_cloud_web, async_io
from .. import namer as N
from ..run import get_run, SwanLabRunState
from ..run.callback import SwanLabRunCallback
from ..store import get_run_store
from ...core_python import *
from ...log.type import LogData
from ...proto.v0 import Log, Runtime, Column, Metric
from ...swanlab_settings import get_settings


class CloudPyCallback(SwanLabRunCallback):

    def __init__(self):
        run_store = get_run_store()
        self.transfer = ProtoV0Transfer(media_dir=run_store.media_dir, file_dir=run_store.file_dir)
        self.device = BackupHandler()
        self.executor = ThreadPoolExecutor(max_workers=1)
        self.exiting = False
        """
        标记是否正在退出云端环境
        """

    def __str__(self):
        return "SwanLabCloudPyCallback"

    def _clean_handler(self):
        if self.exiting:
            return swanlog.debug("SwanLab is exiting, please wait.")
        super()._clean_handler()

    def _except_handler(self, tp, val, tb):
        if self.exiting:
            print("")
            swanlog.error("Aborted uploading by user")
            sys.exit(1)
        super()._except_handler(tp, val, tb)

    def on_init(self, *args, **kwargs):
        try:
            http = get_client()
        except ValueError:
            swanlog.debug("Login info is None, get login info.")
            http = create_client(auth.create_login_info())
        # 检测是否有最新的版本
        check_latest_version()
        run_store = get_run_store()
        with Status("Creating experiment...", spinner="dots"):
            # 创建实验
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

    def on_run(self, *args, **kwargs):
        self.device.start()
        # 注册运行状态
        self.handle_run()
        # 打印实验开始信息，在 cloud 模式下如果没有开启 backup 的话不打印“数据保存在 xxx”的信息
        swanlab_settings = get_settings()
        run_store = get_run_store()
        self._train_begin_print(save_dir=run_store.run_dir if swanlab_settings.backup else None)
        http = get_client()
        swanlog.info("👋 Hi ", Text(http.username, "bold default"), ",welcome to swanlab!", sep="")
        swanlog.info("Syncing run", Text(run_store.run_name, "yellow"), "to the cloud")
        experiment_url = view_cloud_web()
        # 在Jupyter Notebook环境下，显示按钮
        if in_jupyter():
            show_button_html(experiment_url)

    @async_io()
    def _terminal_handler(self, log_data: LogData):
        logs = Log.from_log_data(log_data)
        for log in logs:
            self.device.backup(log)
            self.transfer.publish_log(log)

    @async_io()
    def on_runtime_info_update(self, r: RuntimeInfo, *args, **kwargs):
        runtime = Runtime.from_runtime_info(r)
        self.device.write_runtime_info(r)
        self.device.backup(runtime)
        self.transfer.publish_file(runtime)

    @async_io()
    def on_column_create(self, column_info: ColumnInfo, *args, **kwargs):
        column = Column.from_column_info(column_info)
        self.device.backup(column)
        self.transfer.publish_column(column)

    @async_io()
    def on_metric_create(self, metric_info: MetricInfo, *args, **kwargs):
        # 有错误就不上传
        if metric_info.error:
            return
        if metric_info.column_info.chart_type == metric_info.column_info.chart_type.LINE:
            # 标量
            scalar = Metric.from_metric_info(metric_info)
            self.device.backup(scalar)
            self.transfer.publish_scalar(scalar)
        else:
            # 媒体
            media = Metric.from_metric_info(metric_info)
            self.device.backup(media)
            self.device.write_media_buffer(metric_info)
            self.transfer.publish_media(media)

    def on_stop(self, error: str = None, *args, **kwargs):
        run = get_run()
        # 如果正在退出或者run对象为None或者不在云端环境下，则不执行任何操作
        # 原因是在云端环境下退出时会新建一个线程完成上传日志等操作，此时回调会重复执行
        # 必须要有个标志表明正在退出
        if self.exiting or run is None:
            return swanlog.debug("SwanLab is exiting or run is None, ignore it.")
        # ---------------------------------- 正在退出 ----------------------------------
        self.exiting = True
        # 打印信息
        view_cloud_web()
        state = run.state
        sys.excepthook = self._except_handler
        error_epoch = swanlog.epoch + 1
        self.device.stop(error=error, epoch=error_epoch)
        # 关闭线程池，等待上传线程完成
        self.transfer.join()
        # 上传错误日志
        if error is not None:
            run = get_run()
            assert run is not None, "run must be initialized"
            assert not run.running, "Run must not be running when joining the transfer pool"
            logs = LogModel(
                level="ERROR",
                contents=[{"message": error, "create_time": create_time(), "epoch": error_epoch}],
            )
            upload_logs([logs])
        get_client().update_state(state == SwanLabRunState.SUCCESS)
        reset_client()
        # 取消注册系统回调
        self._unregister_sys_callback()
        self.exiting = False
        # -------------------------------------------------------------------------
