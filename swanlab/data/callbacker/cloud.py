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
from .. import namer as N
from ..run import get_run
from ..run.config import SwanLabConfig
from ..run.metadata.hardware import is_system_key
from ..store import get_run_store, RemoteMetric
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
        http = self._create_client()
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
            try:
                new = http.mount_exp(
                    exp_name=run_store.run_name,
                    colors=run_store.run_colors,
                    description=run_store.description,
                    tags=run_store.tags,
                    cuid=run_store.run_id,
                    must_exist=run_store.resume == 'must',
                )
                run_store.new = new
            except RuntimeError:
                raise RuntimeError(
                    "When resume=must, the experiment must exist in project {}. Please check your parameters.".format(
                        http.proj.name
                    )
                )
            run_store.run_id = http.exp_id
            # 如果不是新实验，需要获取最新的实验配置和指标数据进行本地同步
            if not new:
                # 1. 解析 config，保存到 run_store.config
                config = http.exp.config
                run_store.config = SwanLabConfig.revert_config(config)
                # 2. 获取最新的指标数据，解析为 run_store.metrics
                # 指标总结数据，{log:[{key: '', step: int}] or None, media: [{key: str, step: int}, ...] or None, scalar: [{key: str, step: int}, ...] or None}
                summaries, _ = http.get(
                    "/house/metrics/summaries/{}/{}".format(
                        http.exp.root_proj_cuid or http.proj.cuid, http.exp.root_exp_cuid or http.exp_id
                    ),
                    {"all": True},
                )
                # 设置历史指标条数
                run_store.log_epoch = summaries.get("log", [{}])[0].get("step", 0)
                # 列响应
                columns_resp, _ = http.get("/experiment/{}/column".format(http.exp_id), {"all": True})
                # 列信息, [{error: dict or None, key:str, type:str, class:str}]
                columns = columns_resp.get("list", [])
                # key -> (column_type, column_class, error, latest step)
                metrics: RemoteMetric = {}
                for column in columns:
                    # 从列信息中获取指标信息
                    key = column["key"]
                    column_type = column["type"]
                    column_class = column["class"]
                    error = column.get("error", None)
                    if column_class == "SYSTEM" and not is_system_key(key):
                        # 只记录 sdk 生成的系统指标
                        continue
                    # 从总结数据中获取最新的 step
                    # 这里需要同时查找 media 和 scalar
                    latest_step = None
                    for scalar_summary in summaries.get("scalar") or []:
                        if scalar_summary["key"] == key:
                            latest_step = scalar_summary["step"]
                            break
                    if latest_step is None:
                        for media_summary in summaries.get("media") or []:
                            if media_summary["key"] == key:
                                latest_step = media_summary["step"]
                                break
                    metrics[key] = (column_type, column_class, error, latest_step)
                run_store.metrics = metrics

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
