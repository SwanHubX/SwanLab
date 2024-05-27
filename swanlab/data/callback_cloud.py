#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/5 20:22
@File: callback_cloud.py
@IDE: pycharm
@Description:
    云端回调
"""
from .run.callback import MetricInfo, ColumnInfo
from swanlab.cloud import UploadType
from swanlab.error import ApiError
from swanlab.api.upload.model import ColumnModel
from urllib.parse import quote
from swanlab.api import LoginInfo, create_http, terminal_login
from swanlab.api.upload import upload_logs
from swanlab.log import swanlog
from swanlab.utils.font import FONT
from swanlab.api import get_http
from swanlab.utils.key import get_key
from swanlab.utils.judgment import in_jupyter, show_button_html
from swanlab.package import get_host_web, get_host_api
from swanlab.error import KeyFileError
from swanlab.env import get_swanlab_folder
from .callback_local import LocalRunCallback, get_run, SwanLabRunState
from swanlab.cloud import LogSnifferTask, ThreadPool
from swanlab.db import Experiment
from swanlab.utils import create_time
import sys
import os


class CloudRunCallback(LocalRunCallback):
    login_info: LoginInfo = None
    """
    用户登录信息
    """

    def __init__(self):
        super(CloudRunCallback, self).__init__()
        self.pool = ThreadPool()
        self.exiting = False
        """
        标记是否正在退出云端环境
        """

    @classmethod
    def get_login_info(cls):
        """
        发起登录，获取登录信息，执行此方法会覆盖原有的login_info
        """
        key = None
        try:
            key = get_key(os.path.join(get_swanlab_folder(), ".netrc"), get_host_api())[2]
        except KeyFileError:
            fd = sys.stdin.fileno()
            # 不是标准终端，且非jupyter环境，无法控制其回显
            if not os.isatty(fd) and not in_jupyter():
                raise KeyFileError("The key file is not found, call `swanlab.login()` or use `swanlab login` ")
        return terminal_login(key)

    def _view_web_print(self):
        self._watch_tip_print()
        http = get_http()
        project_url = get_host_web() + f"/@{http.groupname}/{http.projname}"
        experiment_url = project_url + f"/runs/{http.exp_id}"
        swanlog.info("🏠 View project at " + FONT.blue(FONT.underline(project_url)))
        swanlog.info("🚀 View run at " + FONT.blue(FONT.underline(experiment_url)))
        return experiment_url

    def _clean_handler(self):
        run = get_run()
        if run is None:
            return swanlog.debug("SwanLab Runtime has been cleaned manually.")
        if self.exiting:
            return swanlog.debug("SwanLab is exiting, please wait.")
        self._train_finish_print()
        # 如果正在运行
        run.finish() if run.is_running else swanlog.debug("Duplicate finish, ignore it.")

    def _except_handler(self, tp, val, tb):
        if self.exiting:
            # FIXME not a good way to fix '\n' problem
            print("")
            swanlog.error("Aborted uploading by user")
            sys.exit(1)
        self._error_print(tp)
        # 结束运行
        get_run().finish(SwanLabRunState.CRASHED, error=self._traceback_error(tb))
        if tp != KeyboardInterrupt:
            raise tp(val)

    def __str__(self):
        return "SwanLabCloudRunCallback"

    def on_init(self, project: str, workspace: str, logdir: str = None) -> int:
        super(CloudRunCallback, self).on_init(project, workspace, logdir)
        if self.login_info is None:
            swanlog.debug("Login info is None, get login info.")
            self.login_info = self.get_login_info()
        http = create_http(self.login_info)
        return http.mount_project(project, workspace).history_exp_count

    def on_run(self):
        swanlog.install(self.settings.console_dir)
        # 注册实验信息
        try:
            get_http().mount_exp(
                exp_name=self.settings.exp_name,
                colors=self.settings.exp_colors,
                description=self.settings.description,
            )
        except ApiError as e:
            if e.resp.status_code == 409:
                FONT.brush("", 50)
                swanlog.error("The experiment name already exists, please change the experiment name")
                Experiment.purely_delete(run_id=self.settings.run_id)
                sys.exit(409)

        # 资源嗅探器
        sniffer = LogSnifferTask(self.settings.files_dir)
        self.pool.create_thread(sniffer.task, name="sniffer", callback=sniffer.callback)

        # 向swanlog注册输出流回调
        def _write_call_call(message):
            self.pool.queue.put((UploadType.LOG, [message]))

        swanlog.set_write_callback(_write_call_call)

        # 注册系统回调
        self._register_sys_callback()
        # 打印信息
        self._train_begin_print()
        swanlog.info("👋 Hi " + FONT.bold(FONT.default(self.login_info.username)) + ", welcome to swanlab!")
        swanlog.info("Syncing run " + FONT.yellow(self.settings.exp_name) + " to the cloud")
        experiment_url = self._view_web_print()

        # 在Jupyter Notebook环境下，显示按钮
        if in_jupyter():
            show_button_html(experiment_url)

    def on_column_create(self, column_info: ColumnInfo):
        self.pool.queue.put((
            UploadType.COLUMN,
            [ColumnModel(column_info.key, column_info.data_type.upper(), column_info.error)]
        ))

    def on_metric_create(self, metric_info: MetricInfo):
        super(CloudRunCallback, self).on_metric_create(metric_info)
        if metric_info.error:
            return
        new_data = metric_info.metric
        new_data['key'] = metric_info.key
        new_data['index'] = metric_info.step
        new_data['epoch'] = metric_info.epoch
        if metric_info.data_type == "default":
            return self.pool.queue.put((UploadType.SCALAR_METRIC, [new_data]))
        key = quote(metric_info.key, safe="")
        data = (new_data, key, metric_info.data_type, metric_info.static_dir)
        self.pool.queue.put((UploadType.MEDIA_METRIC, [data]))

    def on_stop(self, error: str = None):
        # 打印信息
        self._view_web_print()
        run = get_run()
        # 如果正在退出或者run对象为None或者不在云端环境下
        if self.exiting or run is None:
            return swanlog.debug("SwanLab is exiting or run is None, ignore it.")
        state = run.state
        # 标志正在退出（需要在下面的逻辑之前标志）
        self.exiting = True
        sys.excepthook = self._except_handler

        def _():
            # 关闭线程池，等待上传线程完成
            self.pool.finish()
            # 上传错误日志
            if error is not None:
                msg = [{"message": error, "create_time": create_time(), "epoch": swanlog.epoch + 1}]
                upload_logs(msg, level="ERROR")

        FONT.loading("Waiting for uploading complete", _)
        get_http().update_state(state == SwanLabRunState.SUCCESS)
        # 取消注册系统回调
        self._unregister_sys_callback()
        self.exiting = False
