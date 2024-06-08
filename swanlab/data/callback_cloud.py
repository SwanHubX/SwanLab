#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/5 20:22
@File: callback_cloud.py
@IDE: pycharm
@Description:
    云端回调
"""
from .run.callback import MetricInfo, ColumnInfo, RuntimeInfo
from swanlab.data.cloud import UploadType
from swanlab.api.upload.model import ColumnModel, ScalarModel, MediaModel, FileModel
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
from swanlab.data.cloud import ThreadPool
from swanlab.utils import create_time
from swanlab.package import get_package_version, get_package_latest_version
import json
import sys
import os
import io


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
            try:
                fd = sys.stdin.fileno()
                # 不是标准终端，且非jupyter环境，无法控制其回显
                if not os.isatty(fd) and not in_jupyter():
                    raise KeyFileError("The key file is not found, call `swanlab.login()` or use `swanlab login` ")
            # 当使用capsys、capfd或monkeypatch等fixture来捕获或修改标准输入输出时，会抛出io.UnsupportedOperation
            # 这种情况下为用户自定义情况
            except io.UnsupportedOperation:
                pass
        return terminal_login(key)

    @staticmethod
    def _get_package_latest_version():
        """
        cloud模式训练开始时，检测package是否为最新版本
        """
        latest_version = get_package_latest_version()
        local_version = get_package_version()
        if latest_version is not None and latest_version != local_version:
            swanlog.info(f"swanlab version {latest_version} is available!  Upgrade: `pip install -U swanlab`")

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
        run.finish() if run.running else swanlog.debug("Duplicate finish, ignore it.")

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
        # 检测是否有最新的版本
        self._get_package_latest_version()
        if self.login_info is None:
            swanlog.debug("Login info is None, get login info.")
            self.login_info = self.get_login_info()

        http = create_http(self.login_info)
        return http.mount_project(project, workspace).history_exp_count

    def on_run(self):
        swanlog.install(self.settings.console_dir)
        # 注册实验信息
        get_http().mount_exp(
            exp_name=self.settings.exp_name,
            colors=self.settings.exp_colors,
            description=self.settings.description,
        )

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

    def on_runtime_info_update(self, r: RuntimeInfo):
        # 执行local逻辑，保存文件到本地
        super(CloudRunCallback, self).on_runtime_info_update(r)
        # 添加上传任务到线程池
        rc = r.config.to_dict() if r.config is not None else None
        rr = r.requirements.info if r.requirements is not None else None
        rm = r.metadata.to_dict() if r.metadata is not None else None
        # 不需要json序列化，上传时会自动序列化
        f = FileModel(requirements=rr, config=rc, metadata=rm)
        self.pool.queue.put((UploadType.FILE, [f]))

    def on_column_create(self, column_info: ColumnInfo):
        error = None
        if column_info.error is not None:
            error = {"data_class": column_info.error.got, "excepted": column_info.error.expected}
        column = ColumnModel(
            key=column_info.key,
            column_type=column_info.chart.value.column_type,
            error=error
        )
        self.pool.queue.put((UploadType.COLUMN, [column]))

    def on_metric_create(self, metric_info: MetricInfo):
        super(CloudRunCallback, self).on_metric_create(metric_info)
        # 有错误就不上传
        if metric_info.error:
            return
        metric = metric_info.metric
        key = metric_info.column_info.key
        key_encoded = metric_info.key
        step = metric_info.step
        epoch = metric_info.epoch
        # 标量折线图
        if metric_info.column_info.chart == metric_info.column_info.chart.LINE:
            scalar = ScalarModel(metric, key, step, epoch)
            return self.pool.queue.put((UploadType.SCALAR_METRIC, [scalar]))
        # 媒体指标数据

        # -------------------------- 🤡这里是一点小小的💩 --------------------------
        # 要求上传时的文件路径必须带key_encoded前缀
        if metric_info.buffers is not None:
            metric = json.loads(json.dumps(metric))
            for i, d in enumerate(metric["data"]):
                metric["data"][i] = "{}/{}".format(key_encoded, d)
        # ------------------------------------------------------------------------

        media = MediaModel(metric, key, key_encoded, step, epoch, metric_info.buffers)
        self.pool.queue.put((UploadType.MEDIA_METRIC, [media]))

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
