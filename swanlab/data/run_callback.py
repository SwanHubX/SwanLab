#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/5 20:22
@File: run_callback.py
@IDE: pycharm
@Description:
    SwanLabRun回调函数
"""
from .run.callback import SwanLabRunCallback, NewKeyInfo
from swanlab.cloud import UploadType
from typing import Optional, Dict
from swanlab.api.upload.model import ColumnModel
from swanlab.cloud import ThreadPool
from urllib.parse import quote
from swanlab.api import LoginInfo
from swanlab.log import swanlog
from swanlab.utils.font import FONT
from swanlab.api import get_http
from swanlab.utils.judgment import in_jupyter, show_button_html
from swanlab.package import get_host_web


class CloudRunCallback(SwanLabRunCallback):

    def __init__(self, pool: ThreadPool, login_info: LoginInfo):
        super(CloudRunCallback, self).__init__()
        self.pool = pool
        self.login_info = login_info

    def _view_web_print(self):
        self._command_tip_print()

        http = get_http()
        project_url = get_host_web() + f"/@{http.groupname}/{http.projname}"
        experiment_url = project_url + f"/runs/{http.exp_id}"
        swanlog.info("🏠 View project at " + FONT.blue(FONT.underline(project_url)))
        swanlog.info("🚀 View run at " + FONT.blue(FONT.underline(experiment_url)))
        return experiment_url

    def on_train_begin(self):
        self._train_begin_print()
        swanlog.info("👋 Hi " + FONT.bold(FONT.default(self.login_info.username)) + ", welcome to swanlab!")
        swanlog.info("Syncing run " + FONT.yellow(self.settings.exp_name) + " to the cloud")
        experiment_url = self._view_web_print()

        # 在Jupyter Notebook环境下，显示按钮
        if in_jupyter():
            show_button_html(experiment_url)

    def on_train_end(self):
        self._view_web_print()

    def on_metric_create(self, key: str, key_info: NewKeyInfo, static_dir: str):
        """
        指标创建回调函数,新增指标信息时调用
        :param key: 指标key名称
        :param key_info: 指标信息
        :param static_dir: 媒体文件目录
        """
        if key_info is None:
            return
        new_data, data_type, step, epoch = key_info
        new_data['key'] = key
        new_data['index'] = step
        new_data['epoch'] = epoch
        if data_type == "default":
            return self.pool.queue.put((UploadType.SCALAR_METRIC, [new_data]))
        key = quote(key, safe="")
        data = (new_data, key, data_type, static_dir)
        self.pool.queue.put((UploadType.MEDIA_METRIC, [data]))

    def on_column_create(self, key, data_type: str, error: Optional[Dict] = None):
        self.pool.queue.put((UploadType.COLUMN, [ColumnModel(key, data_type.upper(), error)]))
