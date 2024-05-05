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


class CloudRunCallback(SwanLabRunCallback):

    def __init__(self, pool: ThreadPool):
        self.pool = pool

    @property
    def cloud(self):
        return self.pool is not None

    def on_train_begin(self):
        pass

    def on_train_end(self):
        pass

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
