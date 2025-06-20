"""
@author: cunyue
@file: transfer.py
@time: 2025/6/19 13:58
@description: 在此处定义上传器协议抽象类，这允许上传器不关心如何处理数据细节，只需要返回响应的数据格式即可完成上传
这样可以隔离数据处理与上传逻辑，便于后续向下兼容和维护
"""

import os.path
from abc import ABC, abstractmethod
from typing import Any

from .model import *
from ...data.store import get_run_store


class Transfer(ABC):
    def __init__(self):
        # 媒体目录和文件目录
        run_store = get_run_store()
        self.media_dir = run_store.media_dir
        self.file_dir = run_store.file_dir
        assert os.path.exists(self.media_dir), f"Media directory {self.media_dir} does not exist"
        assert os.path.exists(self.file_dir), f"Files directory {self.file_dir} does not exist"

    @abstractmethod
    def transfer_column(self, data: Any) -> ColumnModel:
        pass

    @abstractmethod
    def transfer_media(self, data: Any) -> MediaModel:
        pass

    @abstractmethod
    def transfer_scalar(self, data: Any) -> ScalarModel:
        """
        转换指标数据为模型
        :param data: 指标数据
        :return: ScalarModel 实例
        """
        pass

    @abstractmethod
    def transfer_file(self, data: Any) -> FileModel:
        pass

    @abstractmethod
    def transfer_log(self, data: Any) -> LogModel:
        pass

    @abstractmethod
    def start(self):
        pass

    @abstractmethod
    def publish_column(self, data: Any):
        """
        发布列任务
        :param data: 列数据列表
        """
        pass

    @abstractmethod
    def publish_scalar(self, data: Any):
        """
        发布指标任务
        :param data: 指标数据列表
        """
        pass

    @abstractmethod
    def publish_media(self, data: Any):
        """
        发布媒体任务
        :param data: 媒体数据列表
        """
        pass

    @abstractmethod
    def publish_file(self, data: Any):
        """
        发布文件任务
        :param data: 文件数据列表
        """
        pass

    @abstractmethod
    def publish_log(self, data: Any):
        """
        发布日志任务
        :param data: 日志数据列表
        """
        pass

    @abstractmethod
    def join(self, error: str = None):
        """
        等待上传线程池完成所有任务
        """
        pass
