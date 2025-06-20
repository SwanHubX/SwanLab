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
from .thread import ThreadPool, UploadType
from .. import get_client
from ...env import get_swanlog_dir


class Transfer(ABC):
    def __init__(self):
        swanlog_dir = get_swanlog_dir()
        # 媒体目录和文件目录
        self.media_dir = os.path.join(swanlog_dir, 'media')
        self.files_dir = os.path.join(swanlog_dir, 'files')
        assert os.path.exists(self.media_dir), f"Media directory {self.media_dir} does not exist"
        assert os.path.exists(self.files_dir), f"Files directory {self.files_dir} does not exist"

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
    def publish_column(self, data: list[Any]):
        """
        发布列任务
        :param data: 列数据列表
        """
        pass

    @abstractmethod
    def publish_scalar(self, data: list[Any]):
        """
        发布指标任务
        :param data: 指标数据列表
        """
        pass

    @abstractmethod
    def publish_media(self, data: list[Any]):
        """
        发布媒体任务
        :param data: 媒体数据列表
        """
        pass

    @abstractmethod
    def publish_file(self, data: list[Any]):
        """
        发布文件任务
        :param data: 文件数据列表
        """
        pass

    @abstractmethod
    def publish_log(self, data: list[Any]):
        """
        发布日志任务
        :param data: 日志数据列表
        """
        pass

    @abstractmethod
    def join(self):
        """
        等待上传线程池完成所有任务
        """
        pass


class PythonTransfer(Transfer, ABC):
    """
    Python 版本上传器，上传服务运行在子线程中
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.media_dir is not None, "Media directory must be specified"
        self._client = get_client(allow_none=True)
        assert self._client is not None, "Client must be initialized before using Transfer"
        # 上传线程池
        self._pool = ThreadPool()

    def start(self):
        pass

    def publish_column(self, data: list[Any]):
        """
        发布列任务
        :param data: 列数据列表
        """
        columns = [self.transfer_column(item) for item in data]
        return self._pool.queue.put((UploadType.COLUMN, columns))

    def publish_scalar(self, data: list[Any]):
        """
        发布指标任务
        :param data: 指标数据列表
        """
        scalars = [self.transfer_scalar(item) for item in data]
        return self._pool.queue.put((UploadType.SCALAR_METRIC, scalars))

    def publish_media(self, data: list[Any]):
        """
        发布媒体任务
        :param data: 媒体数据列表
        """
        medias = [self.transfer_media(item) for item in data]
        self._pool.queue.put((UploadType.MEDIA_METRIC, medias))

    def publish_file(self, data: list[Any]):
        """
        发布文件任务
        :param data: 文件数据列表
        """
        files = [self.transfer_file(item) for item in data]
        self._pool.queue.put((UploadType.FILE, files))

    def publish_log(self, data: list[Any]):
        """
        发布日志任务
        :param data: 日志数据列表
        """
        logs = [self.transfer_log(item) for item in data]
        self._pool.queue.put((UploadType.LOG, logs))

    def join(self):
        """
        等待上传线程池完成所有任务
        """
        self._pool.finish()
