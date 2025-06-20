"""
@author: cunyue
@file: py.py
@time: 2025/6/20 16:39
@description: python版本上传器基类
"""

from abc import ABC
from typing import Any

from swanlab.core_python import *
from swanlab.core_python.uploader.thread import ThreadPool, UploadType
from swanlab.toolkit import create_time
from ..run.main import get_run


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

    def publish_column(self, data: Any):
        """
        发布列任务
        :param data: 列数据列表
        """
        return self._pool.queue.put((UploadType.COLUMN, self.transfer_column(data)))

    def publish_scalar(self, data: Any):
        """
        发布指标任务
        :param data: 指标数据列表
        """
        return self._pool.queue.put((UploadType.SCALAR_METRIC, self.transfer_scalar(data)))

    def publish_media(self, data: Any):
        """
        发布媒体任务
        :param data: 媒体数据列表
        """
        self._pool.queue.put((UploadType.MEDIA_METRIC, self.transfer_media(data)))

    def publish_file(self, data: Any):
        """
        发布文件任务
        :param data: 文件数据列表
        """
        self._pool.queue.put((UploadType.FILE, self.transfer_file(data)))

    def publish_log(self, data: Any):
        """
        发布日志任务
        :param data: 日志数据列表
        """
        self._pool.queue.put((UploadType.LOG, self.transfer_log(data)))

    def join(self, error: str = None):
        """
        等待上传线程池完成所有任务
        """
        self._pool.finish()
        run = get_run()
        assert run is not None, "run must be initialized"
        assert run.running, "Run must be in running state to join transfer"
        # 上传错误日志
        if error is not None:
            logs = LogModel(
                level="ERROR",
                contents=[{"message": error, "create_time": create_time(), "epoch": run.swanlog_epoch + 1}],
            )
            upload_logs([logs])
