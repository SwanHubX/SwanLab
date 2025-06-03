"""
@author: cunyue
@file: backup.py
@time: 2025/6/2 15:07
@description: 日志备份回调
"""

import os.path
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

from swankit.callback import SwanKitCallback, RuntimeInfo

from swanlab.log import swanlog
from swanlab.log.type import LogData
from swanlab.swanlab_settings import get_settings


def async_io():
    """
    类实例异步 IO 方法装饰器
    BackupCallback 实例携带一个线程池，使用此装饰器可以将被装饰的方法放入线程池中执行
    这样能够避免在主线程中执行 IO 密集型操作，提升性能
    """

    def decorator(func):
        def wrapper(self: "BackupCallback", *args, **kwargs):
            feature = self.executor.submit(func, *args, **kwargs)
            return feature

        return wrapper

    return decorator


class BackupCallback(SwanKitCallback):
    def __init__(self, sidecar: bool = False):
        """
        初始化备份回调
        :param sidecar: 是否为边车模式，此模式下作为其他 callback 的备份回调存在，而不是作为主回调存在，举个例子，边车模式下不会主动注册输出流代理，只会将写入处理器注入
        """
        self.sidecar = sidecar
        # 线程池执行器，每次只会有一个线程在执行，这样的设计原因为：
        # 1. io 操作不会特别耗时
        # 2. 避免多线程写入同一文件导致数据混乱
        self.executor = ThreadPoolExecutor(max_workers=1)
        # 日志保存目录
        self.logdir = None
        # 日志文件写入句柄
        self.f = None

    def __str__(self) -> str:
        return "SwanLabBackupCallback"

    # ---------------------------------- 写入函数 ----------------------------------

    @async_io()
    def _backup_terminal(self, log_data: LogData):
        """
        备份终端输出
        """
        pass

    @async_io()
    def _backup_project(self, proj_name: str, workspace: str, public: bool = False):
        """
        备份项目信息
        :param proj_name: 项目名称
        :param workspace: 工作空间路径
        :param public: 是否公开项目
        """
        pass

    @async_io()
    def _backup_column(self):
        """
        备份指标列信息
        """
        pass

    @async_io()
    def _backup_metric(self):
        """
        备份指标信息
        """
        pass

    @async_io()
    def _backup_media(self):
        """
        备份媒体文件信息
        """
        pass

    @async_io()
    def _backup_runtime_info(self):
        """
        备份实验运行信息
        """
        pass

    @async_io()
    def _backup_error(self):
        """
        备份错误信息
        """
        pass

    # ---------------------------------- 事件回调 ----------------------------------

    def _terminal_handler(self, log_data: LogData):
        """
        终端输出写入操作
        """
        self._backup_terminal(log_data)

    def on_init(self, proj_name: str, workspace: str, public: str = None, logdir: str = None, *args, **kwargs):
        # 1. 创建日志备份目录，保证该目录不存在，如果存在则等一会再创建一个新的目录
        while self.logdir is None or os.path.exists(self.logdir):
            self.logdir is not None and time.sleep(1)
            self.logdir = os.path.join(logdir, "backup-{}".format(datetime.now().strftime("%Y%m%d_%H%M%S")))
        os.mkdir(self.logdir)
        # 2. 创建写入文件句柄
        self.f = open(os.path.join(self.logdir, "backup.swanlab"), 'a', encoding='utf-8')
        # 3. 注册终端输出流代理，如果为边车模式只需要注入处理函数，否则需要完整注册整个日志代理过程
        settings = get_settings()
        if settings.log_proxy_type != "none":
            if self.sidecar:
                swanlog.inject_handler(self._terminal_handler)
            else:
                swanlog.start_proxy(
                    proxy_type=settings.log_proxy_type,
                    max_log_length=settings.max_log_length,
                    handler=self._terminal_handler,
                )
        # 4. 写入日志
        self._backup_project(proj_name, workspace, public=public, logdir=logdir)

    def on_runtime_info_update(self, r: RuntimeInfo, *args, **kwargs):
        pass

    def on_stop(self, error: str = None, *args, **kwargs):
        # 等待线程池中的所有任务完成
        self.executor.shutdown(wait=True)
