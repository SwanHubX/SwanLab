"""
@author: cunyue
@file: backup_handler.py
@time: 2025/6/4 12:43
@description: 备份处理器，负责出入日志备份写入操作
"""

import os.path
from concurrent.futures import ThreadPoolExecutor
from typing import Tuple, List, Optional, TextIO

from swanlab.log.type import LogData


def enable_check(enable_attr: str = "enable"):
    """装饰器工厂，根据实例属性决定是否执行函数"""

    def decorator(func):
        def wrapper(obj: "BackupHandler", *args, **kwargs):
            if getattr(obj, enable_attr, False):
                return func(obj, *args, **kwargs)
            return None

        return wrapper

    return decorator


def async_io():
    """
    类实例异步 IO 方法装饰器，判断是否需要备份
    BackupHandler 实例携带一个线程池，使用此装饰器可以将被装饰的方法放入线程池中执行
    这样能够避免在主线程中执行 IO 密集型操作，提升性能
    """

    def decorator(func):
        @enable_check()
        def wrapper(obj: "BackupHandler", *args, **kwargs):
            # 与 https://github.com/SwanHubX/SwanLab/issues/889 相同的问题
            # 在回调中线程池已经关闭，我们需要在主线程中执行
            executor: Optional[ThreadPoolExecutor] = getattr(obj, "executor")
            if executor is None or executor._shutdown:
                return func(*args, **kwargs)
            executor.submit(func, *args, **kwargs)

        return wrapper

    return decorator


class BackupHandler:
    """
    备份处理器，负责处理日志备份相关的操作
    """

    def __init__(self, enable: bool = True, backup_type: str = "DEFAULT"):
        super().__init__()
        # 是否启用备份
        self.enable = enable
        self.backup_type = backup_type
        # 线程执行器
        self.executor: Optional[ThreadPoolExecutor] = None
        # 日志文件写入句柄
        self.f: Optional[TextIO] = None

        # 动态设置包括项目名在内的一些属性，因为在 on_run 之前句柄还未创建，所以需要先缓存，等执行对应的函数的时候再使用
        self.cache_proj_name = None
        self.cache_workspace = None
        self.cache_public = None

    @enable_check()
    def start(self, run_dir: str, exp_name: str, colors: Tuple[str, str], description: str, tags: List[str]):
        """
        开启备份处理器，创建日志文件句柄
        此函数的功能包括：
        1. 创建线程执行器
        2. 创建日志文件句柄
        3. 在日志文件头写入当前备份类型和一些元信息
        4. 写入项目、实验备份
        """
        # 创建线程池执行器，每次只会有一个线程在执行，这样的设计原因为：
        # 1. io 操作不会特别耗时
        # 2. 避免多线程写入同一文件导致数据混乱
        # 3. 部分用户会将 swanlog 文件夹挂载在 NAS 等对写入并发有限制的存储设备上
        self.executor = ThreadPoolExecutor(max_workers=1)
        self.f = open(os.path.join(run_dir, "backup.swanlab"), "a", encoding="utf-8")
        # TODO 在日志头写入当前备份类型和一些元信息

    @enable_check()
    def stop(self, error: str = None):
        """
        停止备份处理器
        :param error: str, 如果有错误信息，则在日志中记录
        """
        # 同步停止
        self.executor.shutdown(wait=True)
        # TODO 如果有错误信息则在日志中记录
        # 关闭日志文件句柄
        self.f.close()

    @async_io()
    def backup_terminal(self, log_data: LogData):
        """
        备份终端输出
        """
        self.f.write(str(log_data) + "\n")

    @async_io()
    def backup_proj(self):
        """
        备份项目信息
        """
        pass

    @async_io()
    def backup_exp(self, exp_name: str, colors: Tuple[str, str], description: str, tags: List[str]):
        """
        备份实验信息
        """
        pass

    @async_io()
    def backup_column(self):
        """
        备份指标列信息
        """
        pass

    @async_io()
    def backup_metric(self):
        """
        备份指标信息
        """
        pass

    @async_io()
    def backup_media(self):
        """
        备份媒体文件信息
        """
        pass


def backup(method: str):
    """
    备份装饰器，用于在方法执行前进行备份操作
    """

    def decorator(func):
        def wrapper(obj, *args, **kwargs):
            # 执行备份操作
            backup_obj = getattr(obj, "backup")
            getattr(backup_obj, f"backup_{method}")(backup_obj, *args, **kwargs)
            # 执行原方法
            func(obj, *args, **kwargs)

        return wrapper

    return decorator
