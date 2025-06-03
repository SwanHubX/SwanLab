"""
@author: cunyue
@file: __init__.py
@time: 2025/6/3 18:11
@description: 日志备份处理模块
"""

from concurrent.futures import ThreadPoolExecutor
from typing import Tuple, List

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
            feature = getattr(obj, "executor").submit(func, *args, **kwargs)
            return feature

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
        # 线程池执行器，每次只会有一个线程在执行，这样的设计原因为：
        # 1. io 操作不会特别耗时
        # 2. 避免多线程写入同一文件导致数据混乱
        # 3. 部分用户会将 swanlog 文件夹挂载在 NAS 等对写入并发有限制的存储设备上
        self.executor = ThreadPoolExecutor(max_workers=1) if enable else None
        # 日志文件写入句柄
        self.f = None

        # 动态设置包括项目名在内的一些属性，因为在 on_run 之前句柄还未创建，所以需要先缓存，等执行对应的函数的时候再使用
        self.cache_proj_name = None
        self.cache_workspace = None
        self.cache_public = None

    @enable_check()
    def start(self, run_dir: str):
        """
        开启备份处理器
        """
        pass

    @enable_check()
    def stop(self):
        """
        停止备份处理器
        """
        pass

    @async_io()
    def backup_terminal(self, log_data: LogData):
        """
        备份终端输出
        """
        pass

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

    @async_io()
    def backup_runtime_info(self):
        """
        备份实验运行信息
        """
        pass

    def backup_error(self):
        """
        备份错误信息
        此函数在 on_stop 主线程中调用，因此为同步函数
        """
        pass
